#!/usr/bin/env python
"""Assert both mobility quantities are written, read, and in the right units.

IM is 1/K0 in Vs/cm^2 -- Bruker's convention, and what DIA-NN's IM column holds.
CCS is square angstroms. They are NOT interchangeable, and a library that
silently swapped them would still look plausible in a histogram, which is why
this checks physical RANGES rather than only presence.

Usage: check_mobility.py <library.parquet> [roundtrip.parquet]
"""
import sys

import pyarrow.parquet as pq

GAS, T = 28.0134, 305.0
# Tryptic peptides on a timsTOF: 1/K0 lands in ~0.6-1.7, CCS in ~250-950 A^2.
# The bands are deliberately wide -- they are here to catch a unit swap or a
# missing conversion, not to police the model.
IM_RANGE, CCS_RANGE = (0.4, 2.5), (150.0, 1400.0)


def col(t, name):
    return [v for v in t.column(name).to_pylist() if v is not None]


def check(path):
    t = pq.read_table(path, columns=["Precursor.Id", "IM", "CCS",
                                     "Precursor.Mz", "Precursor.Charge"])
    n = t.num_rows
    im, ccs = col(t, "IM"), col(t, "CCS")
    if len(im) != n:
        raise SystemExit(f"{path}: IM null for {n - len(im)}/{n} rows -- the "
                         "generator must emit it and the reader must fill it")
    if len(ccs) != n:
        raise SystemExit(f"{path}: CCS null for {n - len(ccs)}/{n} rows")

    lo, hi = min(im), max(im)
    if not (IM_RANGE[0] <= lo and hi <= IM_RANGE[1]):
        raise SystemExit(f"{path}: IM range [{lo:.3f},{hi:.3f}] is outside 1/K0's "
                         f"{IM_RANGE} -- angstroms written into the IM column?")
    lo, hi = min(ccs), max(ccs)
    if not (CCS_RANGE[0] <= lo and hi <= CCS_RANGE[1]):
        raise SystemExit(f"{path}: CCS range [{lo:.1f},{hi:.1f}] is outside "
                         f"{CCS_RANGE} A^2 -- 1/K0 written into the CCS column?")

    # The two must be consistent with each other, not merely both present.
    mz = t.column("Precursor.Mz").to_pylist()
    z = t.column("Precursor.Charge").to_pylist()
    worst = 0.0
    for a, c, m, q in zip(t.column("IM").to_pylist(), t.column("CCS").to_pylist(), mz, z):
        if a is None or c is None or not q:
            continue
        mi = m * q
        mu = mi * GAS / (mi + GAS)
        want = c * (mu * T) ** 0.5 / (18509.0 * q)
        worst = max(worst, abs(a - want) / max(want, 1e-9))
    if worst > 0.02:
        raise SystemExit(f"{path}: IM and CCS disagree by up to {100*worst:.2f}% "
                         "under Mason-Schamp -- they are not describing the same ion")
    print(f"  {path}: {n} rows, IM and CCS both populated, consistent to "
          f"{100*worst:.4f}%")
    return t


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    first = check(sys.argv[1])
    for other in sys.argv[2:]:
        second = check(other)
        # Row counts are NOT comparable: the compact layout is one row per
        # precursor with list-valued fragment columns, the flat one is a row per
        # transition. Compare per precursor.
        by = {}
        for t in (first, second):
            m = {}
            for pid, a, c in zip(t.column("Precursor.Id").to_pylist(),
                                 t.column("IM").to_pylist(),
                                 t.column("CCS").to_pylist()):
                m.setdefault(pid, (a, c))
            by[id(t)] = m
        a, b = by[id(first)], by[id(second)]
        if set(a) != set(b):
            raise SystemExit(f"round trip changed the precursor set: "
                             f"{len(set(a) ^ set(b))} differ")
        bad = [p for p in a
               if abs((a[p][0] or 0) - (b[p][0] or 0)) > 1e-6
               or abs((a[p][1] or 0) - (b[p][1] or 0)) > 1e-6]
        if bad:
            raise SystemExit(f"round trip changed IM or CCS for {len(bad)} of "
                             f"{len(a)} precursors, e.g. {bad[0]}")
        print(f"  round trip preserved IM and CCS exactly for {len(a)} precursors")
