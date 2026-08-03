#!/usr/bin/env python
"""Compare ODIA's predicted CCS against DIA-NN's predicted ion mobility.

The two libraries do not carry the same quantity: ODIA predicts collision
cross-section in square angstroms, DIA-NN predicts the 1/K0 a timsTOF reports.
They are related by Mason-Schamp through the drift gas and a calibration
constant C, and C is not known here independently -- it is fitted from this
very data.

That makes one number meaningless and another decisive:

* **The median agreement is zero by construction.** C is chosen to make it so.
  Quoting it would be quoting the fit back at itself.
* **The spread around it is the real result.** It is the combined disagreement
  of two independently trained models, and nothing in the fit can reduce it.

So the offset is shown only to confirm it was removed, and the analysis is the
dispersion: overall, per charge, and against peptide length, where structure
would mean the functional form is wrong rather than the models merely differing.

A split-half check is included because fitting and evaluating on the same rows
cannot detect overfitting even of a single constant.

Usage: compare_ccs_to_diann.py <odia-with-CCS.tsv> <diann.parquet>
"""
import argparse
import csv
import math
import statistics
import sys
import zlib

DRIFT_GAS_MASS = 28.0


def strip_mods(s):
    out, depth = [], 0
    for c in s:
        if c in "([{":
            depth += 1
        elif c in ")]}":
            depth = max(0, depth - 1)
        elif depth == 0 and c.isalpha():
            out.append(c.upper())
    return "".join(out)


def reduced_mass(mz, z, gas=DRIFT_GAS_MASS):
    m_ion = mz * z
    return m_ion * gas / (m_ion + gas)


def implied_coefficient(ccs, mz, z, im, gas=DRIFT_GAS_MASS):
    return ccs * math.sqrt(reduced_mass(mz, z, gas)) / (z * im)


def load_odia(path):
    out = {}
    with open(path, newline="") as f:
        r = csv.reader(f, delimiter="\t")
        h = next(r)
        i = {n: h.index(n) for n in ("Modified.Sequence", "Precursor.Charge", "Decoy",
                                     "Precursor.Mz", "CCS")}
        for row in r:
            if row[i["Decoy"]] != "0" or not row[i["CCS"]]:
                continue
            key = (strip_mods(row[i["Modified.Sequence"]]),
                   int(float(row[i["Precursor.Charge"]])))
            if key not in out:
                out[key] = (float(row[i["CCS"]]), float(row[i["Precursor.Mz"]]))
    return out


def load_diann(path):
    import pyarrow.parquet as pq
    out = {}
    pf = pq.ParquetFile(path)
    for b in pf.iter_batches(batch_size=300000,
                             columns=["Modified.Sequence", "Precursor.Charge",
                                      "Decoy", "IM", "Precursor.Mz"]):
        d = b.to_pydict()
        for k in range(b.num_rows):
            if d["Decoy"][k]:
                continue
            key = (strip_mods(d["Modified.Sequence"][k]), int(d["Precursor.Charge"][k]))
            if key not in out:
                out[key] = (float(d["IM"][k]), float(d["Precursor.Mz"][k]))
    return out


def quantiles(v, qs=(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)):
    s = sorted(v)
    return {q: s[min(len(s) - 1, int(q * len(s)))] for q in qs}


def spearman(pairs):
    def ranks(vals):
        order = sorted(range(len(vals)), key=lambda i: vals[i])
        r = [0.0] * len(vals)
        for pos, i in enumerate(order):
            r[i] = pos
        return r
    xs, ys = ranks([p[0] for p in pairs]), ranks([p[1] for p in pairs])
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    dx = math.sqrt(sum((a - mx) ** 2 for a in xs))
    dy = math.sqrt(sum((b - my) ** 2 for b in ys))
    return num / (dx * dy) if dx and dy else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("odia")
    ap.add_argument("diann")
    args = ap.parse_args()

    sys.stderr.write("loading...\n")
    ours = load_odia(args.odia)
    theirs = load_diann(args.diann)
    keys = sorted(set(ours) & set(theirs))
    sys.stderr.write(f"  {len(ours):,} ODIA, {len(theirs):,} DIA-NN, {len(keys):,} shared\n")

    rows = []
    for key in keys:
        ccs, mz = ours[key]
        im, _ = theirs[key]
        if ccs > 0 and im > 0:
            rows.append((key[0], key[1], ccs, mz, im))
    if len(rows) < 1000:
        raise SystemExit(f"only {len(rows)} usable pairs")

    coeffs = [implied_coefficient(ccs, mz, z, im) for _, z, ccs, mz, im in rows]
    C = statistics.median(coeffs)

    def residual(row, c):
        _, z, ccs, mz, im = row
        predicted = ccs * math.sqrt(reduced_mass(mz, z)) / (z * c)
        return 100.0 * (predicted - im) / im

    res = [residual(r, C) for r in rows]

    print(f"shared precursors        {len(rows):,}")
    print(f"fitted coefficient C     {C:.4f}")
    print()
    print("ODIA's CCS converted to 1/K0, against DIA-NN's 1/K0, relative %:")
    q = quantiles(res)
    print(f"   1%  {q[0.01]:+7.2f}    5%  {q[0.05]:+7.2f}   25%  {q[0.25]:+7.2f}")
    print(f"  50%  {q[0.5]:+7.2f}   75%  {q[0.75]:+7.2f}   95%  {q[0.95]:+7.2f}"
          f"   99%  {q[0.99]:+7.2f}")
    print(f"  median |difference|    {statistics.median([abs(x) for x in res]):.2f}%")

    def within(t):
        return 100.0 * sum(1 for x in res if abs(x) <= t) / len(res)
    print(f"  within 1%: {within(1):.1f}%    within 2%: {within(2):.1f}%"
          f"    within 5%: {within(5):.1f}%")
    print("  (the median is ~0 by construction -- C was fitted to make it so;")
    print("   the spread is the part the fit cannot touch)")
    print()

    # Raw CCS against raw 1/K0 is a weak correlation and it is supposed to be:
    # 1/K0 also depends on charge and mass, so the two are only monotone at
    # fixed z and m/z. The meaningful correlation is after the conversion.
    converted = [(r[2] * math.sqrt(reduced_mass(r[3], r[1])) / (r[1] * C), r[4])
                 for r in rows]
    print(f"Spearman, raw CCS against raw 1/K0:       "
          f"{spearman([(r[2], r[4]) for r in rows]):.4f}  "
          f"(low by construction -- 1/K0 also depends on z and m/z)")
    print(f"Spearman, converted 1/K0 against DIA-NN's: "
          f"{spearman(converted):.4f}")
    within_z = {}
    for r, (c, t) in zip(rows, converted):
        within_z.setdefault(r[1], []).append((c, t))
    for z in sorted(within_z):
        print(f"  at fixed z={z}: {spearman(within_z[z]):.4f}")
    print()

    print("Per precursor charge:")
    for z in sorted({r[1] for r in rows}):
        sub = [r for r in rows if r[1] == z]
        sc = statistics.median([implied_coefficient(c, m, zz, i)
                                for _, zz, c, m, i in sub])
        sr = [abs(residual(r, C)) for r in sub]
        print(f"  z={z}  n={len(sub):>9,}  own C={sc:.2f}  "
              f"median |diff| against the global C = {statistics.median(sr):.2f}%")
    print()

    print("Median |difference| by peptide length. A trend here would mean the")
    print("Mason-Schamp form is not describing the relation, rather than the two")
    print("models merely disagreeing:")
    buckets = {}
    for r in rows:
        buckets.setdefault(min(28, max(8, len(r[0]))) // 4 * 4, []).append(abs(residual(r, C)))
    for b in sorted(buckets):
        v = buckets[b]
        if len(v) > 200:
            print(f"  length {b:>2}-{b + 3:<2}  n={len(v):>9,}  {statistics.median(v):.2f}%")

    print()
    print("Split-half: C fitted on one half, measured on the other.")
    in_a = [i for i, r in enumerate(rows)
            if (zlib.crc32(f"{r[0]}/{r[1]}".encode()) & 1) == 0]
    a = set(in_a)
    other = [i for i in range(len(rows)) if i not in a]
    c_fit = statistics.median([coeffs[i] for i in in_a])
    held = [residual(rows[i], c_fit) for i in other]
    print(f"  C from half A            {c_fit:.4f}  (vs {C:.4f} on all)")
    print(f"  median difference on B   {statistics.median(held):+.3f}%")
    print(f"  median |difference| on B {statistics.median([abs(x) for x in held]):.2f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
