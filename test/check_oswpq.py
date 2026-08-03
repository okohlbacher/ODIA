#!/usr/bin/env python
"""Check the .oswpq reader against an independent reading of the same bytes.

The independent side is pyarrow plus Python's zipfile -- a different Parquet
implementation and a different ZIP implementation from the ones under test, so
agreement is evidence rather than a tautology.

Totals are not enough. A scatter that puts the right *number* of fragments
under the wrong precursor keeps every count correct, so this compares the full
(precursor -> its fragments) relation row by row.

usage: check_oswpq.py <odia_oswpq_dump> <fixture-dir> [extra bundle ...]
"""
import json
import math
import subprocess
import sys
import zipfile

import pyarrow.parquet as pq
import io

DUMP, FIXTURES = sys.argv[1], sys.argv[2]
EXTRA = sys.argv[3:]
failures = []


def fail(msg):
    failures.append(msg)
    print(f"  FAIL {msg}")


def run(path):
    """Returns (ok, stats dict, [(precursor, [transitions])], stderr)."""
    p = subprocess.run([DUMP, path], capture_output=True, text=True)
    stats, groups = {}, []
    for line in p.stdout.splitlines():
        parts = line.split("\t")
        if parts[0] == "P":
            groups.append(({"seq": parts[1], "prot": parts[2], "mz": float(parts[3]),
                            "charge": int(parts[4]), "decoy": int(parts[5]),
                            "irt": float(parts[6]), "im": float(parts[7]),
                            "count": int(parts[8])}, []))
        elif parts[0] == "T":
            groups[-1][1].append({"mz": float(parts[1]), "intensity": float(parts[2]),
                                  "type": parts[3], "ordinal": int(parts[4]),
                                  "charge": int(parts[5]), "loss": parts[6]})
        elif len(parts) == 2:
            stats[parts[0]] = parts[1]
    return p.returncode == 0, stats, groups, p.stderr


def expected(path):
    """The same relation, read independently."""
    z = zipfile.ZipFile(path)
    prec = pq.read_table(io.BytesIO(z.read("library/precursors.parquet"))).to_pydict()
    trans = pq.read_table(io.BytesIO(z.read("library/transitions.parquet"))).to_pydict()

    # First id wins, mirroring an ambiguous join's only defensible resolution.
    index = {}
    for i, pid in enumerate(prec["precursor_id"]):
        index.setdefault(pid, i)

    groups = [(
        {"seq": prec["modified_sequence"][i],
         "prot": prec["protein_accessions"][i] or "",
         "mz": prec["precursor_mz"][i] or 0.0,
         "charge": prec["charge"][i] or 0,
         "decoy": 1 if prec["decoy"][i] else 0,
         "irt": prec["library_rt"][i],
         "im": prec["library_drift_time"][i]},
        []) for i in range(len(prec["precursor_id"]))]

    orphans = 0
    for r, pid in enumerate(trans["precursor_id"]):
        if pid not in index:
            orphans += 1
            continue
        ann = trans["annotation"][r] or ""
        loss = "noloss"
        if "-" in ann:
            rest = ann.split("-", 1)[1].split("^")[0]
            loss = {"H2O": "H2O", "NH3": "NH3", "H3PO4": "H3PO4", "HPO3": "HPO3",
                    "CO": "CO", "": "noloss"}.get(rest, "?")
        ordinal = trans["ordinal"][r]
        groups[index[pid]][1].append({
            "mz": trans["product_mz"][r] or 0.0,
            "intensity": trans["library_intensity"][r] or 0.0,
            "type": (trans["type"][r] or "?")[:1].lower() or "?",
            "ordinal": 0 if ordinal is None or not 0 <= ordinal <= 255 else ordinal,
            "charge": trans["charge"][r] or 0,
            "loss": loss})
    return groups, orphans


def close(a, b, tol=2e-5):
    return abs(a - b) <= tol


def compare(name, path):
    ok, stats, got, err = run(path)
    if not ok:
        fail(f"{name}: dump failed: {err.strip()}")
        return
    want, orphans = expected(path)

    if len(got) != len(want):
        fail(f"{name}: {len(got)} precursors, expected {len(want)}")
        return
    if int(stats["orphan_transitions"]) != orphans:
        fail(f"{name}: orphan_transitions {stats['orphan_transitions']}, expected {orphans}")

    for i, ((gp, gt), (wp, wt)) in enumerate(zip(got, want)):
        for key in ("seq", "prot", "charge", "decoy"):
            if gp[key] != wp[key]:
                fail(f"{name}: precursor {i} {key} {gp[key]!r}, expected {wp[key]!r}")
        if not close(gp["mz"], wp["mz"]):
            fail(f"{name}: precursor {i} mz {gp['mz']}, expected {wp['mz']}")
        # A null retention time is NaN in the library and prints as 0.0; a
        # library_rt that genuinely is 0.0 is indistinguishable in the dump, so
        # only a present value is compared.
        if wp["irt"] is not None and not close(gp["irt"], wp["irt"], 1e-3):
            fail(f"{name}: precursor {i} irt {gp['irt']}, expected {wp['irt']}")
        # -1 and null both mean absent, and absent must not read as a drift
        # time: a negative 1/K0 matches no frame.
        absent = wp["im"] is None or wp["im"] < 0
        if absent and gp["im"] != -999.0:
            fail(f"{name}: precursor {i} kept an absent drift time as {gp['im']}")
        if not absent and not close(gp["im"], wp["im"], 1e-3):
            fail(f"{name}: precursor {i} im {gp['im']}, expected {wp['im']}")

        if gp["count"] != len(wt) or len(gt) != len(wt):
            fail(f"{name}: precursor {i} ({wp['seq']}) has {len(gt)} fragments "
                 f"(count says {gp['count']}), expected {len(wt)}")
            continue
        # Fragment order within a precursor follows file order, so this catches
        # a scatter that keeps counts right and contents wrong.
        for k, (g, w) in enumerate(zip(gt, wt)):
            for key in ("type", "ordinal", "charge", "loss"):
                if g[key] != w[key]:
                    fail(f"{name}: {wp['seq']} fragment {k} {key} {g[key]!r}, "
                         f"expected {w[key]!r}")
            if not close(g["mz"], w["mz"]):
                fail(f"{name}: {wp['seq']} fragment {k} mz {g['mz']}, expected {w['mz']}")
            if not close(g["intensity"], w["intensity"], 1e-6):
                fail(f"{name}: {wp['seq']} fragment {k} intensity {g['intensity']}, "
                     f"expected {w['intensity']}")
    return stats, got


# --------------------------------------------------------------------------
# Bundles whose full content must survive the round trip unchanged.
for name in ("basic", "f32", "chunked", "shuffled", "orphans", "dupids",
             "decoymismatch", "nulls", "nocensus", "badcensus"):
    compare(name, f"{FIXTURES}/{name}.oswpq")

# The real thing, written by an OpenMS build nobody here controls.
for path in EXTRA:
    compare("upstream", path)

# --------------------------------------------------------------------------
# Reading the same content four ways must give the same library. Chunking,
# intensity width and row order are all representation, not content.
base_ok, base_stats, base_groups, _ = run(f"{FIXTURES}/basic.oswpq")
for variant in ("chunked", "shuffled", "f32"):
    ok, _, groups, err = run(f"{FIXTURES}/{variant}.oswpq")
    if not ok:
        fail(f"{variant}: dump failed: {err.strip()}")
    else:
        # Fragment order within a precursor follows the file, so a reordered
        # file gives a reordered run. What must not change is which fragments
        # belong to which precursor -- compared here as a multiset. f32 is in
        # this list because the fixture's intensities are exact in both widths,
        # so the values must match, not merely round to each other.
        key = lambda fr: tuple(sorted(fr.items()))
        if len(groups) != len(base_groups):
            fail(f"{variant}: {len(groups)} precursors vs basic's {len(base_groups)}")
        for i, ((gp, gt), (bp, bt)) in enumerate(zip(groups, base_groups)):
            if gp != bp:
                fail(f"{variant}: precursor {i} differs from basic: {gp} vs {bp}")
            elif sorted(gt, key=key) != sorted(bt, key=key):
                fail(f"{variant}: precursor {i} ({bp['seq']}) has different fragments "
                     f"from basic: {gt} vs {bt}")

# --------------------------------------------------------------------------
# Specific stats each fixture must produce.
def expect_stat(name, key, value):
    ok, stats, _, err = run(f"{FIXTURES}/{name}.oswpq")
    if not ok:
        fail(f"{name}: dump failed: {err.strip()}")
    elif stats.get(key) != str(value):
        fail(f"{name}: {key} is {stats.get(key)!r}, expected {value!r}")


expect_stat("orphans", "orphan_transitions", 2)
expect_stat("basic", "orphan_transitions", 0)
expect_stat("basic", "childless_precursors", 1)
expect_stat("dupids", "duplicate_precursor_ids", 1)
expect_stat("basic", "duplicate_precursor_ids", 0)
expect_stat("decoymismatch", "decoy_mismatches", 1)
expect_stat("basic", "decoy_mismatches", 0)
expect_stat("basic", "census_present", 1)
expect_stat("basic", "census_agrees", 1)
expect_stat("nocensus", "census_present", 0)
expect_stat("nocensus", "census_agrees", 0)
expect_stat("badcensus", "census_present", 1)
expect_stat("badcensus", "census_agrees", 0)
expect_stat("basic", "decoys", 2)

# --------------------------------------------------------------------------
# Bundles that must be refused, with a reason that names the problem. Loading
# any of these would produce a library rather than an error.
for name, needle in (("future", "schema_version"),
                     ("deflated", "seekable"),
                     ("notalibrary", "no library bundle"),
                     ("notazip", "cannot open archive")):
    p = subprocess.run([DUMP, f"{FIXTURES}/{name}.oswpq"], capture_output=True, text=True)
    if p.returncode == 0:
        fail(f"{name}: loaded, and should have been refused")
    elif needle not in p.stderr:
        fail(f"{name}: refused, but the message does not mention {needle!r}: "
             f"{p.stderr.strip()}")

print(f"oswpq: {len(failures)} failures")
sys.exit(1 if failures else 0)
