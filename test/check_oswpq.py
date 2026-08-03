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
        # A null id joins to nothing: two rows that both lack an id are not the
        # same precursor.
        if pid is not None:
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
                    "CO": "CO", "": "noloss"}.get(rest, "other")
        ordinal = trans["ordinal"][r]
        groups[index[pid]][1].append({
            "mz": trans["product_mz"][r] or 0.0,
            "intensity": trans["library_intensity"][r] or 0.0,
            "type": (trans["type"][r] or "?")[:1].lower() or "?",
            "ordinal": 0 if ordinal is None or not 0 <= ordinal <= 255 else ordinal,
            # Stored in one signed byte; outside that it is recorded as
            # unknown rather than truncated to a different charge.
            "charge": (trans["charge"][r] or 0)
                      if -128 <= (trans["charge"][r] or 0) <= 127 else 0,
            "loss": loss})
    return groups, orphans


def close(a, b, tol=6e-6):
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
        # An absent retention time must arrive as NaN, not as 0.0. A null read
        # as zero puts that precursor at the start of the RT axis, and the dump
        # used to print NaN as 0.0, so nothing could tell the two apart.
        if wp["irt"] is None:
            if not math.isnan(gp["irt"]):
                fail(f"{name}: precursor {i} has a null library_rt but read as {gp['irt']}")
        elif math.isnan(gp["irt"]):
            fail(f"{name}: precursor {i} irt is NaN, expected {wp['irt']}")
        elif not close(gp["irt"], wp["irt"], 1e-3):
            fail(f"{name}: precursor {i} irt {gp['irt']}, expected {wp['irt']}")
        # -1 and null both mean absent, and absent must not read as a drift
        # time: a negative 1/K0 matches no frame.
        absent = wp["im"] is None or wp["im"] < 0
        if absent and not math.isnan(gp["im"]):
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
             "decoymismatch", "nulls", "nocensus", "badcensus", "widths",
             "annotations", "nullpid", "edgecharge", "minified", "nometa",
             "census_falsepass", "census_unnested"):
    compare(name, f"{FIXTURES}/{name}.oswpq")

# The real thing, written by an OpenMS build nobody here controls.
#
# Its metadata census -- 7 precursors, 18 transitions, written by OpenMS itself
# -- is the one genuinely independent statement about this format available
# here. Everything else is a fixture this project wrote. Asserting only
# self-consistency wasted it.
for path in EXTRA:
    compare("upstream", path)
    ok, stats, groups, err = run(path)
    if not ok:
        fail(f"upstream: dump failed: {err.strip()}")
    else:
        for key, want in (("precursor_rows", "7"), ("transition_rows", "18"),
                          ("census_present", "1"), ("census_agrees", "1"),
                          ("orphan_transitions", "0"),
                          ("duplicate_precursor_ids", "0"),
                          ("decoy_mismatches", "0"),
                          # Every transition is unannotated: type "", ordinal
                          # -1. A reader built from the format document would
                          # have dropped all 18 as impossible.
                          ("unusable_ordinals", "18"),
                          # One precursor carries m/z 0 and one has no
                          # transitions at all. Both are legal, and both are
                          # things the document does not mention.
                          ("invalid_precursor_mz", "1"),
                          ("childless_precursors", "1")):
            if stats.get(key) != want:
                fail(f"upstream: {key} is {stats.get(key)!r}, expected {want!r}")

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
        # NaN != NaN, so an absent drift time would make every precursor look
        # different from itself.
        def same(x, y):
            if isinstance(x, float) and isinstance(y, float):
                return (math.isnan(x) and math.isnan(y)) or x == y
            return x == y

        def same_dict(x, y):
            return x.keys() == y.keys() and all(same(x[k], y[k]) for k in x)

        key = lambda fr: tuple(sorted((k, repr(v)) for k, v in fr.items()))
        if len(groups) != len(base_groups):
            fail(f"{variant}: {len(groups)} precursors vs basic's {len(base_groups)}")
        for i, ((gp, gt), (bp, bt)) in enumerate(zip(groups, base_groups)):
            if not same_dict(gp, bp):
                fail(f"{variant}: precursor {i} differs from basic: {gp} vs {bp}")
            elif [sorted(x.items()) for x in sorted(gt, key=key)] != \
                 [sorted(x.items()) for x in sorted(bt, key=key)]:
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
# The census counts rows the writer declared against rows read, NOT against
# rows kept. Comparing against kept makes every bundle with an orphan report
# disagreement, and the detector that exists to catch a truncated read before
# an hour of extraction cries wolf instead. This regressed once, silently.
expect_stat("orphans", "census_agrees", 1)
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

# Wide and unusual column types. Each of these was a mutation that survived
# because every fixture used int64, string and bool.
expect_stat("widths", "precursors", 6)
expect_stat("widths", "transitions", 12)
expect_stat("widths", "orphan_transitions", 0)
expect_stat("widths", "decoys", 2)

# Annotation dialects, and the count of the ones this reader cannot parse.
# LossType::Other carries no mass, so an unrecognised loss silently becomes a
# fragment at the parent mass; "y6-Hex" and "b8-H2O-NH3" are the two here.
expect_stat("annotations", "unrecognised_losses", 2)
expect_stat("basic", "unrecognised_losses", 0)

# A null precursor_id must not join to a real precursor whose id is -1.
expect_stat("nullpid", "orphan_transitions", 1)
expect_stat("nullpid", "duplicate_precursor_ids", 0)

# int8 boundaries: 127 and -128 fit, 128 and -129 do not.
expect_stat("edgecharge", "unusable_transition_charges", 2)

# Metadata that parses the same minified as indented, and none at all.
expect_stat("minified", "census_agrees", 1)
expect_stat("nometa", "census_present", 0)
expect_stat("nometa", "precursors", 6)

# The census must be read from the block it names. A "total" in a neighbouring
# block let a bundle declaring 999 transitions and holding 12 report agreement.
expect_stat("census_falsepass", "census_present", 1)
expect_stat("census_falsepass", "census_agrees", 0)
expect_stat("census_unscoped", "census_present", 0)
expect_stat("census_unscoped", "census_agrees", 0)
expect_stat("census_unnested", "census_present", 0)
expect_stat("census_unnested", "census_agrees", 0)

# --------------------------------------------------------------------------
# Bundles that must be refused, with a reason that names the problem. Loading
# any of these would produce a library rather than an error.
for name, needle in (("future", "schema_version"),
                     ("deflated", "seekable"),
                     ("notalibrary", "no library bundle"),
                     ("notazip", "cannot open archive"),
                     # libzip takes the first entry of a repeated name, Python
                     # and unzip take the last, so resolving rather than
                     # refusing gives a library no other tool would read.
                     ("dupentry", "entries named"),
                     # The declared size is the archive's own unverified claim.
                     ("metabomb", "ceiling"),
                     # A short read must fail rather than return a truncated
                     # metadata document, which would still be valid JSON.
                     ("corrupt_meta", "truncated")):
    p = subprocess.run([DUMP, f"{FIXTURES}/{name}.oswpq"], capture_output=True, text=True)
    if p.returncode == 0:
        fail(f"{name}: loaded, and should have been refused")
    elif needle not in p.stderr:
        fail(f"{name}: refused, but the message does not mention {needle!r}: "
             f"{p.stderr.strip()}")

print(f"oswpq: {len(failures)} failures")
sys.exit(1 if failures else 0)
