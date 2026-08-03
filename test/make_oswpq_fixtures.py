#!/usr/bin/env python
"""Build .oswpq bundles that each encode one way the format can go wrong.

The one real bundle available (OpenSwathWorkflow_tworuns_1_17.output.oswpq) is
7 precursors of unannotated, single-chunk, target-only, float64 data. It cannot
reach the cases that actually matter: multi-chunk columns, decoys, float32
intensities, transitions out of precursor order, or a newer schema version.
Every fixture here exists because reading it wrong produces a library rather
than an error.

Written to argv[1]. Requires only pyarrow, which is an independent
implementation of Parquet from the one under test.
"""
import json
import os
import struct
import sys
import zipfile

import pyarrow as pa
import pyarrow.parquet as pq

PRECURSOR_FIELDS = [
    ("precursor_id", pa.int64()),
    ("precursor_mz", pa.float64()),
    ("charge", pa.int32()),
    ("library_rt", pa.float64()),
    ("library_drift_time", pa.float64()),
    ("decoy", pa.bool_()),
    ("traml_id", pa.string()),
    ("modified_sequence", pa.string()),
    ("unmodified_sequence", pa.string()),
    ("protein_accessions", pa.string()),
]


def transition_fields(intensity_type=pa.float64()):
    return [
        ("transition_id", pa.int64()),
        ("precursor_id", pa.int64()),
        ("traml_id", pa.string()),
        ("product_mz", pa.float64()),
        ("charge", pa.int32()),
        ("type", pa.string()),
        ("annotation", pa.string()),
        ("ordinal", pa.int32()),
        ("detecting", pa.bool_()),
        ("identifying", pa.bool_()),
        ("quantifying", pa.bool_()),
        ("library_intensity", intensity_type),
        ("decoy", pa.bool_()),
    ]


def table(fields, rows):
    schema = pa.schema([pa.field(n, t) for n, t in fields])
    cols = {n: [r.get(n) for r in rows] for n, _ in fields}
    return pa.table([pa.array(cols[n], type=t) for n, t in fields], schema=schema)


def bundle(path, precursors, transitions, *, intensity_type=pa.float64(),
           row_group_size=None, schema_version=1, compress=zipfile.ZIP_STORED,
           census=True, omit_metadata=False):
    ptab = table(PRECURSOR_FIELDS, precursors)
    ttab = table(transition_fields(intensity_type), transitions)

    def parquet_bytes(tab):
        sink = pa.BufferOutputStream()
        kwargs = {}
        if row_group_size:
            kwargs["row_group_size"] = row_group_size
        pq.write_table(tab, sink, **kwargs)
        return sink.getvalue().to_pybytes()

    meta = {"openms": {
        "schema_version": schema_version,
        "generator": "make_oswpq_fixtures",
        "openms_version": "fixture",
    }}
    if census:
        meta["openms"]["counts"] = {
            "proteins": {"total": 1, "target": 1, "decoy": 0},
            "peptides": {"total": len(precursors), "target": len(precursors), "decoy": 0},
            "precursors": {"total": len(precursors),
                           "target": sum(1 for p in precursors if not p["decoy"]),
                           "decoy": sum(1 for p in precursors if p["decoy"])},
            "compounds": {"total": len(precursors), "target": len(precursors), "decoy": 0},
            "transitions": {"total": len(transitions),
                            "target": sum(1 for t in transitions if not t["decoy"]),
                            "decoy": sum(1 for t in transitions if t["decoy"])},
        }

    with zipfile.ZipFile(path, "w", compression=compress) as z:
        z.writestr("library/metadata.json", json.dumps(meta, indent=2))
        z.writestr("library/precursors.parquet", parquet_bytes(ptab))
        z.writestr("library/transitions.parquet", parquet_bytes(ttab))
        if omit_metadata:
            pass
    return path


def precursor(pid, seq, mz, charge=2, rt=10.0, im=-1.0, decoy=False, prot="ProteinA"):
    return {"precursor_id": pid, "precursor_mz": mz, "charge": charge,
            "library_rt": rt, "library_drift_time": im, "decoy": decoy,
            "traml_id": f"{seq}_{charge}", "modified_sequence": seq,
            "unmodified_sequence": seq, "protein_accessions": prot}


def transition(tid, pid, mz, intensity, type_="y", ordinal=3, charge=1,
               annotation=None, decoy=False):
    if annotation is None:
        annotation = f"{type_}{ordinal}"
    return {"transition_id": tid, "precursor_id": pid, "traml_id": str(tid),
            "product_mz": mz, "charge": charge, "type": type_,
            "annotation": annotation, "ordinal": ordinal, "detecting": True,
            "identifying": False, "quantifying": True,
            "library_intensity": intensity, "decoy": decoy}


def basic_content():
    """Six precursors, three of them decoys, with annotated fragments.

    Fragment m/z ascends with transition_id so that a scatter putting fragments
    under the wrong precursor is visible in the dump rather than only in a
    count. Losses and charges vary, because a reader that ignores the
    annotation entirely would otherwise agree on every row.
    """
    precursors = [
        precursor(10, "PEPTIDEK", 450.25, rt=10.5, im=0.85),
        precursor(20, "SAMPLER", 400.75, rt=-5.25, im=1.10),
        precursor(30, "ELVISLIVESK", 620.35, charge=3, rt=88.0),
        precursor(40, "KEDITPEP", 450.25, rt=10.5, im=0.85, decoy=True),
        precursor(50, "RELPMAS", 400.75, rt=-5.25, im=1.10, decoy=True),
        precursor(60, "EMPTYPEPK", 500.0, rt=1.0),   # deliberately childless
    ]
    transitions = [
        transition(1, 10, 300.10, 1.0, "y", 3),
        transition(2, 10, 400.20, 0.5, "b", 4, annotation="b4-H2O"),
        transition(3, 10, 500.30, 0.25, "y", 5, charge=2, annotation="y5^2"),
        transition(4, 20, 310.10, 1.0, "y", 2),
        transition(5, 20, 410.20, 0.75, "b", 6, annotation="b6-NH3"),
        transition(6, 30, 320.10, 1.0, "y", 7),
        transition(7, 30, 420.20, 0.60, "y", 8, annotation="y8-H3PO4"),
        transition(8, 30, 520.30, 0.30, "b", 2, annotation="b2-HPO3"),
        transition(9, 40, 330.10, 1.0, "y", 3, decoy=True),
        transition(10, 40, 430.20, 0.50, "b", 4, decoy=True),
        transition(11, 50, 340.10, 1.0, "y", 2, decoy=True),
        transition(12, 50, 440.20, 0.75, "b", 6, decoy=True),
    ]
    return precursors, transitions


def widths_content():
    """The same library in the widest types a writer might legitimately choose.

    Every one of these was a surviving mutation: dropping LARGE_STRING from
    getString, reading UINT32 as int32, and removing getBool's integer
    fallback all passed a suite whose only fixtures were int64/string/bool.
    Dictionary encoding is here for a second reason -- it is the only way to
    make Arrow return more than one chunk per column at test scale, because it
    is not concatenated across row groups the way a plain string column is.
    """
    precursors, transitions = basic_content()
    # Ids beyond int32, to catch a narrowing cast in the join.
    for i, p in enumerate(precursors):
        p["precursor_id"] = 3_000_000_000 + i * 7
    for t in transitions:
        t["precursor_id"] = 3_000_000_000 + {10: 0, 20: 1, 30: 2, 40: 3, 50: 4}[
            t["precursor_id"]] * 7
    return precursors, transitions


def main(outdir):
    os.makedirs(outdir, exist_ok=True)
    p, t = basic_content()

    bundle(f"{outdir}/basic.oswpq", p, t)

    # float32 intensities: the OpenDIAlyzer patch's breaking on-disk change. A
    # reader that pins float64 rejects every file a patched build writes.
    bundle(f"{outdir}/f32.oswpq", p, t, intensity_type=pa.float32())

    # One row group per row, so every column comes back as many chunks. This is
    # the shape a >2 GB text column is forced into at proteome scale, where
    # reading chunk(0) alone walks off the end of the column.
    bundle(f"{outdir}/chunked.oswpq", p, t, row_group_size=1)

    # Transitions in an order unrelated to the precursor table. Nothing in the
    # format promises grouping, and fragments under the wrong precursor is a
    # worse outcome than a file that will not load.
    shuffled = [t[i] for i in (7, 2, 11, 0, 5, 9, 3, 10, 1, 6, 4, 8)]
    bundle(f"{outdir}/shuffled.oswpq", p, shuffled)

    # Transitions pointing at precursor ids that do not exist. They must be
    # dropped and counted, not attached to whatever index the id happens to
    # collide with.
    orphaned = t + [transition(90, 999, 600.0, 1.0), transition(91, -1, 601.0, 1.0)]
    bundle(f"{outdir}/orphans.oswpq", p, orphaned)

    # A repeated precursor_id makes the join ambiguous. Silence here is
    # indistinguishable from a precursor that genuinely has no fragments.
    dup = p + [precursor(10, "IMPOSTORK", 999.0)]
    bundle(f"{outdir}/dupids.oswpq", dup, t)

    # A transition claiming to be a target under a decoy precursor. The column
    # is redundant, so disagreement means the file is wrong.
    mismatched = list(t)
    mismatched[8] = dict(mismatched[8], decoy=False)
    bundle(f"{outdir}/decoymismatch.oswpq", p, mismatched)

    # Nulls, which the format allows and which are not zero: an absent drift
    # time and a drift time of 0.0 are different claims about the data.
    nulled = [dict(x) for x in p]
    nulled[0]["library_drift_time"] = None
    nulled[1]["library_rt"] = None
    nulled[2]["decoy"] = None
    nulled_t = [dict(x) for x in t]
    nulled_t[0]["library_intensity"] = None
    nulled_t[1]["ordinal"] = None
    bundle(f"{outdir}/nulls.oswpq", nulled, nulled_t)

    # A newer schema version must be refused. The field exists so that the next
    # breaking change is an error rather than a plausible library.
    bundle(f"{outdir}/future.oswpq", p, t, schema_version=2)

    # No census: the reader must load, and must not claim agreement it cannot
    # check.
    bundle(f"{outdir}/nocensus.oswpq", p, t, census=False)

    # A census that disagrees with the rows. This is the cheapest detector for
    # a truncated read there is, and it must fire.
    lying = f"{outdir}/badcensus.oswpq"
    bundle(lying, p, t)
    with zipfile.ZipFile(lying) as z:
        items = {n: z.read(n) for n in z.namelist()}
    meta = json.loads(items["library/metadata.json"])
    meta["openms"]["counts"]["transitions"]["total"] = 999
    items["library/metadata.json"] = json.dumps(meta, indent=2).encode()
    with zipfile.ZipFile(lying, "w", zipfile.ZIP_STORED) as z:
        for n, b in items.items():
            z.writestr(n, b)

    # DEFLATEd entries. Parquet cannot be read without random access, so this
    # must fail loudly rather than return a short table.
    bundle(f"{outdir}/deflated.oswpq", p, t, compress=zipfile.ZIP_DEFLATED)

    # Wide and unusual column types, written with several row groups so the
    # dictionary columns arrive as several Arrow chunks.
    wp, wt = widths_content()
    wide_precursor_fields = [
        ("precursor_id", pa.uint64()), ("precursor_mz", pa.float32()),
        ("charge", pa.int8()), ("library_rt", pa.float32()),
        ("library_drift_time", pa.float64()), ("decoy", pa.int64()),
        ("traml_id", pa.large_string()), ("modified_sequence", pa.large_string()),
        ("unmodified_sequence", pa.string()),
        ("protein_accessions", pa.dictionary(pa.int32(), pa.string())),
    ]
    wide_transition_fields = [
        ("transition_id", pa.uint32()), ("precursor_id", pa.uint64()),
        ("traml_id", pa.string()), ("product_mz", pa.float64()),
        ("charge", pa.int16()), ("type", pa.dictionary(pa.int32(), pa.string())),
        ("annotation", pa.dictionary(pa.int32(), pa.string())),
        ("ordinal", pa.int64()), ("detecting", pa.bool_()),
        ("identifying", pa.bool_()), ("quantifying", pa.bool_()),
        ("library_intensity", pa.float32()), ("decoy", pa.uint8()),
    ]
    # Not `for p in wp` -- p and t are the fixture data for every later bundle,
    # and rebinding them here made the next bundle() receive a single dict.
    for row in wp:
        row["decoy"] = 1 if row["decoy"] else 0
    for row in wt:
        row["decoy"] = 1 if row["decoy"] else 0

    def wide_bundle(path):
        ptab = table(wide_precursor_fields, wp)
        ttab = table(wide_transition_fields, wt)
        meta = {"openms": {"schema_version": 1, "generator": "make_oswpq_fixtures",
                           "counts": {
                               "precursors": {"total": len(wp)},
                               "transitions": {"total": len(wt)}}}}

        def pbytes(tab):
            sink = pa.BufferOutputStream()
            pq.write_table(tab, sink, row_group_size=1)
            return sink.getvalue().to_pybytes()

        with zipfile.ZipFile(path, "w", zipfile.ZIP_STORED) as z:
            z.writestr("library/metadata.json", json.dumps(meta))
            z.writestr("library/precursors.parquet", pbytes(ptab))
            z.writestr("library/transitions.parquet", pbytes(ttab))

    wide_bundle(f"{outdir}/widths.oswpq")

    # Annotation dialects. Every one of these appears in real libraries, and
    # the two lines that handle them -- stripping the charge suffix, and taking
    # the loss from the FIRST hyphen -- were both deletable without failing a
    # test, because no fixture combined a loss with a charge.
    ap, at_ = basic_content()
    dialects = [
        ("y7-H2O", "y", 7, 1), ("b3^2", "b", 3, 2), ("y10-NH3^2", "y", 10, 2),
        ("y7(2+)-H2O", "y", 7, 2), ("y2^2-H3PO4", "y", 2, 2),
        ("b4-HPO3", "b", 4, 1), ("y5-CO", "y", 5, 1),
        # Unrecognised, and it must be counted: LossType::Other carries no mass,
        # so such a fragment silently becomes one at the parent mass.
        ("y6-Hex", "y", 6, 1), ("b8-H2O-NH3", "b", 8, 1),
        ("precursor", "p", 0, 1), ("y3-", "y", 3, 1),
    ]
    at_ = [transition(500 + i, 10, 300.0 + i, 1.0 - 0.05 * i, ty, n, z, annotation=a)
           for i, (a, ty, n, z) in enumerate(dialects)]
    bundle(f"{outdir}/annotations.oswpq", ap[:1], at_)

    # A null precursor_id must not collide with a real precursor whose id is
    # -1. Before the sentinel moved out of the value domain, the transition
    # below was silently attached to REALMINUSONE.
    np_ = [precursor(-1, "REALMINUSONE", 500.0), precursor(7, "NULLID", 600.0)]
    np_[1]["precursor_id"] = None
    nt_ = [transition(1, -1, 300.0, 1.0), transition(2, None, 400.0, 1.0)]
    bundle(f"{outdir}/nullpid.oswpq", np_, nt_)

    # Charges at the edge of the stored width. An off-by-one in the range check
    # is invisible unless a fixture sits on the boundary.
    cp, _ = basic_content()
    ct = [transition(1, 10, 300.0, 1.0, charge=127),
          transition(2, 10, 301.0, 1.0, charge=128),
          transition(3, 10, 302.0, 1.0, charge=-128),
          transition(4, 10, 303.0, 1.0, charge=-129)]
    bundle(f"{outdir}/edgecharge.oswpq", cp[:1], ct)

    # A census whose "total" lives in a block the reader must not wander into.
    # The substring scan this replaced took the first "total" after the first
    # "transitions" at any depth, so it read a neighbouring block's count and
    # reported that a truncated bundle agreed with its census.
    liar = f"{outdir}/census_falsepass.oswpq"
    bundle(liar, p, t)
    with zipfile.ZipFile(liar) as z:
        items = {n: z.read(n) for n in z.namelist()}
    items["library/metadata.json"] = json.dumps({"openms": {
        "schema_version": 1,
        "generator": "make_oswpq_fixtures",
        "fragment_type_counts": {"precursors": {"total": 6}, "transitions": {"total": 12}},
        "counts": {
            "precursors": {"target": 4, "decoy": 2, "total": 6},
            "transitions": {"target": 8, "decoy": 4, "total": 999},
        }}}).encode()
    with zipfile.ZipFile(liar, "w", zipfile.ZIP_STORED) as z:
        for n, b in items.items():
            z.writestr(n, b)

    # Minified metadata, which must parse exactly as the indented form does.
    mini = f"{outdir}/minified.oswpq"
    bundle(mini, p, t)
    with zipfile.ZipFile(mini) as z:
        items = {n: z.read(n) for n in z.namelist()}
    items["library/metadata.json"] = json.dumps(
        json.loads(items["library/metadata.json"]), separators=(",", ":")).encode()
    with zipfile.ZipFile(mini, "w", zipfile.ZIP_STORED) as z:
        for n, b in items.items():
            z.writestr(n, b)

    # No metadata entry at all: legal, and the false branch of zip.has(METADATA)
    # had no fixture.
    nometa = f"{outdir}/nometa.oswpq"
    bundle(nometa, p, t)
    with zipfile.ZipFile(nometa) as z:
        items = {n: z.read(n) for n in z.namelist() if n != "library/metadata.json"}
    with zipfile.ZipFile(nometa, "w", zipfile.ZIP_STORED) as z:
        for n, b in items.items():
            z.writestr(n, b)

    # Two entries with the same name. libzip takes the first, Python's zipfile
    # and unzip take the last, so a reader that resolves rather than refuses
    # gives a different library from every other tool.
    dupentry = f"{outdir}/dupentry.oswpq"
    bundle(dupentry, p, t)
    with zipfile.ZipFile(dupentry) as z:
        items = [(n, z.read(n)) for n in z.namelist()]
    other, _ = basic_content()
    other[0]["modified_sequence"] = "IMPOSTOR"
    sink = pa.BufferOutputStream()
    pq.write_table(table(PRECURSOR_FIELDS, other), sink)
    with zipfile.ZipFile(dupentry, "w", zipfile.ZIP_STORED) as z:
        for n, b in items:
            z.writestr(n, b)
        z.writestr("library/precursors.parquet", sink.getvalue().to_pybytes())

    # A metadata entry declaring far more than it holds. The declared size is
    # the archive's own unverified claim, and allocating on it made a 4 MB
    # bundle take 4.25 GB.
    bomb = f"{outdir}/metabomb.oswpq"
    bundle(bomb, p, t)
    with zipfile.ZipFile(bomb) as z:
        items = {n: z.read(n) for n in z.namelist()}
    items["library/metadata.json"] = json.dumps({"openms": {"schema_version": 1}}).encode() + \
        b" " * (48 << 20)
    with zipfile.ZipFile(bomb, "w", zipfile.ZIP_DEFLATED) as z:
        for n, b in items.items():
            z.writestr(n, b, zipfile.ZIP_STORED if n.endswith(".parquet") else zipfile.ZIP_DEFLATED)

    # A census whose named block has no "total". Scoped lookup finds nothing
    # and must say the census is absent; a scan that wanders on finds the 12 in
    # the block after it and reports a census that was never declared.
    unscoped = f"{outdir}/census_unscoped.oswpq"
    bundle(unscoped, p, t)
    with zipfile.ZipFile(unscoped) as z:
        items = {n: z.read(n) for n in z.namelist()}
    items["library/metadata.json"] = json.dumps({"openms": {
        "schema_version": 1,
        "counts": {"precursors": {"target": 4, "decoy": 2},
                   "transitions": {"target": 8, "decoy": 4}},
        "elsewhere": {"precursors": {"total": 6}, "transitions": {"total": 12}},
    }}).encode()
    with zipfile.ZipFile(unscoped, "w", zipfile.ZIP_STORED) as z:
        for n, b in items.items():
            z.writestr(n, b)

    # A census with no "counts" block at all, but blocks of the right names one
    # level up. A path lookup that skips a missing element instead of failing
    # walks into them and reports a census the file never declared.
    unnested = f"{outdir}/census_unnested.oswpq"
    bundle(unnested, p, t)
    with zipfile.ZipFile(unnested) as z:
        items = {n: z.read(n) for n in z.namelist()}
    items["library/metadata.json"] = json.dumps({"openms": {
        "schema_version": 1,
        "precursors": {"total": 6},
        "transitions": {"total": 12},
    }}).encode()
    with zipfile.ZipFile(unnested, "w", zipfile.ZIP_STORED) as z:
        for n, b in items.items():
            z.writestr(n, b)

    # A metadata entry whose compressed stream is corrupt.
    #
    # Not an overstated size: within a well-formed archive that reads the
    # neighbouring entry's bytes rather than hitting EOF, so the short-read
    # branch is unreachable that way. A broken deflate stream is the realistic
    # corruption, and it must fail rather than yield a truncated document --
    # a JSON census cut short is still valid JSON, with a smaller count.
    corrupt = f"{outdir}/corrupt_meta.oswpq"
    bundle(corrupt, p, t)
    with zipfile.ZipFile(corrupt) as z:
        items = {n: z.read(n) for n in z.namelist()}
    with zipfile.ZipFile(corrupt, "w", zipfile.ZIP_STORED) as z:
        for n, b in items.items():
            if n.endswith("metadata.json"):
                z.writestr(zipfile.ZipInfo(n), b, compress_type=zipfile.ZIP_DEFLATED)
            else:
                z.writestr(n, b)
    raw = bytearray(open(corrupt, "rb").read())
    at = raw.find(b"library/metadata.json")
    # Past the local header's name field is the deflate stream; flipping bytes
    # in the middle of it leaves the sizes intact and the data unrecoverable.
    body = at + len("library/metadata.json")
    for i in range(body + 8, body + 40):
        raw[i] ^= 0xFF
    open(corrupt, "wb").write(bytes(raw))

    # A ZIP that is not a library bundle at all.
    with zipfile.ZipFile(f"{outdir}/notalibrary.oswpq", "w") as z:
        z.writestr("runs/runs.parquet", b"not parquet")

    # Not a ZIP at all.
    with open(f"{outdir}/notazip.oswpq", "wb") as f:
        f.write(b"PK\x03\x04" + struct.pack("<I", 0) + b"truncated nonsense")

    print(f"wrote fixtures to {outdir}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
