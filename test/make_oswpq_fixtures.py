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

    # A ZIP that is not a library bundle at all.
    with zipfile.ZipFile(f"{outdir}/notalibrary.oswpq", "w") as z:
        z.writestr("runs/runs.parquet", b"not parquet")

    # Not a ZIP at all.
    with open(f"{outdir}/notazip.oswpq", "wb") as f:
        f.write(b"PK\x03\x04" + struct.pack("<I", 0) + b"truncated nonsense")

    print(f"wrote fixtures to {outdir}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
