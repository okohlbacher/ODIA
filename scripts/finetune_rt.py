#!/usr/bin/env python
"""Fine-tune AlphaPeptDeep's retention-time model on one run's own identifications.

Why this exists: the stock model's output saturates. Its top two deciles are
compressed 5x and 12x, 12.8% of a human library lands in one output bin, and no
calibration can undo that -- a monotone map can straighten a bent axis, not a
collapsed one. Measured on S08, held-out residual falls 0.701 -> 0.449 min from
500 training peptides in 31 seconds, and afterwards a flexible calibration buys
nothing over a straight line, which is the evidence that the defect was in the
model rather than the mapping. See doc/06-rt-refinement-plan.md.

This runs OUTSIDE ODIA, in a environment carrying peptdeep and torch, and its
output is an ONNX file ODIA reads with `-rt_model`. torch never enters ODIA's
runtime.

**Modifications are resolved through alphabase's own table**, all 2,788 entries
carrying a UniMod id, rather than a hand-written dictionary. The prototype this
replaces knew five modifications and raised on the sixth; that is a fine failure
mode and a poor limit.

**The retention times must come from a search of the run being tuned for.** That
is the point -- the model learns this chromatography -- and it is also the
hazard: if those identifications came from searching with the library this model
will go on to build, the loop is closed and the FDR of the second search is no
longer independent of the first. Feeding this from an external search is safe;
feeding it from ODIA's own first pass is the case that needs the entrapment
measurement doc/06 describes. This script does not know which it has been given,
so it records the source in the sidecar and leaves the judgement to the caller.
"""
import argparse
import json
import hashlib
import os
import re
import sys
import time

MOD_TOKEN = re.compile(r"\(UniMod:(\d+)\)|\[UniMod:(\d+)\]", re.IGNORECASE)


def build_unimod_index():
    """{unimod id: [alphabase mod names]}, from alphabase's own table."""
    from alphabase.constants.modification import MOD_DF
    index = {}
    for name, uid in zip(MOD_DF["mod_name"], MOD_DF["unimod_id"]):
        if uid and uid > 0:
            index.setdefault(int(uid), []).append(name)
    return index


def parse_modified_sequence(seq, index):
    """DIA-NN's Modified.Sequence -> (bare sequence, mods, mod_sites).

    AlphaBase names a modification by residue -- Carbamidomethyl@C -- and by
    terminus for N-terminal ones, so the UniMod id alone is ambiguous and the
    residue it sits on has to pick between the candidates.
    """
    bare, mods, sites = [], [], []
    pos = 0
    for m in MOD_TOKEN.finditer(seq):
        bare.append(seq[pos:m.start()])
        uid = int(m.group(1) or m.group(2))
        site = len("".join(bare))
        residue = bare[-1][-1] if bare and bare[-1] else ""
        candidates = index.get(uid, [])
        chosen = None
        if site == 0:
            # Before any residue: an N-terminal modification.
            for c in candidates:
                if "N-term" in c:
                    chosen = c
                    break
        if chosen is None:
            for c in candidates:
                if c.endswith("@" + residue):
                    chosen = c
                    break
        if chosen is None:
            for c in candidates:
                if "N-term" in c or "C-term" in c:
                    chosen = c
                    break
        if chosen is None:
            raise SystemExit(
                f"UniMod:{uid} on residue {residue!r} in {seq!r} has no alphabase "
                f"name among {candidates}; refusing to guess")
        mods.append(chosen)
        sites.append(str(site))
        pos = m.end()
    bare.append(seq[pos:])
    return "".join(bare), ";".join(mods), ";".join(sites)


def load_identifications(path, q_max):
    import pyarrow.parquet as pq
    index = build_unimod_index()
    per = {}
    pf = pq.ParquetFile(path)
    for b in pf.iter_batches(batch_size=200000,
                             columns=["Modified.Sequence", "Precursor.Charge",
                                      "RT", "Q.Value"]):
        d = b.to_pydict()
        for k in range(b.num_rows):
            if d["Q.Value"][k] > q_max:
                continue
            per.setdefault(d["Modified.Sequence"][k], []).append(float(d["RT"][k]))

    rows = []
    for seq, rts in per.items():
        rts.sort()
        # Median over charge states. They coelute to within 0.0018 min, so this
        # collapses a duplicate rather than averaging away a real difference.
        rt = rts[len(rts) // 2]
        sequence, mods, sites = parse_modified_sequence(seq, index)
        rows.append((sequence, mods, sites, rt))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("identifications", help="a DIA-NN report.parquet for THIS run")
    ap.add_argument("output", help="directory for rt.pth and the sidecar")
    ap.add_argument("--base-model", default="",
                    help="starting checkpoint; default is peptdeep's generic rt.pth")
    ap.add_argument("--max-peptides", type=int, default=2000,
                    help="training peptides. 500 buys 70%% of the gain in 31 s, "
                         "2000 buys 85%% in 78 s, and the full set costs 17 min "
                         "for the rest. Default 2000.")
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--q-value", type=float, default=0.01)
    args = ap.parse_args()

    import numpy as np
    import pandas as pd
    from peptdeep.pretrained_models import ModelManager

    t0 = time.time()
    rows = load_identifications(args.identifications, args.q_value)
    if len(rows) < 200:
        raise SystemExit(f"only {len(rows)} identifications at q <= {args.q_value}; "
                         f"too few to fine-tune on")
    print(f"{len(rows):,} distinct modified sequences at q <= {args.q_value}")

    rng = np.random.default_rng(20260803)
    if args.max_peptides and len(rows) > args.max_peptides:
        pick = rng.choice(len(rows), args.max_peptides, replace=False)
        rows = [rows[i] for i in pick]
    print(f"training on {len(rows):,}")

    df = pd.DataFrame(rows, columns=["sequence", "mods", "mod_sites", "rt"])
    df["nAA"] = df["sequence"].str.len()
    lo, hi = float(df["rt"].min()), float(df["rt"].max())
    if not hi > lo:
        raise SystemExit("every identification has the same retention time")
    # The model's target. Recorded in the sidecar because it is the scale the
    # exported model predicts on, and nothing downstream can recover it.
    df["rt_norm"] = (df["rt"] - lo) / (hi - lo)

    mgr = ModelManager(mask_modloss=False)
    mgr.load_installed_models()
    if args.base_model:
        mgr.rt_model.load(args.base_model)
    mgr.epoch_to_train_rt_ccs = args.epochs
    mgr.train_rt_model(df)

    os.makedirs(args.output, exist_ok=True)
    pth = os.path.join(args.output, "rt.pth")
    mgr.rt_model.save(pth)

    digest = hashlib.sha256(open(pth, "rb").read()).hexdigest()
    sidecar = {
        "model": os.path.abspath(pth),
        "sha256": digest,
        "tuned_from": os.path.abspath(args.identifications),
        "peptides": len(rows),
        "epochs": args.epochs,
        "q_value": args.q_value,
        "rt_norm_min_minutes": lo,
        "rt_norm_max_minutes": hi,
        "seconds": round(time.time() - t0, 1),
        # The caller must decide whether this closes a loop; see the module
        # docstring. Recorded, not judged.
        "warning": "This model is specific to the run it was tuned on and must "
                   "not be reused across runs or gradients.",
    }
    with open(os.path.join(args.output, "rt_provenance.json"), "w") as f:
        json.dump(sidecar, f, indent=2)
    print(json.dumps(sidecar, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
