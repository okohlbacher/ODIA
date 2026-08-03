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

**Two mechanisms, both kept.** `--method direct` retrains the model to predict
normalised retention time. `--method residual` fits a calibration first --
linear, then isotonic on top of it -- and retrains only on what the calibration
could not explain, adding the two back at prediction time.

On S08 direct won by a wide margin: held-out sd 0.429 against 0.635 at 500
peptides and 500 epochs, with residual barely beating the calibration it sits
on. **That is one file, and one file is not a result.** The likely reason is
that retargeting pretrained weights at a small centred residual fights the
initialisation rather than exploiting it -- but a different gradient, instrument
or organism could reverse it, and freezing the trunk (untested) could too. Both
mechanisms are therefore selectable, and `--evaluate` makes every run report
which one won on ITS data rather than inheriting a verdict from ours.

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
    ap.add_argument("--method", choices=("direct", "residual"), default="direct",
                    help="direct retrains on normalised RT; residual calibrates "
                         "first and retrains on the remainder. direct won on S08 "
                         "by 0.429 against 0.635, on one file. Use --evaluate.")
    ap.add_argument("--evaluate", action="store_true",
                    help="hold out a fifth of the peptides, and report held-out "
                         "residual for the stock model and for this method. Costs "
                         "a fifth of the training data and answers the only "
                         "question that matters on a new file: did this help?")
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

    # Held out BY SEQUENCE, so a peptide's charge states cannot straddle the
    # split. Deterministic, so re-running reports the same number.
    import zlib
    held = np.array([(zlib.crc32(s.encode()) & 0xFFFFFFFF) % 5 == 0
                     for s in df["sequence"]]) if args.evaluate else np.zeros(len(df), bool)
    train_df = df[~held].reset_index(drop=True)
    test_df = df[held].reset_index(drop=True)
    if args.evaluate and len(test_df) < 50:
        raise SystemExit("too few peptides to hold any out; drop --evaluate")

    lo, hi = float(train_df["rt"].min()), float(train_df["rt"].max())
    if not hi > lo:
        raise SystemExit("every identification has the same retention time")

    mgr = ModelManager(mask_modloss=False)
    mgr.load_installed_models()
    if args.base_model:
        mgr.rt_model.load(args.base_model)
    mgr.epoch_to_train_rt_ccs = args.epochs

    stock_pred = mgr.rt_model.predict(df.copy())["rt_pred"].values
    calibration = None

    if args.method == "direct":
        # Absolute normalised retention time; the exported model predicts it
        # directly and ODIA consumes it with -rt_model.
        train_df = train_df.assign(rt_norm=(train_df["rt"] - lo) / (hi - lo))
    else:
        # Calibrate first -- linear, then isotonic -- and learn the remainder.
        from sklearn.isotonic import IsotonicRegression
        sp = stock_pred[~held]
        A = np.vstack([sp, np.ones_like(sp)]).T
        m0, c0 = np.linalg.lstsq(A, train_df["rt"].values, rcond=None)[0]
        iso = IsotonicRegression(increasing=True, out_of_bounds="clip").fit(
            m0 * sp + c0, train_df["rt"].values)
        cal_tr = iso.predict(m0 * sp + c0)
        resid = train_df["rt"].values - cal_tr
        rlo, rhi = float(resid.min()), float(resid.max())
        train_df = train_df.assign(rt_norm=(resid - rlo) / (rhi - rlo))
        calibration = dict(slope=float(m0), intercept=float(c0),
                           knots_x=[float(x) for x in iso.X_thresholds_],
                           knots_y=[float(y) for y in iso.y_thresholds_],
                           residual_min=rlo, residual_max=rhi)

    mgr.train_rt_model(train_df)

    os.makedirs(args.output, exist_ok=True)
    pth = os.path.join(args.output, "rt.pth")
    mgr.rt_model.save(pth)

    evaluation = None
    if args.evaluate:
        tuned = mgr.rt_model.predict(df.copy())["rt_pred"].values

        def sd(pred, truth, tr_mask, te_mask):
            p, t = pred[tr_mask], truth[tr_mask]
            A = np.vstack([p, np.ones_like(p)]).T
            m, c = np.linalg.lstsq(A, t, rcond=None)[0]
            return float(np.std(truth[te_mask] - (m * pred[te_mask] + c)))

        truth = df["rt"].values
        if args.method == "residual":
            # A residual model predicts a CORRECTION. Scoring its raw output as
            # though it were a retention time reports a catastrophe that is an
            # artefact of the scoring, and would condemn the method on every
            # dataset. The prediction is calibration(stock) + correction.
            cal_all = iso.predict(m0 * stock_pred + c0)
            combined = cal_all + (tuned * (rhi - rlo) + rlo)
            tuned_sd = float(np.std(truth[held] - combined[held]))
        else:
            tuned_sd = sd(tuned, truth, ~held, held)

        evaluation = {
            "held_out_peptides": int(held.sum()),
            "stock_sd_minutes": round(sd(stock_pred, truth, ~held, held), 4),
            "tuned_sd_minutes": round(tuned_sd, 4),
            "method": args.method,
        }
        evaluation["improvement"] = round(
            1.0 - evaluation["tuned_sd_minutes"] / evaluation["stock_sd_minutes"], 4)
        # Said plainly, because the alternative is a user inheriting our verdict
        # from a different instrument.
        if evaluation["tuned_sd_minutes"] >= evaluation["stock_sd_minutes"]:
            print(f"WARNING: --method {args.method} did NOT improve on this data "
                  f"({evaluation['tuned_sd_minutes']:.3f} against "
                  f"{evaluation['stock_sd_minutes']:.3f} min held out). Try the "
                  f"other method before using this model.", file=sys.stderr)

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
        "method": args.method,
        "seconds": round(time.time() - t0, 1),
        "evaluation": evaluation,
        # Residual mode's model predicts a CORRECTION, not a retention time, so
        # it is not consumable with -rt_model alone; the calibration below has
        # to be applied and added back. ODIA cannot do that yet.
        "calibration": calibration,
        # The caller must decide whether this closes a loop; see the module
        # docstring. Recorded, not judged.
        "warning": "This model is specific to the run it was tuned on and must "
                   "not be reused across runs or gradients.",
    }
    with open(os.path.join(args.output, "rt_provenance.json"), "w") as f:
        json.dump(sidecar, f, indent=2)
    if args.method == "residual":
        print("NOTE: a residual model predicts a correction, not a retention "
              "time. ODIA cannot consume it with -rt_model -- the calibration in "
              "the sidecar must be applied and added back first. Use it to "
              "compare methods, not yet to build a library.", file=sys.stderr)
    print(json.dumps(sidecar, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
