#!/usr/bin/env python
"""Does the FINE-TUNED model, plus the same calibration, beat the map alone?

The fine-tuner's own sidecar compares tuned against STOCK peptdeep with no
calibration, which is not the decision in front of us: ODIA does not use raw
peptdeep output, it uses the library's iRT through a monotone map. So the
comparison that decides whether to wire this in is

    library iRT  -> monotone map   (what ODIA does today)
    tuned iRT    -> monotone map   (what it would do)

both fitted on pass 1's anchors and both evaluated on precursors whose STRIPPED
SEQUENCE is absent from those anchors.

Usage: rt_tuned_vs_map.py <anchors.tsv> <library.tsv> <report.parquet> <rt.onnx>
"""
import csv, subprocess, sys, zlib
import numpy as np, pyarrow.parquet as pq

anchors_path, lib_path, report_path, onnx_path, predictor = sys.argv[1:6]
# Optional 7th argument: the parquet the MODEL was trained on. Its sequences
# must leave the test set too.
#
# This is not optional once bootstrapping starts. Round 1 trains on the anchors,
# so excluding anchor sequences is enough. Round 2 trains on the SEARCH's own
# 6,382 identifications, which overlap DIA-NN's confident set almost entirely --
# so a test set filtered only against the anchors is one the model has largely
# seen, and the number it produces is leakage, not generalisation.
trained_on = sys.argv[6] if len(sys.argv) > 6 else ""

def strip_mods(s):
    out, depth = [], 0
    for c in s:
        if c in "([{": depth += 1
        elif c in ")]}": depth = max(0, depth - 1)
        elif depth == 0 and c.isalpha(): out.append(c)
    return "".join(out)

def predict(seqs, chunk=400):
    """Through ODIA's own odia_predict_rt.

    NOT a hand-rolled encoding. The first version of this script guessed at
    peptdeep's input layout -- 1-based amino-acid index, zero mod features --
    and produced a held-out SD of 440 s, which looked like a catastrophic model
    and was a catastrophic script. ODIA already has an encoder validated
    against test/peptdeep_reference.py; use it.
    """
    out = []
    for b in range(0, len(seqs), chunk):
        part = seqs[b:b+chunk]
        r = subprocess.run([predictor, onnx_path, *part],
                           capture_output=True, text=True, check=True)
        vals = {}
        for line in r.stdout.splitlines():
            f = line.split('\t')
            if len(f) >= 2:
                try: vals[f[0]] = float(f[1])
                except ValueError: pass
        out.extend(vals.get(p, float('nan')) for p in part)
    return np.array(out, dtype=float)

anch = []
with open(anchors_path) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        anch.append((strip_mods(r['Modified.Sequence']), float(r['Library.iRT']),
                     float(r['Observed.RT'])))
seen = {a[0] for a in anch}
if trained_on:
    tt = pq.read_table(trained_on, columns=["Modified.Sequence"]).to_pydict()
    seen |= {strip_mods(x) for x in tt["Modified.Sequence"]}
    print(f"excluding {len(seen)} training sequences (anchors + {trained_on})")

pred = {}
with open(lib_path) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r.get('Decoy','0') != '0' or not r['RT']: continue
        pred.setdefault((strip_mods(r['Modified.Sequence']), int(r['Precursor.Charge'])),
                        float(r['RT']))
t = pq.read_table(report_path, columns=['Stripped.Sequence','Precursor.Charge','RT','Q.Value'])
d = t.to_pydict()
test = []
for i in range(t.num_rows):
    if d['Q.Value'][i] > 0.01: continue
    s, c = d['Stripped.Sequence'][i], int(d['Precursor.Charge'][i])
    if s in seen: continue
    if (s, c) in pred: test.append((s, pred[(s, c)], float(d['RT'][i])*60.0))
# THREE-WAY, not two. The held-out set stops being held out the moment a
# hyper-parameter is chosen on it, and that is exactly what happened: 40 versus
# 100 epochs was selected by comparing 31.14 against 28.30 on THESE rows, then
# 28.30 was reported. That is validation-set model selection and it biases the
# number optimistically.
#
# So the evaluation set is split again by stripped sequence: VALIDATION is where
# recipes may be compared, TEST is reported and must be looked at once. Both are
# printed every run, because a validation number quoted as a test number is the
# same error in a new costume.
val  = [r for r in test if (zlib.crc32(r[0].encode()) & 0xFFFFFFFF) % 2 == 0]
tst  = [r for r in test if (zlib.crc32(r[0].encode()) & 0xFFFFFFFF) % 2 == 1]
print(f"train anchors {len(anch)}   evaluation {len(test)} "
      f"-> validation {len(val)} / TEST {len(tst)} (split by sequence)")

def monotone(xs, ys, knots=120):
    o = np.argsort(xs); xs, ys = np.asarray(xs)[o], np.asarray(ys)[o]
    idx = np.linspace(0, len(xs)-1, knots).astype(int); w = max(1, len(xs)//knots)
    return (lambda kx, ky: (lambda v: np.interp(v, kx, ky)))(
        xs[idx], np.maximum.accumulate([np.median(ys[max(0,i-w):i+w+1]) for i in idx]))

sd = lambda r: float(np.std(r, ddof=1))
ytr = np.array([a[2] for a in anch])
xtr = np.array([a[1] for a in anch])
m = monotone(xtr, ytr)
ptr = predict([a[0] for a in anch])
mt = monotone(ptr, ytr)

for name, rows in (("validation", val), ("TEST", tst)):
    if len(rows) < 100:
        print(f"  {name}: only {len(rows)} rows -- too few to mean anything, not reported")
        continue
    y = np.array([r[2] for r in rows])
    xb = np.array([r[1] for r in rows])
    xt = predict([r[0] for r in rows])
    print(f"  {name:10s} library iRT -> map  SD {sd(y - m(xb)):7.2f} s"
          f"   |  TUNED -> map  SD {sd(y - mt(xt)):7.2f} s   (n={len(rows)})")
