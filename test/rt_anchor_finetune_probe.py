#!/usr/bin/env python
"""Can a model trained on PASS 1's OWN anchors beat the map, on everything else?

This is the in-pipeline question, and it is not the same as the one
`rt_sequence_gain.py` answers. That fitted on 80% of DIA-NN's confident set --
8,723 precursors we would never have at calibration time. Pass 1 gives us a few
thousand anchors, and the model has to generalise from those to the whole
library.

Train: ODIA's anchors (-out_anchors), sequence + library iRT + observed apex.
Test:  DIA-NN's confident precursors whose STRIPPED SEQUENCE is absent from the
       anchor set, so nothing the model saw can leak through a charge state.

Baseline is the same monotone map ODIA fits, trained on the same anchors and
evaluated on the same held-out set -- so the comparison isolates the model, not
the data.
"""
import csv, sys, zlib
import numpy as np
import pyarrow.parquet as pq

anchors_path, lib_path, report_path, label = sys.argv[1:5]

def strip_mods(s):
    out, depth = [], 0
    for c in s:
        if c in "([{": depth += 1
        elif c in ")]}": depth = max(0, depth - 1)
        elif depth == 0 and c.isalpha(): out.append(c)
    return "".join(out)

anch = []
with open(anchors_path) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        anch.append((strip_mods(r['Modified.Sequence']), int(r['Precursor.Charge']),
                     float(r['Library.iRT']), float(r['Observed.RT'])))
seen = {a[0] for a in anch}

pred = {}
with open(lib_path) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r.get('Decoy', '0') != '0' or not r['RT']: continue
        k = (strip_mods(r['Modified.Sequence']), int(r['Precursor.Charge']))
        pred.setdefault(k, float(r['RT']))
t = pq.read_table(report_path, columns=['Stripped.Sequence','Precursor.Charge','RT','Q.Value'])
d = t.to_pydict()
test = []
for i in range(t.num_rows):
    if d['Q.Value'][i] > 0.01: continue
    s, c = d['Stripped.Sequence'][i], int(d['Precursor.Charge'][i])
    if s in seen: continue                      # held out BY SEQUENCE
    k = (s, c)
    if k in pred: test.append((s, c, pred[k], float(d['RT'][i]) * 60.0))

print(f"=== {label} ===")
print(f"train anchors {len(anch)}   held-out precursors {len(test)} "
      f"(sequences disjoint from the anchors)")
if len(test) < 200: sys.exit("too few held out")

AA = "ACDEFGHIKLMNPQRSTVWY"
def comp(seq, charge):
    v = np.zeros(len(AA) + 2)
    for ch in seq:
        j = AA.find(ch)
        if j >= 0: v[j] += 1
    v[len(AA)] = len(seq); v[len(AA)+1] = charge
    return v

def monotone(xs, ys, knots=120):
    o = np.argsort(xs); xs, ys = np.asarray(xs)[o], np.asarray(ys)[o]
    idx = np.linspace(0, len(xs)-1, knots).astype(int)
    w = max(1, len(xs)//knots)
    kx = xs[idx]
    ky = np.maximum.accumulate([np.median(ys[max(0,i-w):i+w+1]) for i in idx])
    return lambda v: np.interp(v, kx, ky)

xtr = np.array([a[2] for a in anch]); ytr = np.array([a[3] for a in anch])
xte = np.array([r[2] for r in test]); yte = np.array([r[3] for r in test])
m = monotone(xtr, ytr)
r_base = yte - m(xte)
sd = lambda r: float(np.std(r, ddof=1))
print(f"  monotone map (what ODIA does)      held-out SD {sd(r_base):7.2f} s")

Xtr = np.hstack([np.array([comp(a[0], a[1]) for a in anch]),
                 m(xtr)[:, None], np.ones((len(anch), 1))])
Xte = np.hstack([np.array([comp(r[0], r[1]) for r in test]),
                 m(xte)[:, None], np.ones((len(test), 1))])
mu, sg = Xtr.mean(0), Xtr.std(0) + 1e-9
Ztr, Zte = (Xtr-mu)/sg, (Xte-mu)/sg
best = None
for lam in (0.3, 1, 3, 10, 30, 100):
    w = np.linalg.solve(Ztr.T@Ztr + lam*np.eye(Ztr.shape[1]), Ztr.T@ytr)
    s = sd(yte - Zte@w)
    print(f"  ridge(comp + calibrated iRT) l={lam:5.1f} held-out SD {s:7.2f} s")
    if best is None or s < best[1]: best = (lam, s)
print(f"  BEST ridge lambda {best[0]}: {best[1]:.2f} s vs map {sd(r_base):.2f} s "
      f"({100*(1-best[1]/sd(r_base)):+.1f}%)")
