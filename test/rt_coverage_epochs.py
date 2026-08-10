#!/usr/bin/env python
"""Does +/-30 s coverage budge with more fine-tuning epochs?

Needs no extraction. For each ONNX model: predict iRT for the anchors and for
the whole library, fit the same monotone map ODIA fits (binned medians + PAVA)
on (predicted anchor iRT, observed anchor RT), apply it to every precursor, and
compare against DIA-NN's observed RT.

Coverage, not residual SD. A model can tighten the bulk without moving the
fraction of peptides whose peak is actually inside the window, and the window is
what the search opens.

Usage: rt_coverage_epochs.py <anchors.tsv> <library.tsv> <report.parquet>
                            <predictor> <label>=<onnx> [...]
"""
import csv, subprocess, sys
import numpy as np, pyarrow.parquet as pq

anchors_path, lib_path, report_path, predictor = sys.argv[1:5]
models = [a.split('=', 1) for a in sys.argv[5:]]

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
        anch.append((strip_mods(r['Modified.Sequence']), float(r['Observed.RT'])))

lib = {}
csv.field_size_limit(10**9)
with open(lib_path) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r.get('Decoy', '0') != '0' or not r.get('RT'): continue
        lib.setdefault((strip_mods(r['Modified.Sequence']), int(r['Precursor.Charge'])),
                       float(r['RT']))

t = pq.read_table(report_path, columns=['Stripped.Sequence','Precursor.Charge','RT','Q.Value'])
d = t.to_pydict()
obs = {}
for i in range(t.num_rows):
    if d['Q.Value'][i] > 0.01: continue
    obs[(d['Stripped.Sequence'][i], int(d['Precursor.Charge'][i]))] = float(d['RT'][i]) * 60.0
keys = sorted(set(lib) & set(obs))
print(f"anchors {len(anch)}   evaluable precursors {len(keys)}")

def predict(onnx, seqs, chunk=400):
    out, uniq = {}, sorted(set(seqs))
    for b in range(0, len(uniq), chunk):
        part = uniq[b:b+chunk]
        r = subprocess.run([predictor, onnx, *part], capture_output=True, text=True, check=True)
        for line in r.stdout.splitlines():
            f = line.split('\t')
            if len(f) >= 2:
                try: out[f[0]] = float(f[1])
                except ValueError: pass
    return [out.get(s, float('nan')) for s in seqs]

def monotone(xs, ys, knots=120):
    o = np.argsort(xs); xs, ys = np.asarray(xs)[o], np.asarray(ys)[o]
    idx = np.linspace(0, len(xs)-1, knots).astype(int); w = max(1, len(xs)//knots)
    kx = xs[idx]
    ky = np.maximum.accumulate([np.median(ys[max(0,i-w):i+w+1]) for i in idx])
    return lambda v: np.interp(v, kx, ky)

WIN = (15, 30, 60, 120)
print(f"{'model':<18} " + "  ".join(f"+/-{w:>4}s" for w in WIN) + "     p50     p95")
# Baseline: the library's own iRT through the same map.
a_seq = [a[0] for a in anch]
a_rt = np.array([a[1] for a in anch])
base_m = monotone(np.array([lib.get((s, 2), np.nan) for s in a_seq]), a_rt)
for label, onnx in [("library iRT", None)] + models:
    if onnx is None:
        xs = np.array([lib[k] for k in keys])
        m = monotone(np.array([lib.get((s, 2), lib.get((s, 3), np.nan)) for s in a_seq]), a_rt)
    else:
        pa = np.array(predict(onnx, a_seq))
        ok = np.isfinite(pa)
        m = monotone(pa[ok], a_rt[ok])
        pk = np.array(predict(onnx, [k[0] for k in keys]))
        xs = pk
    pred = m(xs)
    err = np.abs(pred - np.array([obs[k] for k in keys]))
    err = err[np.isfinite(err)]
    e = np.sort(err)
    cov = [100.0 * np.mean(err <= w) for w in WIN]
    print(f"{label:<18} " + "  ".join(f"{c:>7.2f}%" for c in cov) +
          f"  {e[len(e)//2]:>6.1f}  {e[int(.95*(len(e)-1))]:>6.1f}")
