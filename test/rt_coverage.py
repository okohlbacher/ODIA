#!/usr/bin/env python
"""Does the extraction window CONTAIN the peak? For every precursor, not just ours.

This is the question the calibration exists to answer, and the anchor-based
residual cannot answer it. Anchors are precursors we already identified, so a
p99 over them describes the peptides we find and is silent about the ones we
miss -- which is the entire population at risk.

So: take the library ODIA writes AFTER calibration (-out_lib, whose RT column is
predicted RUN SECONDS), and compare it against DIA-NN's observed RT for the same
precursors. Coverage at half-width W is the fraction whose true apex lies within
+/- W of our prediction. That is a real miss rate, over 11,112 precursors of
which we anchor ~2,500.

Usage: rt_coverage.py <label> <predicted_lib.tsv> <report.parquet> [more triples]
"""
import csv, sys
import pyarrow.parquet as pq

def strip_mods(s):
    out, depth = [], 0
    for c in s:
        if c in "([{": depth += 1
        elif c in ")]}": depth = max(0, depth - 1)
        elif depth == 0 and c.isalpha(): out.append(c)
    return "".join(out)

def load_pred(path):
    pred = {}
    csv.field_size_limit(10**9)
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r.get('Decoy', '0') != '0' or not r.get('RT'): continue
            k = (strip_mods(r['Modified.Sequence']), int(r['Precursor.Charge']))
            pred.setdefault(k, float(r['RT']))
    return pred

def load_obs(path):
    t = pq.read_table(path, columns=['Stripped.Sequence','Precursor.Charge','RT','Q.Value'])
    d = t.to_pydict()
    obs = {}
    for i in range(t.num_rows):
        if d['Q.Value'][i] > 0.01: continue
        obs[(d['Stripped.Sequence'][i], int(d['Precursor.Charge'][i]))] = float(d['RT'][i]) * 60.0
    return obs

args = sys.argv[1:]
print(f"{'config':<16} {'n':>6}  " + "  ".join(f"+/-{w:>4}s" for w in (15,30,60,120,240,480)))
for label, lib, rep in zip(args[0::3], args[1::3], args[2::3]):
    pred, obs = load_pred(lib), load_obs(rep)
    keys = sorted(set(pred) & set(obs))
    if not keys:
        print(f"{label:<16} no overlap"); continue
    err = sorted(abs(pred[k] - obs[k]) for k in keys)
    n = len(err)
    cov = [sum(1 for e in err if e <= w) / n for w in (15,30,60,120,240,480)]
    print(f"{label:<16} {n:>6}  " + "  ".join(f"{100*c:>7.2f}%" for c in cov))
    q = lambda p: err[min(n-1, int(p*(n-1)))]
    print(f"{'':<16} {'':>6}  p50 {q(.50):.1f}s  p90 {q(.90):.1f}s  p95 {q(.95):.1f}s  "
          f"p99 {q(.99):.1f}s  max {err[-1]:.1f}s")
