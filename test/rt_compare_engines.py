#!/usr/bin/env python
"""ODIA against DIA-NN's and OpenSWATH's own recalibrated retention times.

Same footing for all three: the fraction of DIA-NN's confident precursors whose
TRUE apex lies within +/-W of the engine's own predicted RT.

DIA-NN: `Predicted.RT` is its run-refined prediction; the residual is
RT - Predicted.RT directly.

OpenSWATH: it aligns iRT with TransformationModelLowess, and the .osw records
per feature EXP_RT (observed) and NORM_RT (that observation mapped into iRT
space). The pairs therefore TRACE OSW's own nonlinear calibration, and
interpolating them gives its predicted RT for any library iRT -- without
refitting anything of our own.
"""
import sqlite3, sys
import numpy as np, pyarrow.parquet as pq

report, osw = sys.argv[1], (sys.argv[2] if len(sys.argv) > 2 else "")
WIN = (15, 30, 60, 120)

t = pq.read_table(report, columns=['Precursor.Id','RT','Predicted.RT','Q.Value'])
d = t.to_pydict()
obs, err = {}, []
for i in range(t.num_rows):
    if d['Q.Value'][i] > 0.01: continue
    rt = float(d['RT'][i]) * 60.0
    obs[d['Precursor.Id'][i]] = rt
    if d['Predicted.RT'][i] is not None:
        err.append(abs(rt - float(d['Predicted.RT'][i]) * 60.0))

def row(label, e):
    e = np.sort(np.asarray(e))
    n = len(e)
    cov = "  ".join(f"{100*np.mean(e <= w):>7.2f}%" for w in WIN)
    print(f"{label:<22} {n:>6}  {cov}  {e[n//2]:>6.1f}  {e[int(.95*(n-1))]:>6.1f}")

print(f"{'engine':<22} {'n':>6}  " + "  ".join(f"+/-{w:>4}s" for w in WIN) + "     p50     p95")
row("DIA-NN (own refit)", err)

if osw:
    c = sqlite3.connect(f"file:{osw}?mode=ro", uri=True)
    # The calibration itself, traced from the features: (NORM_RT -> EXP_RT).
    pairs = np.array(c.execute(
        "SELECT NORM_RT, EXP_RT FROM FEATURE WHERE NORM_RT IS NOT NULL "
        "AND EXP_RT IS NOT NULL").fetchall(), dtype=float)
    if len(pairs) < 1000:
        print("  OpenSWATH: too few features to trace its calibration"); sys.exit(0)
    o = np.argsort(pairs[:, 0]); xs, ys = pairs[o, 0], pairs[o, 1]
    # Bin to a monotone curve: many features share one iRT and only their median
    # traces the map rather than the interference around it.
    idx = np.linspace(0, len(xs) - 1, 400).astype(int)
    w = max(1, len(xs) // 400)
    kx = xs[idx]
    ky = np.maximum.accumulate([np.median(ys[max(0, i - w):i + w + 1]) for i in idx])
    lib = {r[0]: r[1] for r in c.execute(
        "SELECT p.ID, p.LIBRARY_RT FROM PRECURSOR p WHERE p.DECOY = 0 "
        "AND p.LIBRARY_RT IS NOT NULL")}
    pid = {r[0]: r[1] for r in c.execute(
        "SELECT PRECURSOR_ID, PEPTIDE_ID FROM PRECURSOR_PEPTIDE_MAPPING")}
    seq = {r[0]: r[1] for r in c.execute("SELECT ID, MODIFIED_SEQUENCE FROM PEPTIDE")}
    chg = {r[0]: r[1] for r in c.execute("SELECT ID, CHARGE FROM PRECURSOR")}
    oerr = []
    for prec, lrt in lib.items():
        key = seq.get(pid.get(prec, -1))
        if key is None: continue
        k = f"{key}{chg.get(prec,'')}"
        if k not in obs: continue
        oerr.append(abs(float(np.interp(lrt, kx, ky)) - obs[k]))
    if oerr: row("OpenSWATH (lowess)", oerr)
    else: print("  OpenSWATH: no precursor ids matched DIA-NN's")
