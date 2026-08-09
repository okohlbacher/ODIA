"""Does SEQUENCE information beyond the library iRT beat the calibration ceiling?

The whole premise of fine-tuning is that the residual left after the best
monotone map of library iRT is elution-ORDER error that only a sequence-aware
model can touch. That is testable without peptdeep or torch: fit ridge on
amino-acid composition PLUS the calibrated iRT, held out by stripped sequence,
and see how far below the monotone ceiling it goes.

If it does not move, fine-tuning will not either and the cluster time is saved.
"""
import sys, math, zlib, csv
import numpy as np
import pyarrow.parquet as pq

sys.path.insert(0, 'test')
LIB = '/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/lib/astral_lib_own.tsv'
REP = '/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/lib/astral_truth.parquet'

def strip_mods(s):
    out, depth = [], 0
    for c in s:
        if c in "([{": depth += 1
        elif c in ")]}": depth = max(0, depth-1)
        elif depth == 0 and c.isalpha(): out.append(c)
    return "".join(out)

# library iRT per (stripped seq, charge)
pred = {}
with open(LIB) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        k = (strip_mods(r['Modified.Sequence']), int(r['Precursor.Charge']))
        if r.get('Decoy','0') != '0' or not r['RT']: continue
        if k not in pred: pred[k] = float(r['RT'])
t = pq.read_table(REP, columns=['Stripped.Sequence','Precursor.Charge','RT','Q.Value'])
d = t.to_pydict()
obs = {}
for i in range(t.num_rows):
    if d['Q.Value'][i] > 0.01: continue
    obs[(d['Stripped.Sequence'][i], int(d['Precursor.Charge'][i]))] = float(d['RT'][i])

keys = sorted(set(pred) & set(obs))
print(f"paired precursors: {len(keys)}")
if len(keys) < 1000: sys.exit("too few")

AA = "ACDEFGHIKLMNPQRSTVWY"
def feats(seq, charge, irt):
    v = np.zeros(len(AA)+4)
    for c in seq:
        j = AA.find(c)
        if j >= 0: v[j] += 1
    v[len(AA)]   = len(seq)
    v[len(AA)+1] = charge
    v[len(AA)+2] = irt
    v[len(AA)+3] = 1.0
    return v

X = np.array([feats(k[0], k[1], pred[k]) for k in keys])
y = np.array([obs[k] for k in keys])
h = np.array([(zlib.crc32(k[0].encode()) & 0xFFFFFFFF) % 5 for k in keys])
tr, te = h != 0, h == 0
print(f"train {tr.sum()}  held out {te.sum()}  (split by stripped sequence)")

def sd(r): return float(np.std(r, ddof=1))

# 1. monotone baseline: piecewise-linear through quantile knots on iRT only
def monotone(xs, ys, knots=200):
    o = np.argsort(xs); xs, ys = xs[o], ys[o]
    qs = np.linspace(0, len(xs)-1, knots).astype(int)
    kx = xs[qs]
    ky = np.array([np.median(ys[max(0,i-len(xs)//knots):i+len(xs)//knots+1]) for i in qs])
    ky = np.maximum.accumulate(ky)
    return lambda v: np.interp(v, kx, ky)
m = monotone(X[tr, len(AA)+2], y[tr])
r_mono = y[te] - m(X[te, len(AA)+2])
print(f"monotone on iRT only   held-out SD {sd(r_mono)*60:7.2f} s   ({sd(r_mono):.4f} min)")

# 2. ridge on composition + calibrated iRT
Xc = X.copy()
Xc[:, len(AA)+2] = m(X[:, len(AA)+2])          # feed the CALIBRATED iRT
mu, sg = Xc[tr].mean(0), Xc[tr].std(0) + 1e-9
Z = (Xc - mu) / sg
for lam in (1.0, 10.0, 100.0):
    A = Z[tr].T @ Z[tr] + lam*np.eye(Z.shape[1])
    w = np.linalg.solve(A, Z[tr].T @ y[tr])
    r = y[te] - Z[te] @ w
    print(f"ridge(comp+calib iRT) lam={lam:5.0f} held-out SD {sd(r)*60:7.2f} s   ({sd(r):.4f} min)")
