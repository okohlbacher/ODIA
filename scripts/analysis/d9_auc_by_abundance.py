#!/usr/bin/env python3
"""D9: which sub-scores stop working when the precursor gets faint?

The funnel says the deficit is discrimination at low abundance -- acceptance
runs 7.0% -> 81.1% with abundance among precursors we already find, and 81% of
oracle-admitted precursors still fail. If that is right, the sub-scores should
lose their separation as abundance falls, and WHICH ones lose it says what to
fix.

Positives: DIA-NN-confident precursors whose top-DScore candidate is within
tol of DIA-NN's retention time -- stage D of doc/51, so the peak really is the
right one and a low AUC means the FEATURE failed, not the peak-picking.
Negatives: every decoy precursor's best-scoring row, which is the population the
q-value threshold actually competes against.

Decoys have no abundance, so each bin is "targets this faint against the WHOLE
decoy null". That is the operational question, not a matched-abundance one.
"""
import sys, re, numpy as np, pyarrow.parquet as pq

odia_tsv, dn_pq, lib_pq = sys.argv[1], sys.argv[2], sys.argv[3]
# 'best' (default) or 'random'. Selecting the decoy row by DScore selects decoys
# to look good on every feature the model weighs, which can manufacture an AUC
# below 0.5 for a sound feature. A feature that straightens out under 'random'
# was an artefact of the argmax, not a defect.
DECOY_PICK = sys.argv[4] if len(sys.argv) > 4 else 'best'
TOL = 20.0
ALIAS = {'UniMod:4': 'Carbamidomethyl'}
def norm(p):
    return re.sub(r'\((.*?)\)', lambda m: '(' + ALIAS.get(m.group(1), m.group(1)) + ')', str(p))

dn = pq.read_table(dn_pq, columns=['Precursor.Id', 'RT', 'Q.Value', 'Precursor.Quantity'])
q = dn.column('Q.Value').to_numpy()
ids = [norm(p) for p in dn.column('Precursor.Id').to_pylist()]
rt = dn.column('RT').to_numpy() * 60.0
qt = dn.column('Precursor.Quantity').to_numpy()
conf = {ids[i]: (rt[i], qt[i]) for i in np.where(q <= 0.01)[0]}
lib = pq.read_table(lib_pq, columns=['Precursor.Id', 'Decoy'])
L = {p for p, d in zip(lib.column('Precursor.Id').to_pylist(),
                       lib.column('Decoy').to_numpy().astype(bool)) if not d and p in conf}

rng = np.random.default_rng(0)
nseen = {}
h = None
pos, neg = {}, {}          # id -> (dscore, rt, feature vector)
with open(odia_tsv) as fh:
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if h is None:
            h = {n: i for i, n in enumerate(f)}
            cols = [n for n in f if n.startswith('var_')] + ['DScore']
            idx = [h[n] for n in cols]
            continue
        try:
            ds = float(f[h['DScore']])
        except ValueError:
            continue
        pid = f[h['Precursor.Id']]
        tgt = f[h['Decoy']] in ('0', 'false', 'False')
        if tgt and pid in L:
            if pid not in pos or ds > pos[pid][0]:
                pos[pid] = (ds, float(f[h['RT']]),
                            [float(f[k]) if f[k] not in ('', 'nan') else np.nan for k in idx])
        elif not tgt:
            row = (ds, 0.0,
                   [float(f[k]) if f[k] not in ('', 'nan') else np.nan for k in idx])
            if DECOY_PICK == 'random':
                nseen[pid] = nseen.get(pid, 0) + 1
                if rng.random() < 1.0 / nseen[pid]:   # reservoir of one, single pass
                    neg[pid] = row
            elif pid not in neg or ds > neg[pid][0]:
                neg[pid] = row

keep = [p for p in pos if abs(pos[p][1] - conf[p][0]) <= TOL]
print(f'positives (stage D): {len(keep):,}   decoy negatives: {len(neg):,}   [decoy row = {DECOY_PICK}]')
P = np.array([pos[p][2] for p in keep])
N = np.array([v[2] for v in neg.values()])
A = np.array([conf[p][1] for p in keep])

def auc(a, b):
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if len(a) < 20 or len(b) < 20: return np.nan
    s = np.sort(b)
    r = (np.searchsorted(s, a, 'left') + np.searchsorted(s, a, 'right')) / 2.0
    return float(np.mean(r) / len(s))

edges = np.percentile(A, [0, 20, 40, 60, 80, 100]); edges[-1] *= 1.0001
hdr = ''.join(f'{f"Q{i+1}":>8}' for i in range(len(edges) - 1))
print(f'\n{"sub-score":>22}{hdr}{"all":>8}{"drop":>8}')
rows = []
for c, name in enumerate(cols):
    vals = [auc(P[(A >= edges[i]) & (A < edges[i + 1]), c], N[:, c])
            for i in range(len(edges) - 1)]
    allv = auc(P[:, c], N[:, c])
    if not np.isfinite(allv): continue
    drop = (vals[-1] - vals[0]) if np.isfinite(vals[0]) and np.isfinite(vals[-1]) else np.nan
    rows.append((name, vals, allv, drop))
for name, vals, allv, drop in sorted(rows, key=lambda r: -(r[3] if np.isfinite(r[3]) else -9)):
    v = ''.join(f'{x:8.3f}' if np.isfinite(x) else f'{"-":>8}' for x in vals)
    print(f'{name:>22}{v}{allv:8.3f}{drop:8.3f}')
print(f'\nQ1 is the faintest fifth (DIA-NN quantity < {edges[1]:,.0f}), '
      f'Q5 the brightest (> {edges[-2]:,.0f}).')
print('"drop" is Q5 - Q1: how much of the feature\'s power is abundance.')
