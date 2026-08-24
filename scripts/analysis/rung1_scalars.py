#!/usr/bin/env python3
"""Rung 1: ODIA's shipped 19 sub-scores, on the SAME contrast and split as rung 2.

The comparison only means something if the split, the labels and the contrast
are identical, so this reuses all three. Each precursor contributes its own
window (positive) and its +300 s window (negative), scored by the production
scorer in both arms.
"""
import sys, hashlib, numpy as np, pyarrow.parquet as pq
from sklearn.metrics import roc_auc_score
import xgboost as xgb

S = '/scratch/kohlbach/odia2x2'
def best_rows(path):
    """Best-DScore candidate per (Precursor.Id, Decoy), with its sub-scores."""
    out, h = {}, None
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if h is None:
            h = {n: i for i, n in enumerate(f)}
            cols = [n for n in f if n.startswith('var_')]
            idx = [h[c] for c in cols]
            continue
        try: ds = float(f[h['DScore']])
        except ValueError: continue
        k = (f[h['Precursor.Id']], f[h['Decoy']])
        if k not in out or ds > out[k][0]:
            v = []
            for j in idx:
                try: v.append(float(f[j]))
                except ValueError: v.append(np.nan)
            out[k] = (ds, v)
    return out, cols

own,  cols = best_rows(f'{S}/corpus_scores_s08.tsv')
shift, _   = best_rows(f'{S}/corpus_scores_s08_shift.tsv')
print(f'{len(cols)} sub-scores; scored precursors own {len(own):,} shifted {len(shift):,}')

lab = pq.read_table(f'{S}/tensor_s08_labels.parquet')
ids = np.array(lab.column('Precursor.Id').to_pylist())
dec = lab.column('Decoy').to_numpy()
label = np.array(lab.column('Label').to_pylist())

lib = pq.read_table(f'{S}/corpus_lib.parquet', columns=['Precursor.Id','Decoy','Protein.Group'])
pgm = {(p, int(d)): str(g) for p, d, g in zip(lib.column('Precursor.Id').to_pylist(),
                                              lib.column('Decoy').to_numpy(),
                                              lib.column('Protein.Group').to_pylist())}
fold = lambda p, d: int(hashlib.md5(pgm.get((p, int(d)), p).encode()).hexdigest(), 16) % 10

X, y, te = [], [], []
kept = miss = 0
for p, d, l in zip(ids, dec, label):
    if l != 'pos': continue
    k = (p, str(int(d)))
    if k not in own or k not in shift: miss += 1; continue
    kept += 1
    t = fold(p, d) < 3
    X.append(own[k][1]);   y.append(1); te.append(t)
    X.append(shift[k][1]); y.append(0); te.append(t)
X = np.array(X, dtype=np.float32); y = np.array(y); te = np.array(te)
print(f'positives usable {kept:,}  (dropped {miss:,} lacking a candidate in one arm)')
print(f'train {(~te).sum():,}  test {te.sum():,}')

clf = xgb.XGBClassifier(n_estimators=400, max_depth=6, learning_rate=0.05,
                        subsample=0.8, colsample_bytree=0.8, n_jobs=32,
                        eval_metric='logloss', tree_method='hist')
clf.fit(X[~te], y[~te])
pr = clf.predict_proba(X[te])[:, 1]
print(f'\nRUNG 1  shipped 19 sub-scores + GBT   AUC {roc_auc_score(y[te], pr):.4f}')
imp = sorted(zip(cols, clf.feature_importances_), key=lambda t: -t[1])[:6]
print('  top features: ' + ', '.join(f'{n}={v:.3f}' for n, v in imp))
