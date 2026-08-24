#!/usr/bin/env python3
"""The leakage floor for the ENTRAPMENT contrast, which has to be measured.

On the shifted contrast the metadata floor is 0.5 by construction: a positive
and its own shifted negative are the same precursor. Entrapment negatives are
DIFFERENT precursors of a DIFFERENT organism, so descriptors that cannot leak
there can leak here -- and the library's own relative intensities were already
measured to differ between the classes (top-fragment share 0.2404 human against
0.2165 Arabidopsis).

So: train on the descriptors ALONE, with no trace information whatsoever. That
number is the floor any trace result on this contrast has to beat, and it is the
number doc/56 A2 said had to be measured rather than argued.
"""
import sys, hashlib, numpy as np, pyarrow.parquet as pq
from sklearn.metrics import roc_auc_score
import xgboost as xgb

D = sys.argv[1] if len(sys.argv) > 1 else '/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'
lab = pq.read_table(f'{D}/tensor_s08_labels.parquet')
ids = np.array(lab.column('Precursor.Id').to_pylist())
dec = lab.column('Decoy').to_numpy()
label = np.array(lab.column('Label').to_pylist())
A0 = np.load(f'{D}/apex_s08.npy')
RI = np.load(f'{D}/desc_relint.npy'); PM = np.load(f'{D}/desc_prodmz.npy')
FC = np.load(f'{D}/desc_frcharge.npy'); OR = np.load(f'{D}/desc_ordinal.npy')
SE = np.load(f'{D}/desc_series.npy');  M = np.load(f'{D}/tensor_s08_mask.npy')

lib = pq.read_table(f'{D}/corpus_lib.parquet', columns=['Precursor.Id','Decoy','Protein.Group'])
pgm = {(p, int(d)): str(g) for p, d, g in zip(lib.column('Precursor.Id').to_pylist(),
                                              lib.column('Decoy').to_numpy(),
                                              lib.column('Protein.Group').to_pylist())}
fold = np.array([int(hashlib.md5(pgm.get((p, int(d)), p).encode()).hexdigest(), 16) % 10
                 for p, d in zip(ids, dec)])
ok = A0 >= 0
pos = (label == 'pos') & ok
neg = (label == 'ent') & ok
sel = pos | neg
y = pos[sel].astype(int)
X = np.concatenate([RI[sel], PM[sel] / 1000.0, FC[sel], OR[sel] / 20.0, SE[sel],
                    M[sel].astype(np.float32)], axis=1)
f = fold[sel]
tr, te = f >= 3, f < 3
print(f'{pos.sum():,} positives against {neg.sum():,} entrapment; '
      f'train {tr.sum():,} test {te.sum():,}')
c = xgb.XGBClassifier(n_estimators=400, max_depth=6, learning_rate=0.05, n_jobs=32,
                      subsample=0.8, colsample_bytree=0.8,
                      eval_metric='logloss', tree_method='hist').fit(X[tr], y[tr])
a = roc_auc_score(y[te], c.predict_proba(X[te])[:, 1])
print(f'\nLEAKAGE FLOOR, descriptors only, NO trace information:  AUC {a:.4f}')
print('  any trace model on this contrast must beat this, not 0.5')
