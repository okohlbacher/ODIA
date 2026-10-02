#!/usr/bin/env python3
"""Exact TreeSHAP over the three tree models, via xgboost's own pred_contribs.

Not the `shap` package: xgboost computes exact Shapley values for tree ensembles
internally, so this is the exact quantity rather than a sampling approximation,
and it needs no extra dependency.

Three models, and the third is the interesting one:

  A  the 19 shipped sub-scores on the SHIFTED contrast  -- which scalars carry
     the signal ODIA already has
  B  per-fragment ORDER STATISTICS on the same contrast -- whether the signal
     sits in the best fragments or the worst ones, which is doc/45's question
     about whether interference is concentrated or distributed
  C  descriptors-only on the ENTRAPMENT contrast -- the model that scores 0.7060
     with no trace information at all. Its SHAP says WHICH descriptor leaks
     taxonomy, which decides whether the leak can be removed or only measured.
"""
import sys, hashlib, re, numpy as np, pyarrow.parquet as pq
import xgboost as xgb

D = '/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'
S = '/scratch/kohlbach/odia2x2'

lab = pq.read_table(f'{D}/tensor_ih1_labels.parquet')
ids = np.array(lab.column('Precursor.Id').to_pylist())
dec = lab.column('Decoy').to_numpy()
label = np.array(lab.column('Label').to_pylist())
A0 = np.load(f'{D}/apex_ih1.npy'); A1 = np.load(f'{D}/apex_ih1shift.npy')
lib = pq.read_table(f'{D}/corpus_lib.parquet', columns=['Precursor.Id','Decoy','Protein.Group'])
pgm = {(p, int(d)): str(g) for p, d, g in zip(lib.column('Precursor.Id').to_pylist(),
                                              lib.column('Decoy').to_numpy(),
                                              lib.column('Protein.Group').to_pylist())}
fold = np.array([int(hashlib.md5(pgm.get((p, int(d)), p).encode()).hexdigest(), 16) % 10
                 for p, d in zip(ids, dec)])

def shap_report(X, y, tr, te, names, title, top=12):
    m = xgb.XGBClassifier(n_estimators=300, max_depth=5, learning_rate=0.06, n_jobs=32,
                          subsample=0.8, colsample_bytree=0.8,
                          eval_metric='logloss', tree_method='hist').fit(X[tr], y[tr])
    from sklearn.metrics import roc_auc_score
    auc = roc_auc_score(y[te], m.predict_proba(X[te])[:, 1])
    dm = xgb.DMatrix(X[te], feature_names=[f'f{i}' for i in range(X.shape[1])])
    contrib = m.get_booster().predict(dm, pred_contribs=True)[:, :-1]   # drop bias
    mean_abs = np.abs(contrib).mean(0)
    signed = contrib.mean(0)
    order = np.argsort(-mean_abs)
    print(f'\n=== {title}   (AUC {auc:.4f}, {X.shape[1]} features) ===')
    print(f'{"feature":<34}{"mean|SHAP|":>12}{"mean SHAP":>12}  direction')
    for i in order[:top]:
        d = 'towards POSITIVE' if signed[i] > 0 else 'towards NEGATIVE'
        print(f'{names[i]:<34}{mean_abs[i]:>12.4f}{signed[i]:>12.4f}  {d}')
    tot = mean_abs.sum()
    print(f'top-3 share of total |SHAP|: {100*mean_abs[order[:3]].sum()/tot:.1f}%')
    return mean_abs, names

# ---------- A: the shipped scalars, shifted contrast ----------
SC0 = np.load(f'{D}/scal_ih1.npy'); SC1 = np.load(f'{D}/scal_ih1shift.npy')
OK0 = np.load(f'{D}/scal_ih1_ok.npy'); OK1 = np.load(f'{D}/scal_ih1shift_ok.npy')
cols = [str(c) for c in np.load(f'{D}/scal_cols.npy')]
pos = (label == 'pos') & (A0 >= 0) & (A1 >= 0) & OK0 & OK1
idx = np.flatnonzero(pos)
X = np.concatenate([SC0[idx], SC1[idx]]); y = np.concatenate([np.ones(len(idx)), np.zeros(len(idx))])
f2 = np.concatenate([fold[idx], fold[idx]])
shap_report(X, y, f2 >= 3, f2 < 3, cols, 'A. shipped sub-scores, SHIFTED contrast')

# ---------- C: descriptors only, entrapment contrast (the 0.7060 floor) ----------
RI = np.load(f'{D}/desc_relint.npy'); PM = np.load(f'{D}/desc_prodmz.npy')
FC = np.load(f'{D}/desc_frcharge.npy'); OR = np.load(f'{D}/desc_ordinal.npy')
SE = np.load(f'{D}/desc_series.npy'); M = np.load(f'{D}/tensor_ih1_mask.npy')
ok = A0 >= 0
p2 = (label == 'pos') & ok; n2 = (label == 'ent') & ok
sel = p2 | n2
blocks = [('relint', RI), ('prodmz', PM / 1000.0), ('frcharge', FC),
          ('ordinal', OR / 20.0), ('series', SE), ('mask', M.astype(np.float32))]
Xd = np.concatenate([b[sel] for _, b in blocks], axis=1)
names = [f'{n}[{k}]' for n, b in blocks for k in range(b.shape[1])]
yd = p2[sel].astype(int); fd = fold[sel]
ma, _ = shap_report(Xd, yd, fd >= 3, fd < 3, names,
                    'C. descriptors ONLY, ENTRAPMENT contrast (the 0.7060 taxonomy floor)')
print('\n  by descriptor BLOCK (summed mean|SHAP|):')
o = 0
tot = ma.sum()
for n, b in blocks:
    w = ma[o:o + b.shape[1]].sum(); o += b.shape[1]
    print(f'    {n:<10} {w:>8.4f}   {100*w/tot:>5.1f}%')
