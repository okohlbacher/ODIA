#!/usr/bin/env python3
"""Rungs 1-2 of doc/56's baseline ladder: per-fragment ORDER STATISTICS + GBT.

Both reviewers named this, not the 19 scalars, as the bar a transformer has to
clear. Order statistics of the per-fragment quantities are permutation-invariant
by construction -- fragment order in the library is arbitrary -- and doc/45's
finding that interference is DISTRIBUTED across fragments is exactly why they
should carry most of what attention would find.

PRIMARY CONTRAST: positives against RT-SHIFTED versions of THE SAME precursors.
Same sequence, same charge, same product m/z, same library intensities, same
everything -- only the window moves 300 s. So the metadata-only leakage control
is guaranteed to be AUC 0.5 BY CONSTRUCTION rather than by an argument, and none
of the organism confounds that make entrapment negatives suspect can operate.

Splitting is by PROTEIN, not precursor: paralogues and shared peptides leak
across a precursor-level split (doc/56 A8).
"""
import sys, numpy as np, pyarrow.parquet as pq, hashlib
from sklearn.metrics import roc_auc_score
import xgboost as xgb

S = sys.argv[1] if len(sys.argv) > 1 else '/scratch/kohlbach/odia2x2'
rng = np.random.default_rng(20260824)

lab = pq.read_table(f'{S}/tensor_ih1_labels.parquet')
ids = np.array(lab.column('Precursor.Id').to_pylist())
dec = lab.column('Decoy').to_numpy()
label = np.array(lab.column('Label').to_pylist())

m2 = pq.read_table(f'{S}/tensor_ih1shift_meta.parquet')
ids2 = np.array(m2.column('Precursor.Id').to_pylist())
dec2 = m2.column('Decoy').to_numpy()
assert len(ids) == len(ids2) and (ids == ids2).all() and (dec == dec2).all(), \
    'the two arms are not row-aligned; the paired design depends on it'
print(f'{len(ids):,} rows, both arms row-aligned')

X0 = np.load(f'{S}/tensor_ih1_traces.npy', mmap_mode='r')
X1 = np.load(f'{S}/tensor_ih1shift_traces.npy', mmap_mode='r')
M  = np.load(f'{S}/tensor_ih1_mask.npy')

def features(T, mask, chunk=8192):
    """Per-fragment quantities, then sorted across fragments."""
    P, F, C = T.shape
    out = np.zeros((P, 4 * F + 5), dtype=np.float32)
    for a in range(0, P, chunk):
        b = min(a + chunk, P)
        t = np.asarray(T[a:b], dtype=np.float32)
        mk = mask[a:b][:, :, None]
        t = t * mk
        area = t.sum(2)
        mx   = t.max(2)
        apex = t.argmax(2).astype(np.float32)
        med  = np.median(t, axis=2)
        snr  = mx / np.maximum(med, 1e-3)
        tot  = t.sum(1)                                  # group consensus
        gap  = tot.argmax(1).astype(np.float32)          # group apex cycle
        # Correlation of each fragment to the consensus of the OTHERS, so a
        # fragment is never correlated with itself through the sum.
        oth = tot[:, None, :] - t
        tc = t - t.mean(2, keepdims=True)
        oc = oth - oth.mean(2, keepdims=True)
        num = (tc * oc).sum(2)
        den = np.sqrt((tc ** 2).sum(2) * (oc ** 2).sum(2))
        corr = np.where(den > 0, num / np.maximum(den, 1e-12), 0.0)
        off = np.abs(apex - gap[:, None])
        good = mask[a:b]
        def srt(v, fill):
            w = np.where(good, v, fill)
            return np.sort(w, axis=1)[:, ::-1]
        f = [srt(corr, -1.0), srt(-off, -1e3), srt(np.log1p(area), -1.0),
             srt(np.log1p(snr), -1.0)]
        g = np.stack([good.sum(1), np.log1p(area.sum(1)), gap,
                      np.log1p(tot.max(1)), (tot > 0).sum(1)], axis=1)
        out[a:b] = np.concatenate(f + [g], axis=1).astype(np.float32)
    return out

print('computing rung-2 features (positives arm)...')
F0 = features(X0, M)
print('computing rung-2 features (shifted arm)...')
F1 = features(X1, M)

# Protein-level split. corpus_lib carries Protein.Group; hash it so the split is
# deterministic and independent of any ordering.
lib = pq.read_table(f'{S}/corpus_lib.parquet', columns=['Precursor.Id', 'Decoy',
                                                        'Protein.Group'])
pgm = {(p, int(d)): str(g) for p, d, g in zip(lib.column('Precursor.Id').to_pylist(),
                                              lib.column('Decoy').to_numpy(),
                                              lib.column('Protein.Group').to_pylist())}
def fold(pid, d):
    g = pgm.get((pid, int(d)), pid)
    return int(hashlib.md5(g.encode()).hexdigest(), 16) % 10
folds = np.array([fold(p, d) for p, d in zip(ids, dec)])
test = folds < 3
print(f'protein-level split: train {(~test).sum():,}  test {test.sum():,}')

pos = label == 'pos'
print(f'positives {pos.sum():,}; each contributes one positive (own window) and '
      f'one negative (shifted window)')

def build(sel):
    X = np.concatenate([F0[sel], F1[sel]])
    y = np.concatenate([np.ones(sel.sum()), np.zeros(sel.sum())])
    return X, y

Xtr, ytr = build(pos & ~test)
Xte, yte = build(pos & test)
print(f'train {len(ytr):,}  test {len(yte):,}')

clf = xgb.XGBClassifier(n_estimators=400, max_depth=6, learning_rate=0.05,
                        subsample=0.8, colsample_bytree=0.8, n_jobs=32,
                        eval_metric='logloss', tree_method='hist')
clf.fit(Xtr, ytr)
p = clf.predict_proba(Xte)[:, 1]
print(f'\nRUNG 2  order statistics + GBT   AUC {roc_auc_score(yte, p):.4f}')

# TWO different controls, which a first pass conflated -- and the conflation
# announced itself, because the "metadata-only" control returned 0.594 where
# the paired design guarantees 0.5.
#
#   TRUE metadata: identical between a positive and its own shifted negative by
#   construction (same precursor, same fragments, only the window moves). Any
#   value above 0.5 is a broken pairing, not leakage.
#
#   INTENSITY AND COVERAGE: total area and how many cycles carry signal. These
#   are TRACE-derived and legitimately differ between the arms, so they are not
#   a leakage control at all -- they are the crude floor that any shape-based
#   model has to beat to have shown anything about shape.
true_meta = [4 * 12 + 0]                      # fragment count: identical by pairing
crude     = [4 * 12 + 1, 4 * 12 + 3, 4 * 12 + 4]   # total area, peak height, coverage
for name, cols, expect in (('true-metadata control', true_meta, '0.5 by construction'),
                           ('intensity+coverage floor', crude, 'the bar for "shape"')):
    a = np.concatenate([F0[pos & ~test][:, cols], F1[pos & ~test][:, cols]])
    b = np.concatenate([F0[pos & test][:, cols], F1[pos & test][:, cols]])
    c = xgb.XGBClassifier(n_estimators=200, max_depth=4, n_jobs=32,
                          eval_metric='logloss', tree_method='hist').fit(a, ytr)
    print(f'  {name:<32} AUC {roc_auc_score(yte, c.predict_proba(b)[:, 1]):.4f}   ({expect})')

ysh = ytr.copy(); rng.shuffle(ysh)
cs = xgb.XGBClassifier(n_estimators=200, max_depth=4, n_jobs=32,
                       eval_metric='logloss', tree_method='hist').fit(Xtr, ysh)
print(f'  label-shuffle control            AUC '
      f'{roc_auc_score(yte, cs.predict_proba(Xte)[:, 1]):.4f}   (0.5 expected)')

# By abundance: the deficit is at Q1, and a model strong only at Q5 adds nothing.
dn = pq.read_table('/scratch/kohlbach/xic/dn_xic.parquet',
                   columns=['Precursor.Id', 'Q.Value', 'Precursor.Quantity'])
import re
A = {'UniMod:4': 'Carbamidomethyl'}
nrm = lambda s: re.sub(r'\((.*?)\)', lambda m: '(' + A.get(m.group(1), m.group(1)) + ')', str(s))
qq = dn.column('Q.Value').to_numpy()
qmap = {nrm(a): b for a, b, k in zip(dn.column('Precursor.Id').to_pylist(),
                                     dn.column('Precursor.Quantity').to_numpy(),
                                     qq <= 0.01) if k}
ab = np.array([qmap.get(p, np.nan) for p in ids[pos & test]])
ab2 = np.concatenate([ab, ab])
ed = np.nanpercentile(ab, [0, 20, 40, 60, 80, 100]); ed[-1] *= 1.001
print(f'\n{"abundance":>22} {"n":>8} {"AUC":>7}')
for i in range(5):
    s = (ab2 >= ed[i]) & (ab2 < ed[i + 1])
    if s.sum() > 50 and len(np.unique(yte[s])) == 2:
        print(f'  Q{i+1} {ed[i]:>9,.0f}-{ed[i+1]:>9,.0f} {s.sum():>8,} '
              f'{roc_auc_score(yte[s], p[s]):>7.4f}')
