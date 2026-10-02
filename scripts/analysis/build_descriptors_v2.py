#!/usr/bin/env python3
"""Per-fragment descriptors for BOTH classes, from the library ODIA actually used.

The first version built descriptors from `corpus_lib.parquet` and matched only
targets: every decoy row came out with all-zero library intensities and no
precursor charge. Training any target-vs-decoy contrast on that would have
learned one fact -- descriptors present means target -- scored near-perfect, and
passed every leakage control, because those controls were built for the shifted
contrast where both arms are the SAME precursor.

Two mismatches caused it and both are fixed here:

  * ODIA REGENERATES decoys internally, so they are absent from the input
    library. `-stop_after library -out_lib` writes the library it actually used.
  * That library names a decoy `<target_id>_decoy`, while the trace dump
    reconstructs the id from modified_sequence + charge and so writes the BARE
    target id -- exactly the collision ChromatogramTsv.cpp's own comment warns
    about. Strip the suffix and key on (id, decoy_flag).

The written library is FLAT per-transition, unlike the per-precursor list format
of the input, so transitions are grouped by precursor in the order they appear
(which is the library's own transition order, matching the tensor's slot).
"""
import sys, collections, numpy as np, pyarrow.parquet as pq

SRC = sys.argv[1] if len(sys.argv) > 1 else '/scratch/kohlbach/odia2x2/corpus_lib_withdecoys.parquet'
D   = sys.argv[2] if len(sys.argv) > 2 else '/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'
F = 12

meta = pq.read_table(f'{D}/tensor_ih1_meta.parquet')
ids = meta.column('Precursor.Id').to_pylist()
dec = meta.column('Decoy').to_numpy()
pos_of = {(p, int(d)): i for i, (p, d) in enumerate(zip(ids, dec))}
P = len(ids)

lib = pq.read_table(SRC, columns=['Precursor.Id', 'Decoy', 'Precursor.Charge',
                                  'Precursor.Mz', 'Modified.Sequence', 'Product.Mz',
                                  'Relative.Intensity', 'Fragment.Charge',
                                  'Fragment.Series.Number', 'Fragment.Type'])
lid = lib.column('Precursor.Id').to_pylist()
ldec = lib.column('Decoy').to_numpy().astype(int)
pm = lib.column('Product.Mz').to_numpy(); ri = lib.column('Relative.Intensity').to_numpy()
fc = lib.column('Fragment.Charge').to_numpy(); sn = lib.column('Fragment.Series.Number').to_numpy()
ft = lib.column('Fragment.Type').to_pylist()
chg = lib.column('Precursor.Charge').to_numpy(); pmz = lib.column('Precursor.Mz').to_numpy()
seq = lib.column('Modified.Sequence').to_pylist()

relint = np.zeros((P, F), np.float32); prodmz = np.zeros((P, F), np.float32)
frchg  = np.zeros((P, F), np.float32); ordin  = np.zeros((P, F), np.float32)
series = np.zeros((P, F), np.float32); prec   = np.zeros((P, 3), np.float32)
slot = collections.Counter(); ft_seen = {}; hit = set()

for k in range(lib.num_rows):
    name = lid[k]
    if ldec[k] and name.endswith('_decoy'):
        name = name[:-len('_decoy')]
    key = (name, int(ldec[k]))
    i = pos_of.get(key)
    if i is None:
        continue
    j = slot[key]; slot[key] += 1
    if j >= F:
        continue
    hit.add(key)
    relint[i, j] = ri[k]; prodmz[i, j] = pm[k] / 1000.0
    frchg[i, j] = fc[k];  ordin[i, j] = sn[k] / 20.0
    t = str(ft[k])
    if t not in ft_seen: ft_seen[t] = float(len(ft_seen) + 1)
    series[i, j] = ft_seen[t]
    if prec[i, 0] == 0:
        n = sum(1 for c in str(seq[k]) if c.isalpha() and c.isupper())
        prec[i] = (float(chg[k]), float(pmz[k]) / 1000.0, n / 30.0)

s = relint.sum(1, keepdims=True)
relint = np.where(s > 0, relint / np.maximum(s, 1e-9), 0.0).astype(np.float32)

for n_, a in (('relint', relint), ('prodmz', prodmz), ('frcharge', frchg),
              ('ordinal', ordin), ('series', series)):
    np.save(f'{D}/desc_{n_}.npy', a)
np.save(f'{D}/desc_precursor.npy', prec)

t_m = dec == 0; d_m = dec == 1
print(f'matched {len(hit):,} of {P:,} tensor precursors')
for nm, m in (('targets', t_m), ('decoys', d_m)):
    print(f'  {nm:8s} n={m.sum():>7,}  with library intensity {int((relint[m].sum(1)>0).sum()):>7,}'
          f' ({100*(relint[m].sum(1)>0).mean():5.1f}%)   with charge {int((prec[m,0]>0).sum()):>7,}')
print(f'ion series encountered: {ft_seen}')
