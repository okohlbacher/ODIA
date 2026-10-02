#!/usr/bin/env python3
"""ODIA's own DScore per corpus row, so the transfer comparison ranks against
what the tool actually ships rather than against a single-column proxy."""
import sys, numpy as np, pyarrow.parquet as pq
S = sys.argv[1] if len(sys.argv) > 1 else '/scratch/kohlbach/odia2x2'
D = sys.argv[2] if len(sys.argv) > 2 else '/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'
best, h = {}, None
for line in open(f'{S}/corpus_scores_ih1.tsv'):
    f = line.rstrip('\n').split('\t')
    if h is None:
        h = {n: i for i, n in enumerate(f)}; continue
    try: ds, q = float(f[h['DScore']]), float(f[h['QValue']])
    except ValueError: continue
    k = (f[h['Precursor.Id']], f[h['Decoy']])
    if k not in best or ds > best[k][0]: best[k] = (ds, q)
meta = pq.read_table(f'{D}/tensor_ih1_meta.parquet')
ids = meta.column('Precursor.Id').to_pylist(); dec = meta.column('Decoy').to_numpy()
ds = np.full(len(ids), -1e9, dtype=np.float32); ok = np.zeros(len(ids), dtype=bool)
for i, (p, d) in enumerate(zip(ids, dec)):
    v = best.get((p, str(int(d))))
    if v is not None: ds[i], ok[i] = v[0], True
np.save(f'{D}/dscore_ih1.npy', ds); np.save(f'{D}/dscore_ih1_ok.npy', ok)
print(f'{ok.sum():,} of {len(ids):,} rows carry a DScore ({100*ok.mean():.1f}%)')
print(f'  range {ds[ok].min():.2f} .. {ds[ok].max():.2f}  median {np.median(ds[ok]):.2f}')
