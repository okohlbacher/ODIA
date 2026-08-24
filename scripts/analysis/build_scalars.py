#!/usr/bin/env python3
"""The 19 shipped sub-scores as arrays aligned to the trace tensor, per arm.

Unlike the per-fragment descriptors, these are NOT arm-invariant: each arm's
scalars are computed by the production scorer on that arm's own picked
candidate, so they differ between a precursor's own window and its shifted one
and can carry real signal. That also means they are NOT leakage-free the way the
descriptors are, and a model given them is no longer protected by the pairing
argument alone -- the metadata control has to be re-read as a scalars control.

NaN where an arm had no candidate; the model gets a validity mask alongside.
"""
import sys, numpy as np, pyarrow.parquet as pq

S = sys.argv[1] if len(sys.argv) > 1 else '/scratch/kohlbach/odia2x2'
OUT = sys.argv[2] if len(sys.argv) > 2 else '/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'

def best(path):
    out, h, cols = {}, None, None
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

meta = pq.read_table(f'{OUT}/tensor_s08_meta.parquet')
ids = meta.column('Precursor.Id').to_pylist()
dec = meta.column('Decoy').to_numpy()

for tag, path in (('s08', f'{S}/corpus_scores_s08.tsv'),
                  ('s08shift', f'{S}/corpus_scores_s08_shift.tsv')):
    b, cols = best(path)
    A = np.full((len(ids), len(cols)), np.nan, dtype=np.float32)
    ok = np.zeros(len(ids), dtype=bool)
    for i, (p, d) in enumerate(zip(ids, dec)):
        v = b.get((p, str(int(d))))
        if v is not None:
            A[i] = v[1]; ok[i] = True
    # A tree can take NaN; a network cannot. Zero-fill and hand the model an
    # explicit validity flag rather than letting it read the fill as a value.
    A = np.nan_to_num(A, nan=0.0, posinf=0.0, neginf=0.0)
    np.save(f'{OUT}/scal_{tag}.npy', A)
    np.save(f'{OUT}/scal_{tag}_ok.npy', ok)
    print(f'{tag}: {ok.sum():,} of {len(ids):,} rows have scalars ({100*ok.mean():.1f}%)')
np.save(f'{OUT}/scal_cols.npy', np.array(cols))
print(f'{len(cols)} columns: {cols}')
