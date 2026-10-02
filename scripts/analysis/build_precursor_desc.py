#!/usr/bin/env python3
"""Precursor-level descriptors the transformer was missing.

It already had per-FRAGMENT library intensity, fragment charge, product m/z,
series and ordinal. It had nothing about the PRECURSOR: not its charge, not its
m/z, not the peptide's length. Precursor charge governs which fragment series
are populated and how crowded the isolation window is, so it belongs beside the
per-fragment descriptors rather than nowhere.

All of these are arm-invariant on the shifted contrast -- the same precursor in
both arms -- so they cannot leak; they can only act through interaction with the
traces, exactly like the fragment descriptors.
"""
import sys, numpy as np, pyarrow.parquet as pq

S = sys.argv[1] if len(sys.argv) > 1 else '/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'
meta = pq.read_table(f'{S}/tensor_ih1_meta.parquet')
ids = meta.column('Precursor.Id').to_pylist()
dec = meta.column('Decoy').to_numpy()
pos_of = {(p, int(d)): i for i, (p, d) in enumerate(zip(ids, dec))}

lib = pq.read_table(f'{S}/corpus_lib.parquet',
                    columns=['Precursor.Id', 'Decoy', 'Precursor.Charge',
                             'Precursor.Mz', 'Modified.Sequence'])
P = len(ids)
out = np.zeros((P, 3), dtype=np.float32)
hit = 0
for p, d, ch, mz, seq in zip(lib.column('Precursor.Id').to_pylist(),
                             lib.column('Decoy').to_numpy(),
                             lib.column('Precursor.Charge').to_numpy(),
                             lib.column('Precursor.Mz').to_numpy(),
                             lib.column('Modified.Sequence').to_pylist()):
    i = pos_of.get((p, int(d)))
    if i is None:
        continue
    hit += 1
    # Stripped length: modification names in parentheses are not residues, and
    # counting them would make a carbamidomethylated peptide look 15 residues
    # longer than the same peptide unmodified.
    n = sum(1 for c in str(seq) if c.isalpha() and c.isupper())
    out[i] = (float(ch), float(mz) / 1000.0, n / 30.0)
np.save(f'{S}/desc_precursor.npy', out)
print(f'{hit:,} of {P:,} rows matched ({100*hit/P:.1f}%)')
print(f'  charge: median {np.median(out[out[:,0]>0,0]):.0f}  '
      f'range {out[out[:,0]>0,0].min():.0f}-{out[:,0].max():.0f}')
print(f'  m/z:    median {1000*np.median(out[out[:,1]>0,1]):.1f}')
print(f'  length: median {30*np.median(out[out[:,2]>0,2]):.0f}')
