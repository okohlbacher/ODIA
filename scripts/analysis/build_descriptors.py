#!/usr/bin/env python3
"""Per-fragment DESCRIPTORS aligned to the trace tensor.

All three reviewers found the same defect in doc/58: the trace models never saw
which fragments SHOULD be strong. Rung 1's var_library_corr, var_library_dotprod
and var_library_rmsd compare observed against EXPECTED, and an expected-intensity
vector has AUC 0.5 alone while being highly predictive through its relationship
with the observed traces. So the comparison was scalars-with-a-library-prior
against anonymous traces, and the negative it produced is not licensed.

Emits, per (precursor, fragment slot): library relative intensity, product m/z,
fragment charge, ion series and ordinal. The tensor's slot is
`transition_index - transition_begin`, which is the library's own order, so the
list columns index straight across.
"""
import sys, numpy as np, pyarrow.parquet as pq

S = sys.argv[1] if len(sys.argv) > 1 else '/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'
F = 12
meta = pq.read_table(f'{S}/tensor_ih1_meta.parquet')
ids = meta.column('Precursor.Id').to_pylist()
dec = meta.column('Decoy').to_numpy()
pos_of = {(p, int(d)): i for i, (p, d) in enumerate(zip(ids, dec))}

lib = pq.read_table(f'{S}/corpus_lib.parquet',
                    columns=['Precursor.Id', 'Decoy', 'Relative.Intensity', 'Product.Mz',
                             'Fragment.Charge', 'Fragment.Series.Number', 'Fragment.Type'])
P = len(ids)
relint = np.zeros((P, F), dtype=np.float32)
prodmz = np.zeros((P, F), dtype=np.float32)
frchg  = np.zeros((P, F), dtype=np.float32)
ordin  = np.zeros((P, F), dtype=np.float32)
series = np.zeros((P, F), dtype=np.float32)
hit = 0
ft_seen = {}
for p, d, ri, pm, fc, sn, ft in zip(lib.column('Precursor.Id').to_pylist(),
                                    lib.column('Decoy').to_numpy(),
                                    lib.column('Relative.Intensity').to_pylist(),
                                    lib.column('Product.Mz').to_pylist(),
                                    lib.column('Fragment.Charge').to_pylist(),
                                    lib.column('Fragment.Series.Number').to_pylist(),
                                    lib.column('Fragment.Type').to_pylist()):
    i = pos_of.get((p, int(d)))
    if i is None:
        continue
    hit += 1
    n = min(F, len(ri))
    relint[i, :n] = np.asarray(ri[:n], dtype=np.float32)
    prodmz[i, :n] = np.asarray(pm[:n], dtype=np.float32)
    frchg[i, :n]  = np.asarray(fc[:n], dtype=np.float32)
    ordin[i, :n]  = np.asarray(sn[:n], dtype=np.float32)
    for k in range(n):
        t = str(ft[k])
        if t not in ft_seen: ft_seen[t] = float(len(ft_seen) + 1)
        series[i, k] = ft_seen[t]

# Library intensities are relative and arbitrarily scaled; the model needs the
# SHAPE of the expectation, so normalise each precursor's vector to sum 1. The
# absolute scale carries nothing a library can be trusted on.
s = relint.sum(1, keepdims=True)
relint_n = np.where(s > 0, relint / np.maximum(s, 1e-9), 0.0).astype(np.float32)

np.save(f'{S}/desc_relint.npy', relint_n)
np.save(f'{S}/desc_prodmz.npy', prodmz)
np.save(f'{S}/desc_frcharge.npy', frchg)
np.save(f'{S}/desc_ordinal.npy', ordin)
np.save(f'{S}/desc_series.npy', series)
print(f'{hit:,} of {P:,} tensor rows matched a library entry '
      f'({100*hit/P:.1f}%; the rest are ODIA-regenerated decoys, unused in training)')
print(f'ion series encountered: {ft_seen}')
print(f'library intensity: nonzero slots {100*(relint_n>0).mean():.1f}%, '
      f'median top share {np.median(relint_n.max(1)[relint_n.max(1)>0]):.3f}')
