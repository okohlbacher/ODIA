#!/usr/bin/env python3
"""D2: is the A->B loss real, or bookkeeping?

Two artefact hypotheses for the 39.9% that produce no candidate:
  H1  the 10% "not in our library" is an ID-STRING mismatch, not a content gap
  H2  the zero-candidate set is the fragment floor (-min_library_fragments 3)
      and the m/z range filters removing them BEFORE the search
Anything left after those two is a genuine detection loss.
"""
import sys, numpy as np, pyarrow.parquet as pq, pyarrow.compute as pc

odia_tsv, dn_pq, lib_pq = sys.argv[1], sys.argv[2], sys.argv[3]

dn = pq.read_table(dn_pq, columns=['Precursor.Id', 'RT', 'Q.Value'])
q = dn.column('Q.Value').to_numpy()
ids = np.array(dn.column('Precursor.Id').to_pylist())
dn_ids = set(ids[q <= 0.01])

lib = pq.read_table(lib_pq, columns=['Precursor.Id', 'Decoy', 'Precursor.Mz',
                                     'Product.Mz', 'Precursor.Charge'])
nfr = pc.list_value_length(lib.column('Product.Mz')).to_numpy()
lid = np.array(lib.column('Precursor.Id').to_pylist())
ldec = lib.column('Decoy').to_numpy()
lmz = lib.column('Precursor.Mz').to_numpy()
lch = lib.column('Precursor.Charge').to_numpy()
tgt = ~ldec.astype(bool)
libmap = {p: (int(n), float(m), int(c))
          for p, n, m, c, d in zip(lid, nfr, lmz, lch, tgt) if d}

missing = sorted(dn_ids - set(libmap))
print(f'H1  DIA-NN confident not found by exact Precursor.Id: {len(missing)}')
print('    examples:', missing[:6])
# strip charge suffix and compare on the modified-sequence part only
def base(p):
    s = str(p)
    return s[:-1] if s and s[-1].isdigit() else s
lib_base = set(base(p) for p in libmap)
soft = sum(1 for p in missing if base(p) in lib_base)
print(f'    of those, sequence present but CHARGE differs: {soft} '
      f'({100*soft/max(len(missing),1):.1f}%)')
# case/format insensitive
lib_norm = set(str(p).replace('(', '[').replace(')', ']').upper() for p in libmap)
fmt = sum(1 for p in missing
          if str(p).replace('(', '[').replace(')', ']').upper() in lib_norm)
print(f'    of those, recovered by bracket/case normalisation: {fmt}')

present = {p: libmap[p] for p in dn_ids if p in libmap}
seen = set()
h = None
with open(odia_tsv) as fh:
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if h is None:
            h = {n: i for i, n in enumerate(f)}
            continue
        if f[h['Decoy']] in ('0', 'false', 'False'):
            seen.add(f[h['Precursor.Id']])
zero = [p for p in present if p not in seen]
print(f'\nH2  in library but zero candidates: {len(zero)} '
      f'of {len(present)} ({100*len(zero)/len(present):.1f}%)')

nz = np.array([present[p][0] for p in zero])
na = np.array([v[0] for v in present.values()])
print(f'    fragments in library   zero-cand: mean {nz.mean():.2f} '
      f'min {nz.min()} <3 {100*(nz<3).mean():.1f}%')
print(f'                           all conf.: mean {na.mean():.2f} '
      f'min {na.min()} <3 {100*(na<3).mean():.1f}%')

mz = np.array([present[p][1] for p in zero])
ma = np.array([v[1] for v in present.values()])
print(f'    precursor m/z          zero-cand: median {np.median(mz):.1f} '
      f'outside 300-1800 {100*((mz<300)|(mz>1800)).mean():.1f}%')
print(f'                           all conf.: median {np.median(ma):.1f} '
      f'outside 300-1800 {100*((ma<300)|(ma>1800)).mean():.1f}%')

cz = np.array([present[p][2] for p in zero])
ca = np.array([v[2] for v in present.values()])
for c in sorted(set(ca.tolist())):
    print(f'    charge {c}: zero-cand {100*(cz==c).mean():5.1f}%   '
          f'all confident {100*(ca==c).mean():5.1f}%')
