#!/usr/bin/env python3
"""D3: the staging of D1, with the modification-notation alias, plus a test of
WHY the zero-candidate precursors are lost.

The alias matters: DIA-NN writes C(UniMod:4), our library writes
C(Carbamidomethyl). It is the only modification on either side, so exact
Precursor.Id matching silently dropped EVERY cysteine peptide -- 3,902 of
DIA-NN's 39,149 confident precursors, 10.0% -- from both the numerator and the
denominator of every recovery number computed this way.

For the precursors that yield no candidate at all, two explanations need
separating, because they have opposite fixes:

  "we looked and saw nothing"  -- the gate rejected real but weak signal
  "we never looked there"      -- the ion-mobility band or the RT window was
                                  centred somewhere the peptide is not

DIA-NN reports the OBSERVED RT and 1/K0 for every precursor it identifies, so
the second is directly measurable against what our library predicted.
"""
import sys, re, numpy as np, pyarrow.parquet as pq, pyarrow.compute as pc

odia_tsv, dn_pq, lib_pq = sys.argv[1], sys.argv[2], sys.argv[3]
TOL = float(sys.argv[4]) if len(sys.argv) > 4 else 20.0

ALIAS = {'UniMod:4': 'Carbamidomethyl'}
def norm(p):
    return re.sub(r'\((.*?)\)', lambda m: '(' + ALIAS.get(m.group(1), m.group(1)) + ')', str(p))

dn = pq.read_table(dn_pq, columns=['Precursor.Id', 'RT', 'IM', 'Q.Value',
                                   'Precursor.Mz', 'Precursor.Quantity'])
q = dn.column('Q.Value').to_numpy()
k = q <= 0.01
ids = np.array([norm(p) for p in dn.column('Precursor.Id').to_pylist()])
dn_rt = dn.column('RT').to_numpy() * 60.0
dn_im = dn.column('IM').to_numpy()
dn_mz = dn.column('Precursor.Mz').to_numpy()
dn_qt = dn.column('Precursor.Quantity').to_numpy()
conf = {ids[i]: (dn_rt[i], dn_im[i], dn_mz[i], dn_qt[i]) for i in np.where(k)[0]}
print(f'DIA-NN confident (alias applied): {len(conf)}')

lib = pq.read_table(lib_pq, columns=['Precursor.Id', 'Decoy', 'RT', 'IM',
                                     'Precursor.Mz', 'Product.Mz'])
nfr = pc.list_value_length(lib.column('Product.Mz')).to_numpy()
lid = lib.column('Precursor.Id').to_pylist()
ldec = lib.column('Decoy').to_numpy().astype(bool)
lrt, lim = lib.column('RT').to_numpy(), lib.column('IM').to_numpy()
L = {}
for i, p in enumerate(lid):
    if not ldec[i] and p in conf:
        L[p] = (lrt[i], lim[i], int(nfr[i]))
print(f'A  in our library: {len(L)}  ({100*len(L)/len(conf):.1f}%)')

seen, best, near = set(), {}, {}
h = None
with open(odia_tsv) as fh:
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if h is None:
            h = {n: i for i, n in enumerate(f)}; continue
        pid = f[h['Precursor.Id']]
        if pid not in L or f[h['Decoy']] not in ('0', 'false', 'False'):
            continue
        try:
            ds, r, qv = float(f[h['DScore']]), float(f[h['RT']]), float(f[h['QValue']])
        except ValueError:
            continue
        seen.add(pid)
        near.setdefault(pid, []).append((ds, r, qv))
        if pid not in best or ds > best[pid][0]:
            best[pid] = (ds, r, qv)

C = D = E = 0
for pid, rows in near.items():
    t = conf[pid][0]
    if not any(abs(x[1] - t) <= TOL for x in rows):
        continue
    C += 1
    b = best[pid]
    if abs(b[1] - t) <= TOL:
        D += 1
        if b[2] <= 0.01:
            E += 1
n = len(L)
print(f'B  candidate proposed:  {len(seen):7d}  ({100*len(seen)/n:5.1f}% of A)')
print(f'C  true peak available: {C:7d}  ({100*C/n:5.1f}% of A)   [+-{TOL:.0f}s]')
print(f'D  true peak best:      {D:7d}  ({100*D/max(C,1):5.1f}% of C)')
print(f'E  accepted q<=0.01:    {E:7d}  ({100*E/max(D,1):5.1f}% of D, {100*E/n:5.1f}% of A)')

zero = [p for p in L if p not in seen]
found = [p for p in L if p in seen]
print(f'\n-- the {len(zero)} that yield NO candidate, against the {len(found)} that do --')
def stat(name, fz, ff, fmt='{:8.3f}'):
    z = np.array([fz(p) for p in zero]); g = np.array([ff(p) for p in found])
    print(f'  {name:26s} zero ' + fmt.format(np.median(z)) +
          '   found ' + fmt.format(np.median(g)))
    return z, g

dz, dg = stat('|dIM| library vs DIA-NN', lambda p: abs(L[p][1] - conf[p][1]),
                                          lambda p: abs(L[p][1] - conf[p][1]))
print(f'    beyond the +-0.025 band:   zero {100*(dz>0.025).mean():5.1f}%   '
      f'found {100*(dg>0.025).mean():5.1f}%')
stat('DIA-NN 1/K0',       lambda p: conf[p][1], lambda p: conf[p][1])
stat('precursor m/z',     lambda p: conf[p][2], lambda p: conf[p][2])
stat('DIA-NN quantity',   lambda p: conf[p][3], lambda p: conf[p][3], '{:12.0f}')
stat('fragments in lib',  lambda p: L[p][2],    lambda p: L[p][2], '{:8.1f}')
qz = np.array([conf[p][3] for p in zero]); qg = np.array([conf[p][3] for p in found])
lo = np.percentile(qg, 10)
print(f'    below the found set\'s 10th pct quantity ({lo:.0f}): '
      f'zero {100*(qz<lo).mean():.1f}%   found 10.0%')
