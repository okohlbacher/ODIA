#!/usr/bin/env python3
"""D5: if the gate admitted the precursors it now rejects, what would we GET?

The set we lose is 2.7x fainter than the set we keep, and faint peptides score
worse, so admitting them need not identify them -- they could simply die one
stage later at q <= 0.01 instead. That is the difference between a fix worth
building and a fix worth skipping, and it is predictable WITHOUT running
anything: bin the precursors we DO find by DIA-NN's abundance, read the
acceptance rate off each bin, and apply those rates to the abundance
distribution of the ones we lose.

The prediction is an UPPER bound and is stated as one. Within an abundance bin
the precursors we currently admit are the ones that passed a co-elution gate,
so they are the better-behaved members of that bin; the ones we lose would
score below their bin-mates, not at the bin average.
"""
import sys, re, numpy as np, pyarrow.parquet as pq

odia_tsv, dn_pq, lib_pq = sys.argv[1], sys.argv[2], sys.argv[3]
TOL = 20.0
ALIAS = {'UniMod:4': 'Carbamidomethyl'}
def norm(p):
    return re.sub(r'\((.*?)\)', lambda m: '(' + ALIAS.get(m.group(1), m.group(1)) + ')', str(p))

dn = pq.read_table(dn_pq, columns=['Precursor.Id', 'RT', 'Q.Value', 'Precursor.Quantity'])
q = dn.column('Q.Value').to_numpy()
ids = np.array([norm(p) for p in dn.column('Precursor.Id').to_pylist()])
rt = dn.column('RT').to_numpy() * 60.0
qt = dn.column('Precursor.Quantity').to_numpy()
conf = {ids[i]: (rt[i], qt[i]) for i in np.where(q <= 0.01)[0]}

lib = pq.read_table(lib_pq, columns=['Precursor.Id', 'Decoy'])
ldec = lib.column('Decoy').to_numpy().astype(bool)
L = {p for p, d in zip(lib.column('Precursor.Id').to_pylist(), ldec)
     if not d and p in conf}

best, h = {}, None
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
        if pid not in best or ds > best[pid][0]:
            best[pid] = (ds, r, qv)

found = np.array(sorted(best))
lost  = np.array(sorted(L - set(best)))
fq = np.array([conf[p][1] for p in found])
lq = np.array([conf[p][1] for p in lost])
acc = np.array([1 if (abs(best[p][1] - conf[p][0]) <= TOL and best[p][2] <= 0.01) else 0
                for p in found])

edges = np.percentile(np.concatenate([fq, lq]), [0, 10, 25, 50, 75, 90, 100])
edges[-1] *= 1.0001
print(f'{"DIA-NN quantity":>26}  {"found":>7} {"acc%":>6}  {"lost":>7}  {"predicted":>9}')
tot = 0.0
for i in range(len(edges) - 1):
    lo, hi = edges[i], edges[i + 1]
    mf = (fq >= lo) & (fq < hi)
    ml = (lq >= lo) & (lq < hi)
    rate = acc[mf].mean() if mf.sum() else 0.0
    gain = rate * ml.sum()
    tot += gain
    print(f'{lo:11.0f} .. {hi:11.0f}  {mf.sum():7d} {100*rate:5.1f}%  '
          f'{ml.sum():7d}  {gain:9.0f}')
print(f'\nfound {len(found)}, of which accepted {acc.sum()} ({100*acc.mean():.1f}%)')
print(f'lost  {len(lost)}')
print(f'UPPER-BOUND yield if every lost precursor were admitted and behaved '
      f'like its abundance bin: {tot:.0f} '
      f'(+{100*tot/max(acc.sum(),1):.1f}% on {acc.sum()})')
flat = acc.mean() * len(lost)
print(f'  for contrast, assuming NO abundance dependence at all: {flat:.0f}')
