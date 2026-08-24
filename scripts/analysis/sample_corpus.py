#!/usr/bin/env python3
"""Select the precursors whose traces form the training corpus, and write a
library subset containing exactly them.

Four groups, deliberately a SUPERSET of any single label design, so that a
redesign after review does not require re-extracting anything:

  pos_diann   targets DIA-NN identifies at q <= 0.01           (positives)
  ent         ENTRAPMENT targets -- Arabidopsis accessions in  (negatives)
              the searched library: real sequences, real
              fragment masses, absent from human liver
  unid        ordinary human targets DIA-NN does NOT identify  (negatives?)
              -- absent is UNVERIFIABLE for these, which is
              why they are kept separate rather than pooled
  decoy       decoys, HELD OUT of training entirely, as the
              check that a presence model ranks them low

Entrapment negatives are a different ORGANISM, so composition can separate the
classes without any trace information. Sampling is therefore STRATIFIED on
(charge, precursor m/z, peptide length, fragment count) to match the positive
distribution -- and the strata are written out so the match can be audited
rather than asserted.
"""
import sys, re, collections, numpy as np, pyarrow.parquet as pq, pyarrow as pa
import pyarrow.compute as pc

lib_pq, dn_pq, arab_txt, out_ids, out_lib = sys.argv[1:6]
N_ENT   = int(sys.argv[6]) if len(sys.argv) > 6 else 20000
N_UNID  = int(sys.argv[7]) if len(sys.argv) > 7 else 12000
N_DECOY = int(sys.argv[8]) if len(sys.argv) > 8 else 12000
rng = np.random.default_rng(20260824)

ALIAS = {'UniMod:4': 'Carbamidomethyl'}
def norm(p):
    return re.sub(r'\((.*?)\)', lambda m: '(' + ALIAS.get(m.group(1), m.group(1)) + ')', str(p))

dn = pq.read_table(dn_pq, columns=['Precursor.Id', 'Q.Value'])
q = dn.column('Q.Value').to_numpy()
conf = {norm(p) for p, k in zip(dn.column('Precursor.Id').to_pylist(), q <= 0.01) if k}
print(f'DIA-NN confident: {len(conf):,}')

arab = {l.strip() for l in open(arab_txt) if l.strip()}
lib = pq.read_table(lib_pq, columns=['Precursor.Id', 'Decoy', 'Precursor.Charge',
                                     'Precursor.Mz', 'Protein.Group', 'Product.Mz',
                                     'Modified.Sequence'])
nfr = pc.list_value_length(lib.column('Product.Mz')).to_numpy()
ids  = lib.column('Precursor.Id').to_pylist()
dec  = lib.column('Decoy').to_numpy().astype(bool)
chg  = lib.column('Precursor.Charge').to_numpy()
mz   = lib.column('Precursor.Mz').to_numpy()
pg   = lib.column('Protein.Group').to_pylist()
seq  = lib.column('Modified.Sequence').to_pylist()

def strat(i):
    """Coarse cell the composition match is made within."""
    return (int(chg[i]), int(mz[i] // 50), min(len(str(seq[i])) // 4, 8),
            min(int(nfr[i]), 14))

pos, ent, unid, dcy = [], [], [], []
for i, pid in enumerate(ids):
    if dec[i]:
        dcy.append(i); continue
    is_ent = any(a in arab for a in str(pg[i]).replace(';', '|').split('|') if a)
    if is_ent:
        ent.append(i)
    elif pid in conf:
        pos.append(i)
    else:
        unid.append(i)
print(f'library targets: pos {len(pos):,}  entrapment {len(ent):,}  '
      f'unidentified {len(unid):,}  decoys {len(dcy):,}')

# Stratified match: draw negatives cell-by-cell in the positives' proportions,
# so amino-acid composition, length, charge and m/z cannot separate the classes.
want = collections.Counter(strat(i) for i in pos)
tot  = sum(want.values())
def matched(pool, n, name):
    by = collections.defaultdict(list)
    for i in pool: by[strat(i)].append(i)
    out, short = [], 0
    for cell, c in want.items():
        k = int(round(n * c / tot))
        have = by.get(cell, [])
        if len(have) <= k: out.extend(have); short += k - len(have)
        else: out.extend(rng.choice(have, size=k, replace=False).tolist())
    print(f'  {name}: {len(out):,} drawn, {short:,} short of the matched target '
          f'({100*short/max(n,1):.1f}% of cells underpopulated)')
    return out

print('composition-matched sampling:')
ent_s  = matched(ent,  N_ENT,   'entrapment')
unid_s = matched(unid, N_UNID,  'unidentified')
dcy_s  = rng.choice(dcy, size=min(N_DECOY, len(dcy)), replace=False).tolist()
print(f'  decoys: {len(dcy_s):,} drawn UNMATCHED (held out of training; they are '
      f'the check, not a class)')

keep = sorted(set(pos) | set(ent_s) | set(unid_s) | set(dcy_s))
group = {}
for i in pos: group[i] = 'pos_diann'
for i in ent_s: group[i] = 'ent'
for i in unid_s: group[i] = 'unid'
for i in dcy_s: group[i] = 'decoy'

with open(out_ids, 'w') as f:
    f.write('Precursor.Id\tDecoy\tGroup\tPrecursor.Mz\tPrecursor.Charge\tN.Fragments\n')
    for i in keep:
        f.write(f'{ids[i]}\t{1 if dec[i] else 0}\t{group[i]}\t{mz[i]:.4f}\t'
                f'{int(chg[i])}\t{int(nfr[i])}\n')
print(f'\n{len(keep):,} precursors -> {out_ids}')

sub = pq.read_table(lib_pq).take(pa.array(keep))
pq.write_table(sub, out_lib)
print(f'library subset -> {out_lib}  ({sub.num_rows:,} rows)')
