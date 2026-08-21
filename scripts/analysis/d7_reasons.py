#!/usr/bin/env python3
"""D7: cross-tabulate the terminal-reason table against DIA-NN's confident set.

This is the measurement that replaces the inference in doc/51. Every library
precursor carries exactly one reason, so "N% of DIA-NN's confident precursors
die at stage X" becomes a count rather than an argument from counter magnitudes.

Two things it still does not settle, stated in the output because they are easy
to forget once the table looks authoritative:

  * DIA-NN's q <= 0.01 calls are NOMINAL, not verified true positives. On this
    library DIA-NN's own entrapment FDP is 7.42%, so ~7% of the denominator is
    false and its reasons are not diagnostic of anything.
  * The key space is (modified sequence + charge, Decoy). Duplicate library rows
    sharing that key are counted once and reported.
"""
import sys, re, collections, numpy as np, pyarrow.parquet as pq

reasons_tsv, dn_pq = sys.argv[1], sys.argv[2]
ALIAS = {'UniMod:4': 'Carbamidomethyl'}
def norm(p):
    return re.sub(r'\((.*?)\)', lambda m: '(' + ALIAS.get(m.group(1), m.group(1)) + ')', str(p))

dn = pq.read_table(dn_pq, columns=['Precursor.Id', 'Q.Value'])
q = dn.column('Q.Value').to_numpy()
ids = [norm(p) for p in dn.column('Precursor.Id').to_pylist()]
conf = {ids[i] for i in np.where(q <= 0.01)[0]}
print(f'DIA-NN confident on this run: {len(conf)}')

reason, dup = {}, 0
allc = collections.Counter()
with open(reasons_tsv) as fh:
    h = fh.readline()
    for line in fh:
        pid, dec, r = line.rstrip('\n').split('\t')
        if dec != '0':
            continue
        allc[r] += 1
        if pid in reason and reason[pid] != r:
            dup += 1
        reason[pid] = r
print(f'target precursors in the table: {len(reason)}'
      + (f'   ({dup} duplicate keys disagreed and the last won)' if dup else ''))

hit = collections.Counter(reason.get(p, '<not in library>') for p in conf)
n = len(conf)
print(f'\n{"terminal reason":>24} {"DIA-NN confident":>17} {"share":>7}'
      f' {"whole library":>14} {"share":>7}')
tot_lib = sum(allc.values()) or 1
for r, c in hit.most_common():
    lib = allc.get(r, 0)
    print(f'{r:>24} {c:>17,} {100*c/n:6.1f}% {lib:>14,} {100*lib/tot_lib:6.1f}%')

good = hit.get('scored', 0)
print(f'\nreached candidate formation: {good:,} of {n:,} ({100*good/n:.1f}%)')
print('the rest, by what stopped them:')
for r, c in hit.most_common():
    if r != 'scored':
        print(f'   {r:<24} {c:>7,}  {100*c/(n-good):5.1f}% of the loss')

print(f'\nCAVEATS')
print(f'  DIA-NN q<=0.01 is nominal; its entrapment FDP on this library is')
print(f'  7.42%, so ~{0.0742*n:.0f} of the {n:,} are false and their reasons mean nothing.')
print(f'  Key space: modified sequence + charge, targets only.')
