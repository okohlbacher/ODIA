#!/usr/bin/env python3
"""D8: diff two terminal-reason tables over DIA-NN's confident set.

A gate change that admits 1,000 more precursors has not necessarily admitted
the ones we were losing -- it can admit a different 1,000 and lose some it used
to keep. Totals cannot see that; a paired diff can.

    d8_reason_diff.py before.tsv after.tsv diann.parquet [before_ids.tsv after_ids.tsv]

With the two rank files (bench.sh writes rank_bench_<tag>.tsv) it also reports
whether the newly admitted precursors were actually IDENTIFIED, which is the
only number that matters -- admission is a means.
"""
import sys, re, collections, numpy as np, pyarrow.parquet as pq

before_tsv, after_tsv, dn_pq = sys.argv[1], sys.argv[2], sys.argv[3]
rank_before = sys.argv[4] if len(sys.argv) > 5 else None
rank_after  = sys.argv[5] if len(sys.argv) > 5 else None

ALIAS = {'UniMod:4': 'Carbamidomethyl'}
def norm(p):
    return re.sub(r'\((.*?)\)', lambda m: '(' + ALIAS.get(m.group(1), m.group(1)) + ')', str(p))

dn = pq.read_table(dn_pq, columns=['Precursor.Id', 'Q.Value'])
q = dn.column('Q.Value').to_numpy()
ids = [norm(p) for p in dn.column('Precursor.Id').to_pylist()]
conf = {ids[i] for i in np.where(q <= 0.01)[0]}

def load(path):
    r = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            pid, dec, reason = line.rstrip('\n').split('\t')
            if dec == '0' and pid in conf:
                r[pid] = reason
    return r

a, b = load(before_tsv), load(after_tsv)
print(f'DIA-NN confident: {len(conf)}   in both tables: {len(set(a) & set(b))}')

moved = collections.Counter()
for p in set(a) | set(b):
    ra, rb = a.get(p, '<absent>'), b.get(p, '<absent>')
    if ra != rb:
        moved[(ra, rb)] += 1
sa = sum(1 for v in a.values() if v == 'scored')
sb = sum(1 for v in b.values() if v == 'scored')
print(f'reached candidate formation: before {sa:,}  after {sb:,}  ({sb-sa:+,})')
print('\nthe net is made of:')
for (ra, rb), c in moved.most_common(12):
    arrow = 'GAINED' if rb == 'scored' else ('LOST  ' if ra == 'scored' else '      ')
    print(f'  {arrow}  {ra:>22} -> {rb:<22} {c:>7,}')

if rank_before and rank_after:
    def acc(path):
        s = set()
        for line in open(path):
            f = line.rstrip('\n').split('\t')
            try:
                if float(f[2]) <= 0.01: s.add(f[0])
            except (ValueError, IndexError):
                pass
        return s
    ia, ib = acc(rank_before), acc(rank_after)
    gained = {p for p in b if b[p] == 'scored' and a.get(p) != 'scored'}
    print(f'\nidentified at q<=0.01 among DIA-NN confident: '
          f'before {len(ia & conf):,}  after {len(ib & conf):,}  '
          f'({len(ib & conf) - len(ia & conf):+,})')
    if gained:
        got = len(gained & ib)
        print(f'of the {len(gained):,} newly ADMITTED, {got:,} were identified '
              f'({100*got/len(gained):.1f}%) -- admission is a means, not the result')
    lost = (ia & conf) - ib
    if lost:
        print(f'and {len(lost):,} that WERE identified no longer are')
