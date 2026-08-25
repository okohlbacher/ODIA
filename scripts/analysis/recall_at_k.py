#!/usr/bin/env python3
"""Does ODIA EMIT a candidate at the right retention time, and merely rank it
below a wrong one -- or does it never emit one at all?

This decides what the boundary work can be worth. Sub-scores only reorder a
list; if the correct position is not on the list, no scoring change can
recover it and the deficit is the picker's. The earlier replay measured only
the TOP-ranked candidate and found ODIA's best apex more than 3 cycles from
DIA-NN's for 37% of confident precursors. That number cannot distinguish the
two cases, and they need completely different work.
"""
import sys, re, numpy as np, pyarrow.parquet as pq
from collections import defaultdict
S = sys.argv[1] if len(sys.argv) > 1 else '/scratch/kohlbach/odia2x2'
SP = 1.385
ALIAS = {'UniMod:4': 'Carbamidomethyl'}
nrm = lambda p: re.sub(r'\((.*?)\)', lambda m: '(' + ALIAS.get(m.group(1), m.group(1)) + ')', str(p))

dn = pq.read_table('/scratch/kohlbach/xic/dn_xic.parquet', columns=['Precursor.Id','RT','Q.Value'])
q = dn.column('Q.Value').to_numpy()
TRT = {nrm(p): r*60.0 for p, r, k in zip(dn.column('Precursor.Id').to_pylist(),
                                         dn.column('RT').to_numpy(), q <= 0.01) if k}
cands = defaultdict(list)
h = None
for line in open(f'{S}/corpus_scores_s08.tsv'):
    f = line.rstrip('\n').split('\t')
    if h is None:
        h = {n: i for i, n in enumerate(f)}; continue
    if f[h['Decoy']] != '0': continue
    try: ds, rt = float(f[h['DScore']]), float(f[h['RT']])
    except ValueError: continue
    cands[f[h['Precursor.Id']]].append((ds, rt))

TOL = 3 * SP   # the same +-3 cycles the boundary replay used
NEAR = []      # signed distance of the CLOSEST candidate, in cycles
n_tot = n_any = 0
rank_of_correct, n_cand = [], []
top1 = 0
for pid, cs in cands.items():
    t = TRT.get(pid)
    if t is None: continue
    n_tot += 1
    cs.sort(key=lambda x: -x[0])
    n_cand.append(len(cs))
    d = [(rt - t) / SP for _, rt in cs]
    NEAR.append(min(d, key=abs))
    hit = [i for i, (_, rt) in enumerate(cs) if abs(rt - t) <= TOL]
    if hit:
        n_any += 1
        rank_of_correct.append(hit[0] + 1)
        if hit[0] == 0: top1 += 1

print(f'{n_tot:,} DIA-NN-confident precursors that ODIA emitted candidates for')
print(f'  candidates per precursor: mean {np.mean(n_cand):.2f}  p50 {np.percentile(n_cand,50):.0f}'
      f'  p90 {np.percentile(n_cand,90):.0f}  max {max(n_cand)}\n')
print(f'  a candidate exists within +-3 cycles of DIA-NN (recall@all): '
      f'{100*n_any/n_tot:5.1f}%   ({n_any:,})')
print(f'  and it is ranked FIRST (top-1):                             '
      f'{100*top1/n_tot:5.1f}%   ({top1:,})')
r = np.array(rank_of_correct)
for k in (1, 2, 3, 5):
    print(f'  recall@{k}: {100*(r<=k).sum()/n_tot:5.1f}%')
nd = np.array(NEAR)
print(f'\n  CONTROL -- is +-3 cycles hiding a systematic RT offset?')
print(f'    signed distance of the CLOSEST candidate, cycles: '
      f'median {np.median(nd):+.2f}  mean {nd.mean():+.2f}')
print(f'    |closest| within: 1 cy {100*(np.abs(nd)<=1).mean():.1f}%   '
      f'3 cy {100*(np.abs(nd)<=3).mean():.1f}%   5 cy {100*(np.abs(nd)<=5).mean():.1f}%   '
      f'10 cy {100*(np.abs(nd)<=10).mean():.1f}%   20 cy {100*(np.abs(nd)<=20).mean():.1f}%')
print(f'    so widening the tolerance to 10 cycles adds only '
      f'{100*((np.abs(nd)<=10).mean()-(np.abs(nd)<=3).mean()):.1f} points: the misses are FAR, '
      f'not just off by calibration.')
print(f'\n  Of the precursors that DO have a correct candidate, it ranks first '
      f'{100*(r==1).mean():.1f}% of the time.')
print(f'  So the ceiling that better sub-scores could reach is {100*n_any/n_tot:.1f}%,')
print(f'  and {100*(n_any-top1)/n_tot:.1f} points of that are currently lost to RANKING,')
print(f'  {100*(n_tot-n_any)/n_tot:.1f} points to the candidate never being emitted.')
