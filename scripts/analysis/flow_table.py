#!/usr/bin/env python3
"""One flow table, because the emission numbers came from different denominators
and a reviewer was right that they cannot be read together without it.

Reconciles: DIA-NN's headline count, the confident set used as the ruler, what
ODIA searched, what it emitted, and what it emitted in the right place.
"""
import re, numpy as np, pyarrow.parquet as pq
from collections import defaultdict
S='/scratch/kohlbach/odia2x2'; SP=1.385
AL={'UniMod:4':'Carbamidomethyl'}
nrm=lambda p: re.sub(r'\((.*?)\)',lambda m:'('+AL.get(m.group(1),m.group(1))+')',str(p))
dn=pq.read_table('/scratch/kohlbach/xic/dn_xic.parquet',columns=['Precursor.Id','RT','Q.Value'])
q=dn.column('Q.Value').to_numpy(); pid=dn.column('Precursor.Id').to_pylist(); rt=dn.column('RT').to_numpy()
allp=set(nrm(p) for p in pid)
TRT={nrm(p):r*60.0 for p,r,k in zip(pid,rt,q<=0.01) if k}
cands=defaultdict(list); h=None; seen=set()
for line in open(f'{S}/corpus_scores_ih1.tsv'):
    f=line.rstrip('\n').split('\t')
    if h is None: h={n:i for i,n in enumerate(f)}; continue
    if f[h['Decoy']]!='0': continue
    seen.add(f[h['Precursor.Id']])
    try: ds,r=float(f[h['DScore']]),float(f[h['RT']])
    except ValueError: continue
    cands[f[h['Precursor.Id']]].append((ds,r))
TOL=3*SP
both=[p for p in TRT if p in cands]
near=[p for p in both if any(abs(r-TRT[p])<=TOL for _,r in cands[p])]
top1=[p for p in both if max(cands[p],key=lambda x:x[0])[1:] and abs(max(cands[p],key=lambda x:x[0])[1]-TRT[p])<=TOL]
rows=[
 ('rows in the DIA-NN report (any q)',            len(allp)),
 ('DIA-NN confident, q <= 0.01 -- THE RULER',      len(TRT)),
 ('...that ODIA also searched and scored',         len(both)),
 ('...with any candidate within +-3 cycles',       len(near)),
 ('...and that candidate ranked first',            len(top1)),
]
print(f'{"stage":<44}{"n":>9}{"% of ruler":>12}{"% of prev":>11}')
prev=None
for nm,n in rows:
    pr = f'{100*n/prev:>10.1f}%' if prev else f'{"":>11}'
    print(f'  {nm:<42}{n:>9,}{100*n/len(TRT):>11.1f}%{pr}')
    prev=n
print(f'\n  precursors ODIA emitted ANY candidate for (targets): {len(seen):,}')
print(f'  of the confident ruler, ODIA scored {100*len(both)/len(TRT):.1f}%; the rest were')
print(f'  never searched (not in the library subset) or yielded no candidate at all.')
print(f'\n  The 33,330 headline is DIA-NN on the .d file with the production library.')
print(f'  This table is the corpus run: a different library subset, so the ruler here')
print(f'  is {len(TRT):,}, and none of these percentages should be compared to it.')
