#!/usr/bin/env python3
"""D1: where do DIA-NN-confident precursors die inside ODIA?

Five stages, each a strict subset of the one above:

  A  in library                 -- the precursor exists in the library we search
  B  candidate proposed         -- ODIA emitted at least one row for it
  C  true peak available        -- some candidate sits within tol of DIA-NN's RT
  D  true peak ranked best      -- that candidate has the highest DScore
  E  accepted                   -- its q <= 0.01

A->B is DETECTION.  C->D is RANKING.  D->E is the FDR threshold.
The B->C step separates "we found peaks but not the right one" from "we found
nothing", which no earlier analysis has separated on the FULL library.
"""
import sys, numpy as np, pyarrow.parquet as pq

odia_tsv = sys.argv[1]
dn_pq    = sys.argv[2]
lib_pq   = sys.argv[3]
tols     = [float(x) for x in (sys.argv[4:] or [10, 20, 30, 60])]

dn = pq.read_table(dn_pq, columns=['Precursor.Id', 'RT', 'Q.Value'])
q  = dn.column('Q.Value').to_numpy()
ids = np.array(dn.column('Precursor.Id').to_pylist())
rt  = dn.column('RT').to_numpy() * 60.0          # DIA-NN reports minutes
keep = q <= 0.01
dn_rt = dict(zip(ids[keep], rt[keep]))
print(f'DIA-NN confident: {len(dn_rt)}')

lib = set(pq.read_table(lib_pq, columns=['Precursor.Id']).column('Precursor.Id').to_pylist())
in_lib = {p: r for p, r in dn_rt.items() if p in lib}
print(f'  of which in our library: {len(in_lib)}  '
      f'({100*len(in_lib)/len(dn_rt):.1f}%)')

# one pass over the ODIA TSV, keeping only rows for precursors we care about
best = {}     # id -> (dscore, rt, qvalue) of the top-DScore candidate
near = {}     # id -> list of (dscore, rt, qvalue) rows
h = None
with open(odia_tsv) as fh:
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if h is None:
            h = {n: i for i, n in enumerate(f)}
            continue
        pid = f[h['Precursor.Id']]
        if pid not in in_lib:
            continue
        d = f[h['Decoy']]
        if d not in ('0', 'false', 'False'):
            continue
        try:
            ds, r, qv = float(f[h['DScore']]), float(f[h['RT']]), float(f[h['QValue']])
        except ValueError:
            continue
        near.setdefault(pid, []).append((ds, r, qv))
        if pid not in best or ds > best[pid][0]:
            best[pid] = (ds, r, qv)

print(f'\nB  candidate proposed:      {len(near):7d}  '
      f'({100*len(near)/len(in_lib):5.1f}% of A)')
mult = np.array([len(v) for v in near.values()])
print(f'   candidates per precursor: mean {mult.mean():.3f}  '
      f'at cap(3) {100*(mult>=3).mean():.1f}%  =1 {100*(mult==1).mean():.1f}%')

print(f'\n{"tol":>5} {"C avail":>9} {"%A":>6} {"D best":>9} {"%C":>6} '
      f'{"E q<=.01":>9} {"%D":>6} {"%A":>6}')
for tol in tols:
    C = D = E = 0
    for pid, rows in near.items():
        t = in_lib[pid]
        hit = [x for x in rows if abs(x[1] - t) <= tol]
        if not hit:
            continue
        C += 1
        b = best[pid]
        if abs(b[1] - t) <= tol:
            D += 1
            if b[2] <= 0.01:
                E += 1
    print(f'{tol:5.0f} {C:9d} {100*C/len(in_lib):5.1f}% {D:9d} '
          f'{100*D/max(C,1):5.1f}% {E:9d} {100*E/max(D,1):5.1f}% '
          f'{100*E/len(in_lib):5.1f}%')

# how far away IS the nearest candidate, for the ones that miss at 60 s?
res = []
for pid, rows in near.items():
    t = in_lib[pid]
    res.append(min(abs(x[1] - t) for x in rows))
res = np.array(res)
print(f'\nnearest-candidate |dRT| over all {len(res)} proposed: '
      f'median {np.median(res):.1f}s  p75 {np.percentile(res,75):.1f}s  '
      f'p90 {np.percentile(res,90):.1f}s  p95 {np.percentile(res,95):.1f}s')
