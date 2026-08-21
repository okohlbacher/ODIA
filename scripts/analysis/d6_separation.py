#!/usr/bin/env python3
"""D6: the 7,138 that are correctly ranked and still rejected -- why?

Two possibilities with different fixes:

  SEPARATION   their DScore is interleaved with the decoy null. The features do
               not distinguish them. Needs better scoring.
  CALIBRATION  their DScore sits above most decoys and the threshold is
               conservative anyway -- an inflated null, e.g. from the decoy
               argmax over as many candidates as the targets get.

Reported as the decoy percentile each group's DScore lands at. If the rejected
true peaks sit at, say, the 97th percentile of decoys, the null is the problem;
if they sit at the 60th, the features are.
"""
import sys, re, numpy as np, pyarrow.parquet as pq

odia_tsv, dn_pq, lib_pq = sys.argv[1], sys.argv[2], sys.argv[3]
TOL = 20.0
ALIAS = {'UniMod:4': 'Carbamidomethyl'}
def norm(p):
    return re.sub(r'\((.*?)\)', lambda m: '(' + ALIAS.get(m.group(1), m.group(1)) + ')', str(p))

dn = pq.read_table(dn_pq, columns=['Precursor.Id', 'RT', 'Q.Value'])
q = dn.column('Q.Value').to_numpy()
ids = np.array([norm(p) for p in dn.column('Precursor.Id').to_pylist()])
rt = dn.column('RT').to_numpy() * 60.0
conf = {ids[i]: rt[i] for i in np.where(q <= 0.01)[0]}
lib = pq.read_table(lib_pq, columns=['Precursor.Id', 'Decoy'])
L = {p for p, d in zip(lib.column('Precursor.Id').to_pylist(),
                       lib.column('Decoy').to_numpy().astype(bool))
     if not d and p in conf}

best_t, best_d, ncand_t, ncand_d, h = {}, {}, {}, {}, None
with open(odia_tsv) as fh:
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if h is None:
            h = {n: i for i, n in enumerate(f)}; continue
        try:
            ds = float(f[h['DScore']])
        except ValueError:
            continue
        pid = f[h['Precursor.Id']]
        if f[h['Decoy']] in ('0', 'false', 'False'):
            ncand_t[pid] = ncand_t.get(pid, 0) + 1
            if pid in L:
                try:
                    r, qv = float(f[h['RT']]), float(f[h['QValue']])
                except ValueError:
                    continue
                if pid not in best_t or ds > best_t[pid][0]:
                    best_t[pid] = (ds, r, qv)
        else:
            ncand_d[pid] = ncand_d.get(pid, 0) + 1
            if pid not in best_d or ds > best_d[pid]:
                best_d[pid] = ds

dec = np.sort(np.array(list(best_d.values())))
print(f'decoy precursors with a best score: {len(dec)}')

acc, rej = [], []
for p, (ds, r, qv) in best_t.items():
    if abs(r - conf[p]) > TOL:
        continue
    (acc if qv <= 0.01 else rej).append(ds)
acc, rej = np.array(acc), np.array(rej)
print(f'stage D total {len(acc)+len(rej)}: accepted {len(acc)}, rejected {len(rej)}')

def pct(v):
    return 100.0 * np.searchsorted(dec, v) / len(dec)
for name, v in (('accepted at q<=0.01', acc), ('REJECTED', rej)):
    p = pct(v)
    print(f'  {name:22s} DScore median {np.median(v):7.3f}   '
          f'decoy percentile: p25 {np.percentile(p,25):5.1f} '
          f'median {np.median(p):5.1f} p75 {np.percentile(p,75):5.1f}')
print(f'  {"decoy null":22s} DScore median {np.median(dec):7.3f}   '
      f'p99 {np.percentile(dec,99):7.3f}  max {dec.max():7.3f}')
above = (rej[:, None] > np.percentile(dec, [90, 99, 99.9])).mean(axis=0) * 100
print(f'  rejected true peaks above decoy p90/p99/p99.9: '
      f'{above[0]:.1f}% / {above[1]:.1f}% / {above[2]:.1f}%')

ct = np.array([ncand_t.get(p, 0) for p in best_t])
cd = np.array([ncand_d.get(p, 0) for p in best_d])
print(f'\ncandidates per precursor: targets {ct.mean():.3f}  decoys {cd.mean():.3f}'
      f'   (the decoy argmax is taken over as many draws as the target one)')
