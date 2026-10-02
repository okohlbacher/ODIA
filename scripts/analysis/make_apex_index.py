#!/usr/bin/env python3
"""Per-arm cycle index of ODIA's OWN picked apex, for each precursor.

Rung 1 scores a candidate the picker already located; rungs 2-4 get a 96-cycle
window and must find the peak themselves. That is a confound in rung 1's favour
and it is the likeliest reason three independent trace models converged just
below it.

Centring both arms on their own picked apex removes it -- and is at the same
time the evaluation doc/56 A1 makes mandatory, because a class-conditioned
centring artefact can only be caught on the deployment distribution, where every
candidate is centred on its own apex. There is no asymmetry here: each arm's
apex comes from the same picker applied to that arm's own window.

No re-extraction: the dump already spans +-90 s, so the apex cycle is
(apex_RT - RT0) / spacing within the window we have.
"""
import sys, numpy as np, pyarrow.parquet as pq

S = sys.argv[1] if len(sys.argv) > 1 else '/scratch/kohlbach/odia2x2'
OUT = sys.argv[2] if len(sys.argv) > 2 else '/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'
SPACING = 1.385

def best_rt(path):
    out, h = {}, None
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if h is None:
            h = {n: i for i, n in enumerate(f)}; continue
        try: ds, rt = float(f[h['DScore']]), float(f[h['RT']])
        except ValueError: continue
        k = (f[h['Precursor.Id']], f[h['Decoy']])
        if k not in out or ds > out[k][0]: out[k] = (ds, rt)
    return out

for tag, scores in (('ih1', f'{S}/corpus_scores_ih1.tsv'),
                    ('ih1shift', f'{S}/corpus_scores_ih1_shift.tsv')):
    b = best_rt(scores)
    m = pq.read_table(f'{S}/tensor_{"ih1" if tag=="ih1" else "ih1shift"}_meta.parquet')
    ids = m.column('Precursor.Id').to_pylist()
    dec = m.column('Decoy').to_numpy()
    rt0 = m.column('RT0').to_numpy()
    idx = np.full(len(ids), -1, dtype=np.int16)
    for i, (p, d) in enumerate(zip(ids, dec)):
        v = b.get((p, str(int(d))))
        if v is None: continue
        c = int(round((v[1] - rt0[i]) / SPACING))
        if 0 <= c < 128: idx[i] = c
    np.save(f'{OUT}/apex_{tag}.npy', idx)
    print(f'{tag}: apex located for {(idx >= 0).sum():,} of {len(idx):,} '
          f'({100*(idx>=0).mean():.1f}%)   median cycle {int(np.median(idx[idx>=0]))}')
