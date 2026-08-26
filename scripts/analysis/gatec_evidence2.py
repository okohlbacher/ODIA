#!/usr/bin/env python3
"""Gate C is the bottleneck, not the picker. Can its statistic be made better
WITHOUT being made looser?

Measured: of DIA-NN's confident precursors that ODIA loses, gate_c accounts for
98.2% and the picker returning nothing accounts for 0.6%. Relaxing gate_c has
already been measured end-to-end at -20.4% at matched entrapment FDP, so the
threshold is not the lever. The statistic might be.

coelutionEvidence today gives every transition an EQUAL vote: sqrt-transform,
robust z-score against the transition's own median/MAD, sum across transitions,
smooth, take the max over positions. Library intensity is never consulted -- so
a fragment the library says should be the brightest in the spectrum counts
exactly as much as one it says is barely there.

Arms, all at MATCHED ADMITTED VOLUME so that "better" cannot mean "more":
  equal      what ships today
  libw       each transition's z-score weighted by its library intensity
  libw_sqrt  weighted by sqrt(intensity), a softer version
  topn       only the library's top-6 fragments vote at all (DIA-NN's TopF)

The winner's-curse objection that killed library weighting at candidate
SELECTION does not transfer: that was a coupling with LIBRARY_CORR, a downstream
classifier feature. Gate C's statistic is not a feature, and its max over
positions applies to targets and decoys alike.
"""
import numpy as np, pyarrow.parquet as pq
D='/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'
HALF=2   # gate_smooth_half default

def evidence(tr, msk, rel, mode):
    F,T = tr.shape
    s = np.zeros(T)
    contributing = 0
    for f in range(F):
        if not msk[f]: continue
        y = np.sqrt(np.clip(tr[f], 0, None))
        med = np.median(y)
        mad = np.median(np.abs(y-med))
        if not (mad > 0): continue
        w = 1.0
        if mode == 'libw':      w = rel[f]
        elif mode == 'libw_sqrt': w = np.sqrt(max(rel[f], 0.0))
        elif mode == 'topn':
            w = 1.0   # membership handled by the caller's mask
        if w <= 0: continue
        contributing += 1
        s += (y-med) * (1.0/(1.4826*mad)) * w
    if contributing == 0: return 0.0, 0
    k = 2*HALF+1
    c = np.convolve(s, np.ones(k), mode='same')/k
    return float(c.max()), contributing

lab=pq.read_table(f'{D}/tensor_s08_labels.parquet')
label=np.array(lab.column('Label').to_pylist())
meta=pq.read_table(f'{D}/tensor_s08_meta.parquet')
dec=np.array(meta.column('Decoy').to_pylist())
X=np.load(f'{D}/tensor_s08_traces.npy',mmap_mode='r'); M=np.load(f'{D}/tensor_s08_mask.npy')
REL=np.load(f'{D}/desc_relint.npy')

rng=np.random.default_rng(0)
pos=rng.permutation(np.flatnonzero(label=='pos'))[:9000]
dcy=rng.permutation(np.flatnonzero(dec==1))[:9000]
sel=np.concatenate([pos,dcy])
istgt=np.concatenate([np.ones(len(pos),bool),np.zeros(len(dcy),bool)])

MODES=['equal','libw','libw_sqrt','topn']
S={m:[] for m in MODES}; keep=[]
for a0 in range(0,len(sel),200):
    idx=sel[a0:a0+200]; T=np.asarray(X[idx],dtype=np.float32)
    for j in range(len(idx)):
        i=idx[j]; msk=M[i].copy(); tr=T[j].astype(float); rel=REL[i].astype(float)
        if msk.sum()<2: continue
        keep.append(istgt[a0+j])
        for m in MODES:
            mm = msk.copy()
            if m=='topn':
                order=np.argsort(-rel)
                top=set(order[:6].tolist())
                mm = np.array([msk[f] and f in top for f in range(len(msk))])
                if mm.sum()<2: mm = msk.copy()
            v,_=evidence(tr,mm,rel,m)
            S[m].append(v)
keep=np.array(keep)
print(f'{len(keep):,} precursors ({keep.sum():,} DIA-NN-confident targets, '
      f'{(~keep).sum():,} decoys)\n')
print('Matched admitted VOLUME: the threshold for each arm is set so the same')
print('total number of precursors is admitted. Better therefore cannot mean more.\n')
PCTS=[20,40,60,70,80,90]
print(f'{"statistic":>12}' + ''.join(f'{"adm "+str(p)+"%":>12}' for p in PCTS))
for m in MODES:
    v=np.array(S[m]); row=f'{m:>12}'
    for p in PCTS:
        thr=np.percentile(v,100-p)
        adm=v>=thr
        tgt=keep[adm].sum(); dcyn=(~keep[adm]).sum()
        row+=f'{100*tgt/max(keep.sum(),1):>11.1f}%'
    print(row)
print()
print(f'{"statistic":>12}' + ''.join(f'{"T/D "+str(p)+"%":>12}' for p in PCTS))
for m in MODES:
    v=np.array(S[m]); row=f'{m:>12}'
    for p in PCTS:
        thr=np.percentile(v,100-p); adm=v>=thr
        t=keep[adm].sum(); d=(~keep[adm]).sum()
        row+=f'{(t/max(d,1)):>12.2f}'
    print(row)
print('\n  Left table: what fraction of DIA-NN-confident targets survive the gate.')
print('  Right table: target/decoy ratio among the admitted -- the reference-free half.')
