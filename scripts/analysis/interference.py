#!/usr/bin/env python3
"""DIA-NN's fragment-level interference correction, measured before it is built.

The idea: within one candidate, the fragment least affected by interference
(highest summed correlation to the others among the library's brightest) carries
the TRUE elution shape. Every other fragment is then that shape plus whatever
else happens to be co-eluting in its m/z channel. Project each fragment onto the
reference profile and the projection is the part that belongs to this peptide;
the residual is interference.

Two things that could be worth having, measured separately:
  explained   the fraction of each fragment's variance the common profile
              explains, averaged -- a NEW per-candidate interference score
  lib_clean   library correlation computed on the PROJECTED areas instead of
              the raw ones -- the existing score with interference removed

Judged reference-free: target fraction among the top-N with decoys as control.
This is label-symmetric by construction -- a per-candidate transform applied
identically to both classes -- unlike DIA-NN's cross-precursor remove_ifs, which
deletes targets only and must not be ported symmetrically.
"""
import numpy as np, pyarrow.parquet as pq
D='/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'
HALF=2
def scores(tr, msk, rel, k):
    n=tr.shape[1]; lo=max(k-HALF,0); hi=min(k+HALF,n-1)
    seg=tr[:,lo:hi+1]
    fl=np.concatenate([tr[:,max(lo-20,0):lo],tr[:,hi+1:hi+21]],axis=1)
    base=np.median(fl,axis=1) if fl.shape[1] else np.zeros(tr.shape[0])
    seg=np.clip(seg-base[:,None],0,None)
    f=msk&(rel>0)
    if f.sum()<3: return None
    idx=np.flatnonzero(f)
    # reference = highest summed correlation among the library's brightest six
    bright=idx[np.argsort(-rel[idx])[:6]]
    S=seg[bright]; Sc=S-S.mean(1,keepdims=True); nn=np.linalg.norm(Sc,axis=1)
    ok=nn>0
    if ok.sum()<3: return None
    U=Sc[ok]/nn[ok][:,None]; C=U@U.T; np.fill_diagonal(C,0.0)
    ref=bright[ok][int(np.argmax(C.sum(1)))]
    y=seg[ref].astype(float); yy=float(y@y)
    if yy<=0: return None
    raw=seg[idx].sum(1)
    proj=np.array([max(0.0,float(seg[j]@y)/yy) for j in idx])*float(y.sum())
    tot=np.array([float(seg[j]@seg[j]) for j in idx])
    expl=np.array([ (float(seg[j]@y)**2/yy)/t if t>0 else 0.0 for j,t in zip(idx,tot) ])
    l=rel[idx]
    def corr(a):
        if a.std()<=0 or l.std()<=0: return -1.0
        r=float(np.corrcoef(a,l)[0,1]); return r if np.isfinite(r) else -1.0
    return corr(raw), corr(proj), float(np.mean(np.clip(expl,0,1)))

lab=pq.read_table(f'{D}/tensor_ih1_labels.parquet')
label=np.array(lab.column('Label').to_pylist())
dec=np.array(pq.read_table(f'{D}/tensor_ih1_meta.parquet').column('Decoy').to_pylist())
X=np.load(f'{D}/tensor_ih1_traces.npy',mmap_mode='r'); M=np.load(f'{D}/tensor_ih1_mask.npy')
REL=np.load(f'{D}/desc_relint.npy'); A0=np.load(f'{D}/apex_ih1.npy')
rng=np.random.default_rng(0)
pos=rng.permutation(np.flatnonzero((label=='pos')&(A0>=0)))[:9000]
dcy=rng.permutation(np.flatnonzero((dec==1)&(A0>=0)))[:9000]
sel=np.concatenate([pos,dcy]); istgt=np.concatenate([np.ones(len(pos),bool),np.zeros(len(dcy),bool)])
R_={'lib_raw':[],'lib_clean':[],'explained':[]}; keep=[]
for a0 in range(0,len(sel),300):
    idx=sel[a0:a0+300]; T=np.asarray(X[idx],dtype=np.float32)
    for j in range(len(idx)):
        i=idx[j]; r=scores(T[j].astype(float),M[i],REL[i].astype(float),int(A0[i]))
        if r is None: continue
        keep.append(istgt[a0+j])
        R_['lib_raw'].append(r[0]); R_['lib_clean'].append(r[1]); R_['explained'].append(r[2])
keep=np.array(keep)
print(f'{len(keep):,} candidates ({keep.sum():,} DIA-NN-confident, {(~keep).sum():,} decoy)\n')
print(f'{"score":>12}' + ''.join(f'{"top-"+str(k):>11}' for k in (500,1000,2000,5000)))
for k_ in ('lib_raw','lib_clean','explained'):
    v=np.array(R_[k_]); o=np.argsort(-v); row=f'{k_:>12}'
    for n in (500,1000,2000,5000): row+=f'{100*keep[o[:n]].mean():>10.1f}%'
    print(row)
print('\n  Target fraction among the top-N, decoys as control.')
