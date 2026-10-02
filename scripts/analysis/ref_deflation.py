#!/usr/bin/env python3
"""The margin test compares corr_sum values that may come from different
reference fragments. Does that matter?

findCandidatesByCorrelation walks fragments in descending corr_sum and takes the
FIRST that passes the shape gates, storing that fragment's sum as the position's
score. So a position whose best fragment fails the local-max or apex-evidence
test is scored by a weaker fragment -- and then judged against `max_corr_diff`
on that deflated value. A position can lose the margin not because its evidence
is weak but because its strongest fragment was slightly off-apex.

Two numbers decide whether this is a real effect or a curiosity: how often the
passing reference is not the argmax, and how large the deflation is relative to
the margin it is then tested against.
"""
import numpy as np, pyarrow.parquet as pq
D='/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'
S, SMOOTH, MINCORR, APEXEV, MAXDIFF = 4, 2, 0.5, 0.99, 2.0
def movavg(v,h):
    n=len(v); c=np.cumsum(np.insert(v,0,0.0))
    lo=np.maximum(np.arange(n)-h,0); hi=np.minimum(np.arange(n)+h,n-1)
    return (c[hi+1]-c[lo])/(hi-lo+1)
def rollmax(a,h):
    out=a.copy()
    for d in range(1,h+1):
        out=np.maximum(out,np.pad(a[:,d:],((0,0),(0,d)),constant_values=-np.inf))
        out=np.maximum(out,np.pad(a[:,:-d],((0,0),(d,0)),constant_values=-np.inf))
    return out
meta=pq.read_table(f'{D}/tensor_ih1_meta.parquet')
dec=np.array(meta.column('Decoy').to_pylist())
X=np.load(f'{D}/tensor_ih1_traces.npy',mmap_mode='r'); M=np.load(f'{D}/tensor_ih1_mask.npy')
rng=np.random.default_rng(0)
sel=rng.permutation(np.flatnonzero(dec==0))[:6000]
n_hit=0; n_notmax=0; defl=[]; flip=0; n_prec=0; flipped_prec=0
for a0 in range(0,len(sel),200):
    idx=sel[a0:a0+200]; T=np.asarray(X[idx],dtype=np.float32)
    for j in range(len(idx)):
        i=idx[j]; msk=M[i]; tr=T[j][msk].astype(float); F=tr.shape[0]
        if F<2: continue
        n=tr.shape[1]
        sm=np.stack([movavg(tr[f],SMOOTH) for f in range(F)])
        ks=np.arange(S,n-S)
        if len(ks)==0: continue
        ii=ks[:,None]+np.arange(-S,S+1)[None,:]
        seg=tr[:,ii]; seg=seg-seg.mean(2,keepdims=True)
        nn=np.linalg.norm(seg,axis=2)
        with np.errstate(invalid='ignore',divide='ignore'): u=seg/nn[:,:,None]
        u[~np.isfinite(u)]=0.0
        C=np.einsum('ikw,jkw->kij',u,u)
        for f in range(F): C[:,f,f]=0.0
        sc=C.sum(2)
        half=max(S//3,1); e=S-1 if S>1 else 1
        rm_h=rollmax(sm,half)[:,ks]; rm_e=rollmax(sm,e)[:,ks]; smk=sm[:,ks]
        valid=((smk>0)&(smk>=rm_h)&((rm_e<=0)|(smk>=APEXEV*rm_e))).T   # (K,F)
        okc=(sc>=MINCORR)&valid
        has=okc.any(1)
        if not has.any(): continue
        n_prec+=1
        stored=np.where(okc,sc,-np.inf).max(1)          # first PASSING == max over passing
        argmaxall=sc.max(1)                            # what the position could have scored
        h=np.flatnonzero(has)
        n_hit+=len(h)
        d=argmaxall[h]-stored[h]
        n_notmax+=int((d>1e-9).sum())
        defl.extend(d[d>1e-9].tolist())
        # would the margin winner change if positions were scored by their max?
        if len(h)>1:
            w1=h[np.argmax(stored[h])]; w2=h[np.argmax(argmaxall[h])]
            if w1!=w2: flipped_prec+=1
            keep1=set(h[stored[h]>=stored[h].max()-MAXDIFF].tolist())
            keep2=set(h[argmaxall[h]>=argmaxall[h].max()-MAXDIFF].tolist())
            if keep1!=keep2: flip+=1
d=np.array(defl)
print(f'{n_prec:,} precursors, {n_hit:,} hit positions\n')
print(f'  the scored reference is NOT the highest-correlating fragment: '
      f'{100*n_notmax/n_hit:.1f}% of positions')
if len(d):
    print(f'  deflation when it happens (corr_sum units, margin is {MAXDIFF}):')
    print(f'    median {np.median(d):.3f}   p90 {np.percentile(d,90):.3f}   '
          f'max {d.max():.3f}')
    print(f'    exceeds the margin: {100*(d>MAXDIFF).mean():.2f}% of deflated positions')
print(f'\n  precursors whose top-ranked position would change: {100*flipped_prec/n_prec:.2f}%')
print(f'  precursors whose margin-surviving SET would change: {100*flip/n_prec:.2f}%')
