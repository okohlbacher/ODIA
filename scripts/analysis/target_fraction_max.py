#!/usr/bin/env python3
"""The reference-free ruler, replacing the one the vault disowns.

Every selection number in this round was "selection accuracy" -- how often the
chosen candidate agrees with DIA-NN. The vault is explicit that this is a
ceiling and not a result: recall of another tool's ID set is positive-unlabelled
by construction, and it records a case where that ceiling rose from 85.1% to
93.8% while realised identifications FELL by 319 peptides. Recall rising is
never sufficient evidence.

The prescribed replacement uses no external tool at all: TARGET FRACTION AMONG
THE TOP-N by the score being proposed, at the intended operating point. Decoys
are the control, so it measures discrimination rather than agreement, and a
change that helps only on the precursors DIA-NN happens to find cannot hide in
it.

Selection statistic varies between arms. The discriminant that ranks winners is
held FIXED, so what is measured is whether choosing a better candidate improves
target/decoy separation -- not whether a different score sorts differently.

Reported against THREE discriminants, because using one that the selection also
optimises is a winner's curse trap: an arm that picks the position maximising
`lc` and is then ranked by something containing `lc` lets DECOYS shop for their
best `lc` too, inflating the decoy tail symmetrically and eroding separation for
reasons that have nothing to do with whether the selection was better. `coel`
alone is coupled to corr_sum instead, so it is biased the other way. Neither is
clean on its own; the pair brackets the answer.
"""
import numpy as np, pyarrow.parquet as pq
D='/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'
S, SMOOTH, MINCORR, MAXDIFF, APEXEV, HALF, CAP = 4, 2, 0.5, 2.0, 0.99, 2, 3
NPREC = 9000   # per class

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
def win(tr,rel,k,half,n):
    lo=max(k-half,0); hi=min(k+half,n-1)
    seg=tr[:,lo:hi+1]
    fl=np.concatenate([tr[:,max(lo-20,0):lo],tr[:,hi+1:hi+21]],axis=1)
    base=np.median(fl,axis=1) if fl.shape[1] else np.zeros(tr.shape[0])
    area=np.clip(seg-base[:,None],0,None).sum(1)
    f=rel>0
    lc=-1.0
    if f.sum()>=3:
        a,l=area[f],rel[f]
        if a.std()>0 and l.std()>0:
            r=np.corrcoef(a,l)[0,1]
            if np.isfinite(r): lc=float(r)
    s=seg-seg.mean(1,keepdims=True); nn=np.linalg.norm(s,axis=1); ok=nn>0
    co=-1.0
    if ok.sum()>=3:
        u=s[ok]/nn[ok][:,None]; C=u@u.T; iu=np.triu_indices(len(u),1)
        co=float(C[iu].mean())
    return lc,co

meta=pq.read_table(f'{D}/tensor_s08_meta.parquet')
dec=np.array(meta.column('Decoy').to_pylist())
X=np.load(f'{D}/tensor_s08_traces.npy',mmap_mode='r'); M=np.load(f'{D}/tensor_s08_mask.npy')
REL=np.load(f'{D}/desc_relint.npy')
rng=np.random.default_rng(0)
tgt=rng.permutation(np.flatnonzero(dec==0))[:NPREC]
dcy=rng.permutation(np.flatnonzero(dec==1))[:NPREC]
sel=np.concatenate([tgt,dcy]); istgt=np.concatenate([np.ones(len(tgt),bool),np.zeros(len(dcy),bool)])

ARMS=['corr_sum','corr+lib','corr_max']
DISC=['coel+lib','coel','lib']
score={(a,d):[] for a in ARMS for d in DISC}; lab={(a,d):[] for a in ARMS for d in DISC}
for a0 in range(0,len(sel),200):
    idx=sel[a0:a0+200]; lb=istgt[a0:a0+200]
    T=np.asarray(X[idx],dtype=np.float32)
    for j in range(len(idx)):
        i=idx[j]; msk=M[i]; tr=T[j][msk].astype(float); F=tr.shape[0]
        if F<2: continue
        n=tr.shape[1]; rel=REL[i][msk].astype(float)
        sm=np.stack([movavg(tr[f],SMOOTH) for f in range(F)])
        ks=np.arange(S+1,n-S-2)
        if len(ks)==0: continue
        ii=ks[:,None]+np.arange(-S,S+1)[None,:]
        seg=tr[:,ii]; seg=seg-seg.mean(2,keepdims=True)
        nn=np.linalg.norm(seg,axis=2)
        with np.errstate(invalid='ignore',divide='ignore'): u=seg/nn[:,:,None]
        u[~np.isfinite(u)]=0.0
        C=np.einsum('ikw,jkw->kij',u,u)
        for f in range(F): C[:,f,f]=0.0
        sc=C.sum(2)
        pres=np.array([(np.count_nonzero(tr[:,k-1:k+2].sum(1)>0))>=2 for k in ks])
        half=max(S//3,1); e=S-1 if S>1 else 1
        rm_h=rollmax(sm,half)[:,ks]; rm_e=rollmax(sm,e)[:,ks]; smk=sm[:,ks]
        valid=(smk>0)&(smk>=rm_h)&((rm_e<=0)|(smk>=APEXEV*rm_e))
        ok=(sc>=MINCORR)&valid.T
        hit=np.where(ok,sc,-np.inf).max(1)
        hit=np.where(pres&np.isfinite(hit),hit,np.nan)
        # The position's BEST correlation, whichever fragment carried it.
        # Admission is unchanged -- a position still needs a witness passing the
        # shape gates -- only the value it is then RANKED by differs. Today that
        # value is the first passing fragment's sum, which is not the argmax for
        # 67.7% of positions and is deflated by a median of 0.709 against a
        # margin of 2.0.
        hmax=np.where(np.isfinite(hit), sc.max(1), np.nan)
        h=np.flatnonzero(np.isfinite(hit))
        if len(h)==0: continue
        best=np.nanmax(hit); keep=h[hit[h]>=best-MAXDIFF]
        if len(keep)==0: continue
        W=[win(tr,rel,int(ks[p]),HALF,n) for p in keep]
        lc=np.array([w[0] for w in W]); co=np.array([w[1] for w in W])
        hs=hit[keep]; rg=hs.max()-hs.min()
        hsn=(hs-hs.min())/rg if rg>0 else np.zeros_like(hs)
        hm=hmax[keep]
        for arm in ARMS:
            key = hs if arm=='corr_sum' else (hm if arm=='corr_max' else hsn+lc)
            order=np.argsort(-key,kind='stable')[:CAP]
            for dn,dv in (('coel+lib',co[order]+lc[order]),
                          ('coel',co[order]),
                          ('lib',lc[order])):
                score[(arm,dn)].append(float(dv.max())); lab[(arm,dn)].append(bool(lb[j]))

k0=('corr_sum','coel')
print(f'{len(score[k0]):,} precursors ({int(np.sum(lab[k0])):,} target, '
      f'{int(len(lab[k0])-np.sum(lab[k0])):,} decoy)\n')
for dn in DISC:
    print(f'  discriminant = {dn}')
    print(f'{"selection":>14}' + ''.join(f'{"top-"+str(k):>10}' for k in (500,1000,2000,5000)))
    for arm in ARMS:
        sv=np.array(score[(arm,dn)]); l=np.array(lab[(arm,dn)])
        o=np.argsort(-sv); row=f'{arm:>14}'
        for k in (500,1000,2000,5000):
            row+=f'{100*l[o[:k]].mean():>9.1f}%'
        print(row)
    print()
print('\n  Target fraction among the top-N, decoys as the control. No DIA-NN in this\n'
      '  table at all. The final column is the pool composition, i.e. chance.')
