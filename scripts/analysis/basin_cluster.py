#!/usr/bin/env python3
"""BASIN CLUSTERING. The gates are fine and the thresholds are fine. The SELECTION STATISTIC is
not: corr_sum puts the true position first only 19.8% of the time, and every
attempt to compensate by admitting more candidates makes selection accuracy
worse. So rank the same hit positions by something better and see if the top 3
then contains the right one.

  corr_sum   summed pairwise Pearson of the RAW traces over +-4 cycles, which
             is what the detector uses today
  coel       mean pairwise correlation over the new fixed 5-cycle scoring
             window -- a real sub-score, computed at every hit position
  both       corr_sum normalised plus coel, an equal-weight blend

Raising recall', AND gives the scorer more chances to
pick a wrong candidate. One number contains both: how often does the SUB-SCORE
select the true position out of the top-K by corr_sum?

Selection accuracy over ALL confident positives, so a K that adds recall but
loses precision shows up as a fall, not a rise. This is the measurement the
-20.4% admission result needs re-run under, because that result predates the
boundary fix and was obtained with sub-scores that could not distinguish
candidates at all (median paired difference exactly 0).

The score used is the co-elution correlation over the new fixed scoring window,
which is one of ODIA's real sub-scores -- not the classifier, which cannot be
run here. So this is a lower bound on what the trained model would achieve, and
a mechanism probe, not an FDP result.
"""
import re, numpy as np, pyarrow.parquet as pq
D='/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'; SP=1.385
S, SMOOTH, MINCORR, MAXDIFF, APEXEV, HALF = 4, 2, 0.5, 2.0, 0.99, 2

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
def libcorr(tr, rel, k, half, n):
    """Do the fragments have the RELATIVE INTENSITIES the library predicts?
    The detector's corr_sum sees only shape; two co-eluting species have equally
    good shape and different intensity patterns, so this is evidence the
    detector does not currently have at selection time."""
    lo=max(k-half,0); hi=min(k+half,n-1)
    seg=tr[:,lo:hi+1]
    fl=np.concatenate([tr[:,max(lo-20,0):lo],tr[:,hi+1:hi+21]],axis=1)
    base=np.median(fl,axis=1) if fl.shape[1] else np.zeros(tr.shape[0])
    area=np.clip(seg-base[:,None],0,None).sum(1)
    f=rel>0
    if f.sum()<3: return -1.0
    a,l=area[f],rel[f]
    if a.std()<=0 or l.std()<=0: return -1.0
    r=float(np.corrcoef(a,l)[0,1])
    return r if np.isfinite(r) else -1.0

def coel(tr, k, half, n):
    lo=max(k-half,0); hi=min(k+half,n-1)
    seg=tr[:,lo:hi+1]
    seg=seg-seg.mean(1,keepdims=True)
    nn=np.linalg.norm(seg,axis=1); ok=nn>0
    if ok.sum()<3: return -np.inf
    u=seg[ok]/nn[ok][:,None]; C=u@u.T; iu=np.triu_indices(len(u),1)
    return float(C[iu].mean())

lab=pq.read_table(f'{D}/tensor_s08_labels.parquet')
ids=np.array(lab.column('Precursor.Id').to_pylist()); label=np.array(lab.column('Label').to_pylist())
RT0=pq.read_table(f'{D}/tensor_s08_meta.parquet').column('RT0').to_numpy()
X=np.load(f'{D}/tensor_s08_traces.npy',mmap_mode='r'); M=np.load(f'{D}/tensor_s08_mask.npy')
REL=np.load(f'{D}/desc_relint.npy')
AL={'UniMod:4':'Carbamidomethyl'}
nrm=lambda p: re.sub(r'\((.*?)\)',lambda m:'('+AL.get(m.group(1),m.group(1))+')',str(p))
dn=pq.read_table('/scratch/kohlbach/xic/dn_xic.parquet',columns=['Precursor.Id','RT','Q.Value'])
q=dn.column('Q.Value').to_numpy()
TRT={nrm(p):r*60.0 for p,r,k in zip(dn.column('Precursor.Id').to_pylist(),dn.column('RT').to_numpy(),q<=0.01) if k}

EVS=[('corr_sum',0),('corr+lib',0),('corr+lib',3),('corr+lib',5),('corr+lib',7),('corr+lib',10),('corr_sum',5)]
sel=[i for i in np.flatnonzero(label=='pos')[:12000] if ids[i] in TRT]
tot=0; hitK={e:0 for e in EVS}; pickK={e:0 for e in EVS}; ncand={e:[] for e in EVS}
for a0 in range(0,len(sel),200):
    s=np.array(sel[a0:a0+200]); T=np.asarray(X[s],dtype=np.float32)
    for j in range(len(s)):
        i=s[j]; msk=M[i]; tr=T[j][msk].astype(float); F=tr.shape[0]
        rel=REL[i][msk].astype(float)
        if F<2: continue
        n=tr.shape[1]; c=int(round((TRT[ids[i]]-RT0[i])/SP))
        if not (S+1<=c<n-S-2): continue
        tot+=1
        sm=np.stack([movavg(tr[f],SMOOTH) for f in range(F)])
        ks=np.arange(S+1,n-S-2); idx=ks[:,None]+np.arange(-S,S+1)[None,:]
        seg=tr[:,idx]; seg=seg-seg.mean(2,keepdims=True)
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
        # library-intensity-weighted variant of the same statistic
        wv=rel/(rel.sum() if rel.sum()>0 else 1.0)
        scw_full=np.einsum('kij,i,j->k',C,wv,wv)*float(F)
        ok=(sc>=MINCORR)&valid.T
        hit=np.where(ok,sc,-np.inf).max(1)
        hit=np.where(pres&np.isfinite(hit),hit,np.nan)
        h=np.flatnonzero(np.isfinite(hit))
        if len(h)==0: continue
        best=np.nanmax(hit); keep=h[hit[h]>=best-MAXDIFF]
        if len(keep)==0: continue
        scw_keep=scw_full[keep]
        cs_all=np.array([coel(tr,int(ks[p]),HALF,n) for p in keep])
        cs_all=np.where(np.isfinite(cs_all),cs_all,-1.0)
        lc_all=np.array([libcorr(tr,rel,int(ks[p]),HALF,n) for p in keep])
        hs=hit[keep]
        rng=hs.max()-hs.min()
        hs_n=(hs-hs.min())/rng if rng>0 else np.zeros_like(hs)
        # a library-intensity-weighted corr_sum: fragments the library says are
        # bright get more of a vote on where the peak is
        w=rel/ (rel.sum() if rel.sum()>0 else 1.0)
        wc=(scw[keep]*1.0) if False else None
        for EV in EVS:
            stat, sep = EV
            key = hs if stat=='corr_sum' else hs_n+lc_all
            ranked=keep[np.argsort(-key,kind='stable')]
            if sep<=0:
                order=ranked[:3]
            else:
                # non-maximum suppression: a cap of 3 should mean three distinct
                # PEAKS, not three samples of possibly one. Without it a single
                # broad basin can occupy every slot and evict the true peak --
                # which is what the non-monotone recall under a relaxed
                # apex_evidence was already saying.
                order=[]
                for p in ranked:
                    if all(abs(int(ks[p])-int(ks[o]))>=sep for o in order):
                        order.append(p)
                    if len(order)>=3: break
                order=np.array(order,dtype=int)
            ncand[EV].append(len(order))
            if len(order)==0: continue
            sub=np.array([coel(tr,int(ks[p]),HALF,n) for p in order])
            if np.any(np.abs(ks[order]-c)<=3): hitK[EV]+=1
            if len(sub) and np.isfinite(sub).any():
                pick=order[int(np.nanargmax(np.where(np.isfinite(sub),sub,-np.inf)))]
                if abs(ks[pick]-c)<=3: pickK[EV]+=1

print(f'{tot:,} confident positives\n')
print(f'{"statistic":>12}{"min sep":>9}{"cands":>8}{"recall@3":>11}{"selection accuracy":>21}')
for EV in EVS:
    st,sp=EV
    print(f'{st:>12}{sp:>9}{np.mean(ncand[EV]):>8.2f}{100*hitK[EV]/tot:>10.1f}%{100*pickK[EV]/tot:>20.1f}%')
print(f'\n  Selection accuracy is the one that matters: it already contains the '
      f'cost of\n  giving the scorer more wrong candidates to choose between.')
