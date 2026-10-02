#!/usr/bin/env python3
"""41.7% of DIA-NN-confident precursors get no candidate near the true RT.
Is the true position rejected by a gate, or found and then dropped by the cap?

The detector's reject counters cannot answer this: they count over all 128
positions, and these precursors DO receive candidates, just elsewhere.

Every position is evaluated, not only the true one -- otherwise the true
position's RANK cannot be computed, and rank is the whole question. Gates, in
the C++ order (findCandidatesByCorrelation):

  presence   >=2 fragments non-zero somewhere in {k-1, k, k+1}
  then, walking fragments in DESCENDING corr_sum and taking the first that
  passes everything -- which is the same as the highest-scoring fragment that
  passes, since the walk is ordered:
     corr      summed pairwise correlation over [k-S, k+S] >= min_corr_score
     ref>0     that fragment's smoothed value at k is positive
     local max and is the maximum over +-max(S/3,1)
     apex ev.  and is >= apex_evidence x the maximum over +-(S-1)
  then across surviving positions:
     margin    corr_sum within max_corr_diff of the precursor's best
     cap       and inside the top max_candidates
"""
import re, numpy as np, pyarrow.parquet as pq
D='/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'; SP=1.385
S, SMOOTH, MINCORR, MAXDIFF, APEXEV, CAP = 4, 2, 0.5, 2.0, 0.99, 3

def movavg(v, h):
    n=len(v); c=np.cumsum(np.insert(v,0,0.0))
    lo=np.maximum(np.arange(n)-h,0); hi=np.minimum(np.arange(n)+h,n-1)
    return (c[hi+1]-c[lo])/(hi-lo+1)

def rollmax(a, h):
    """Max over +-h, per row, edges clamped."""
    F,n = a.shape
    out = a.copy()
    for d in range(1, h+1):
        out = np.maximum(out, np.pad(a[:,d:], ((0,0),(0,d)), constant_values=-np.inf))
        out = np.maximum(out, np.pad(a[:,:-d], ((0,0),(d,0)), constant_values=-np.inf))
    return out

lab=pq.read_table(f'{D}/tensor_ih1_labels.parquet')
ids=np.array(lab.column('Precursor.Id').to_pylist()); label=np.array(lab.column('Label').to_pylist())
RT0=pq.read_table(f'{D}/tensor_ih1_meta.parquet').column('RT0').to_numpy()
X=np.load(f'{D}/tensor_ih1_traces.npy',mmap_mode='r'); M=np.load(f'{D}/tensor_ih1_mask.npy')
AL={'UniMod:4':'Carbamidomethyl'}
nrm=lambda p: re.sub(r'\((.*?)\)',lambda m:'('+AL.get(m.group(1),m.group(1))+')',str(p))
dn=pq.read_table('/scratch/kohlbach/xic/dn_xic.parquet',columns=['Precursor.Id','RT','Q.Value'])
q=dn.column('Q.Value').to_numpy()
TRT={nrm(p):r*60.0 for p,r,k in zip(dn.column('Precursor.Id').to_pylist(),dn.column('RT').to_numpy(),q<=0.01) if k}

sel=[i for i in np.flatnonzero(label=='pos')[:12000] if ids[i] in TRT]
tot=0; lost={k:0 for k in ('presence','corr','shape','margin','cap','kept')}
ranks=[]; nhits=[]
for a0 in range(0,len(sel),200):
    s=np.array(sel[a0:a0+200]); T=np.asarray(X[s],dtype=np.float32)
    for j in range(len(s)):
        i=s[j]; msk=M[i]; tr=T[j][msk].astype(float)
        F=tr.shape[0]
        if F<2: continue
        n=tr.shape[1]
        c=int(round((TRT[ids[i]]-RT0[i])/SP))
        if not (S+1<=c<n-S-2): continue
        tot+=1
        sm=np.stack([movavg(tr[f],SMOOTH) for f in range(F)])
        ks=np.arange(S+1,n-S-2)
        idx=ks[:,None]+np.arange(-S,S+1)[None,:]
        seg=tr[:,idx]; seg=seg-seg.mean(2,keepdims=True)
        nn=np.linalg.norm(seg,axis=2)
        with np.errstate(invalid='ignore',divide='ignore'): u=seg/nn[:,:,None]
        u[~np.isfinite(u)]=0.0
        C=np.einsum('ikw,jkw->kij',u,u)
        for f in range(F): C[:,f,f]=0.0
        sc=C.sum(2)                                   # (K, F)
        pres=np.array([(np.count_nonzero(tr[:,k-1:k+2].sum(1)>0))>=2 for k in ks])
        half=max(S//3,1); e=S-1 if S>1 else 1
        rm_h=rollmax(sm,half)[:,ks]; rm_e=rollmax(sm,e)[:,ks]; smk=sm[:,ks]
        valid=(smk>0)&(smk>=rm_h)&((rm_e<=0)|(smk>=APEXEV*rm_e))   # (F, K)
        ok=(sc>=MINCORR)&valid.T
        hit=np.where(ok, sc, -np.inf).max(1)
        hit=np.where(pres&np.isfinite(hit), hit, np.nan)
        kpos=int(np.searchsorted(ks,c))
        if kpos>=len(ks) or ks[kpos]!=c: continue
        h=np.flatnonzero(np.isfinite(hit))
        if len(h)==0: continue
        nhits.append(len(h))
        if not np.isfinite(hit[kpos]):
            if not pres[kpos]: lost['presence']+=1
            elif sc[kpos].max()<MINCORR: lost['corr']+=1
            else: lost['shape']+=1
            continue
        best=np.nanmax(hit)
        if hit[kpos]<best-MAXDIFF: lost['margin']+=1; continue
        keep=h[hit[h]>=best-MAXDIFF]
        order=keep[np.argsort(-hit[keep],kind='stable')]
        rank=int(np.flatnonzero(order==kpos)[0])+1
        ranks.append(rank)
        if rank>CAP: lost['cap']+=1
        else: lost['kept']+=1

print(f'{tot:,} confident positives with the true apex inside the window')
print(f'  positions that become hits, per precursor: median {np.median(nhits):.0f}, '
      f'p90 {np.percentile(nhits,90):.0f}\n')
print(f'{"where the TRUE position is lost":<38}{"n":>8}{"%":>9}')
for k,nm in (('presence','<2 fragments present'),('corr','corr_sum below 0.5'),
             ('shape','shape gates (local max / apex evid.)'),
             ('margin','outside the margin of 2.0'),('cap','outside the top-3 cap'),
             ('kept','SURVIVES -- emitted as a candidate')):
    print(f'  {nm:<36}{lost[k]:>8,}{100*lost[k]/tot:>8.1f}%')
r=np.array(ranks)
print(f'\n  where it survives the gates and the margin, its rank by corr_sum:')
for k in (1,2,3,5,10):
    print(f'    rank <= {k:<3} {100*(r<=k).sum()/tot:5.1f}% of all positives')
print(f'\n  raising the cap from 3 to 10 would add '
      f'{100*((r<=10).sum()-(r<=3).sum())/tot:.1f} points of recall.')
