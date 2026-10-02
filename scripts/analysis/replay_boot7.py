#!/usr/bin/env python3
"""RERUN AT THE SHIPPED DEFAULTS (peak_min_cycles=7, not the 5 this was first
run at). Is "fixed beats adaptive" real, or is 0.012 of AUC inside the noise?

Codex's methodological point, taken: the observations are clustered by
precursor, so a confidence interval has to resample PRECURSORS, not candidates.
Paired bootstrap on the AUC difference, the same candidates under both rules.

Split by walked width as well, because the two rules can only differ where the
walk actually walks. Where the peak is narrow the adaptive rule pads to the
minimum and IS the fixed window -- identical bounds, identical scores -- so an
average over all precursors is diluted by cases that carry no information about
the comparison at all.
"""
import re, numpy as np, pyarrow.parquet as pq
D='/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'; SP=1.385
def smooth(v,h): return np.convolve(v,np.ones(2*h+1)/(2*h+1),mode='same')
def snap(sm,a):
    lo=max(a-2,0); hi=min(a+2,len(sm)-1); return int(lo+np.argmax(sm[lo:hi+1]))
def adaptive(sm,apex,minc=7,frac=0.10,sigmas=1.0,maxh=20):
    n=len(sm); apex=snap(sm,apex); a=sm[apex]
    if not np.isfinite(a) or a<=0: return apex,apex
    f=frac*a; core=6
    far=np.concatenate([sm[:max(apex-core,0)],sm[min(apex+core+1,n):]])
    if far.size>=8:
        base=np.median(far); s=1.4826*np.median(np.abs(far-base))
        f=max(f,min(base+sigmas*s,np.nextafter(a,-np.inf)))
    reb=0.25*a
    def side(st,p):
        run=a; g=0
        while 0<=p+st<n and g<maxh:
            v=sm[p+st]
            if v<=f or v>run+reb: break
            run=min(run,v); p+=st; g+=1
        return p
    l=side(-1,apex); r=side(1,apex)
    while r-l+1<minc and (l>0 or r+1<n):
        if l>0 and (r+1>=n or sm[l-1]>=sm[r+1]): l-=1
        elif r+1<n: r+=1
        else: break
    return l,r
def fixed(sm,apex,half=2):
    a=snap(sm,apex); return max(a-half,0),min(a+half,len(sm)-1)
def sc(tr,msk,rel,lo,hi):
    lo,hi=int(lo),int(hi); seg=tr[:,lo:hi+1]
    fl=np.concatenate([tr[:,max(lo-20,0):lo],tr[:,hi+1:hi+21]],axis=1)
    base=np.median(fl,axis=1) if fl.shape[1] else np.zeros(tr.shape[0])
    area=np.clip(seg-base[:,None],0,None).sum(1); f=msk&(rel>0)
    if f.sum()<3: return np.nan
    s=seg[f]; s=s-s.mean(1,keepdims=True); nn=np.linalg.norm(s,axis=1); ok=nn>0
    if ok.sum()<3: return np.nan
    s=s[ok]/nn[ok][:,None]; C=s@s.T; iu=np.triu_indices(len(s),1)
    return float(C[iu].mean())
def auc(p,n):
    p=np.asarray(p,float); n=np.asarray(n,float); m=np.isfinite(p)&np.isfinite(n)
    p,n=p[m],n[m]
    if len(p)<20: return np.nan
    a=np.concatenate([p,n]); r=a.argsort().argsort().astype(float)+1
    return (r[:len(p)].sum()-len(p)*(len(p)+1)/2)/(len(p)*len(n))
lab=pq.read_table(f'{D}/tensor_ih1_labels.parquet')
ids=np.array(lab.column('Precursor.Id').to_pylist()); label=np.array(lab.column('Label').to_pylist())
RT0=pq.read_table(f'{D}/tensor_ih1_meta.parquet').column('RT0').to_numpy()
X=np.load(f'{D}/tensor_ih1_traces.npy',mmap_mode='r'); M=np.load(f'{D}/tensor_ih1_mask.npy')
REL=np.load(f'{D}/desc_relint.npy')
AL={'UniMod:4':'Carbamidomethyl'}
nrm=lambda p: re.sub(r'\((.*?)\)',lambda m:'('+AL.get(m.group(1),m.group(1))+')',str(p))
dn=pq.read_table('/scratch/kohlbach/xic/dn_xic.parquet',columns=['Precursor.Id','RT','Q.Value'])
q=dn.column('Q.Value').to_numpy()
TRT={nrm(p):r*60.0 for p,r,k in zip(dn.column('Precursor.Id').to_pylist(),dn.column('RT').to_numpy(),q<=0.01) if k}
sel=[i for i in np.flatnonzero(label=='pos')[:12000] if ids[i] in TRT]
rows=[]   # (adaptive_pos, adaptive_neg, fixed_pos, fixed_neg, walked_width, identical)
for a0 in range(0,len(sel),400):
    s=np.array(sel[a0:a0+400]); T=np.asarray(X[s],dtype=np.float32)
    for j in range(len(s)):
        i=s[j]; tr=T[j]; msk=M[i]; rel=REL[i].astype(float)
        c=int(round((TRT[ids[i]]-RT0[i])/SP))
        if not (6<=c<tr.shape[1]-6): continue
        tot=(tr*msk[:,None]).sum(0).astype(float); sm=smooth(tot,1)
        idx=np.arange(len(sm)); ok=(idx>=6)&(idx<len(sm)-6)&(np.abs(idx-c)>=7)
        if not ok.any(): continue
        w=int(idx[ok][np.argmax(sm[ok])])
        la,ha=adaptive(sm,c); lf,hf=fixed(sm,c)
        la2,ha2=adaptive(sm,w); lf2,hf2=fixed(sm,w)
        rows.append((sc(tr,msk,rel,la,ha), sc(tr,msk,rel,la2,ha2),
                     sc(tr,msk,rel,lf,hf), sc(tr,msk,rel,lf2,hf2),
                     ha-la+1, (la,ha)==(lf,hf) and (la2,ha2)==(lf2,hf2)))
R=np.array([r[:5] for r in rows],float); same=np.array([r[5] for r in rows],bool)
print(f'{len(R):,} precursors, wrong apex >=7 cycles away\n')
print(f'  the two rules give IDENTICAL bounds for both candidates: '
      f'{100*same.mean():.1f}% of precursors')
def diff(mask):
    a=auc(R[mask,0],R[mask,1]); f=auc(R[mask,2],R[mask,3]); return a,f,f-a
for nm, mask in (('all precursors', np.ones(len(R),bool)),
                 ('where the rules differ', ~same),
                 ('walked width >= 9', R[:,4]>=9)):
    a,f,d = diff(mask)
    print(f'  {nm:<26} n={mask.sum():>6}  adaptive {a:.3f}  fixed {f:.3f}  diff {d:+.3f}')
rng=np.random.default_rng(0); n=len(R); ds=[]
for _ in range(2000):
    b=rng.integers(0,n,n)
    a=auc(R[b,0],R[b,1]); f=auc(R[b,2],R[b,3])
    if np.isfinite(a) and np.isfinite(f): ds.append(f-a)
ds=np.array(ds)
print(f'\n  paired bootstrap over precursors (2,000 resamples), fixed - adaptive:')
print(f'    mean {ds.mean():+.4f}   95% CI [{np.percentile(ds,2.5):+.4f}, '
      f'{np.percentile(ds,97.5):+.4f}]   P(fixed better) {100*(ds>0).mean():.1f}%')
