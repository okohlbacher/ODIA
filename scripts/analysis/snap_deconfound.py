#!/usr/bin/env python3
"""How much of doc/64's headline is the apex snap rather than the boundaries?

A reviewer flagged the confound. It does NOT apply to the +0.0136 fixed-vs-
walked bootstrap -- both arms there snap, so that comparison was already
controlled. It DOES apply to doc/64's headline 0.545 -> 0.78, because the old
rule had no snap and the new one does, so three things changed at once: the
floor became noise-relative, the width collapsed, and the centre moved.

Separating them, 2x2 over the same paired candidates:
   width   in {walked bounds, fixed +-2}
   centre  in {seed k, snapped to the local max within +-2}
"""
import re, numpy as np, pyarrow.parquet as pq
D='/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'; SP=1.385
def smooth(v,h): return np.convolve(v,np.ones(2*h+1)/(2*h+1),mode='same')
def snap(sm,a):
    lo=max(a-2,0); hi=min(a+2,len(sm)-1); return int(lo+np.argmax(sm[lo:hi+1]))
def walked(sm,apex,minc=7,frac=0.10,sigmas=1.0,maxh=20):
    n=len(sm); a=sm[apex]
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
def old_rule(sm,apex,frac=0.10):
    n=len(sm); f=frac*sm[apex]; l=apex; r=apex
    while l>0 and sm[l-1]>f: l-=1
    while r+1<n and sm[r+1]>f: r+=1
    return l,r
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
    if len(p)<30: return np.nan
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
ARMS=[('old rule (no snap)','old',False),('old rule + snap','old',True),
      ('walked bounds, no snap','walk',False),('walked bounds + snap','walk',True),
      ('fixed +-2, no snap','fix',False),('fixed +-2 + snap','fix',True)]
P={a[0]:[] for a in ARMS}; N={a[0]:[] for a in ARMS}; W={a[0]:[] for a in ARMS}
sel=[i for i in np.flatnonzero(label=='pos')[:12000] if ids[i] in TRT]
for a0 in range(0,len(sel),400):
    s=np.array(sel[a0:a0+400]); T=np.asarray(X[s],dtype=np.float32)
    for j in range(len(s)):
        i=s[j]; tr=T[j]; msk=M[i]; rel=REL[i].astype(float)
        c=int(round((TRT[ids[i]]-RT0[i])/SP))
        if not (6<=c<tr.shape[1]-6): continue
        tot=(tr*msk[:,None]).sum(0).astype(float); sm1=smooth(tot,1); sm2=smooth(tot,2)
        idx=np.arange(len(sm1)); ok=(idx>=6)&(idx<len(sm1)-6)&(np.abs(idx-c)>=20)
        if not ok.any(): continue
        w=int(idx[ok][np.argmax(sm1[ok])])
        for nm,kind,dosnap in ARMS:
            base_sm = sm2 if kind=='old' else sm1
            for apex,acc in ((c,P),(w,N)):
                a2 = snap(base_sm,apex) if dosnap else apex
                if kind=='old': lo,hi=old_rule(base_sm,a2)
                elif kind=='walk': lo,hi=walked(base_sm,a2)
                else: lo,hi=max(a2-2,0),min(a2+2,len(base_sm)-1)
                acc[nm].append(sc(tr,msk,rel,lo,hi))
                if acc is P: W[nm].append(hi-lo+1)
print(f'{len(P[ARMS[0][0]]):,} paired candidates\n')
print(f'{"arm":<26}{"median width":>14}{"AUC":>9}')
for nm,_,_ in ARMS:
    print(f'  {nm:<24}{np.median(W[nm]):>14.0f}{auc(P[nm],N[nm]):>9.3f}')
a=lambda nm: auc(P[nm],N[nm])
print(f'\n  snap alone, holding the old rule:      {a("old rule + snap")-a("old rule (no snap)"):+.3f}')
print(f'  snap alone, holding walked bounds:     {a("walked bounds + snap")-a("walked bounds, no snap"):+.3f}')
print(f'  snap alone, holding the fixed window:  {a("fixed +-2 + snap")-a("fixed +-2, no snap"):+.3f}')
print(f'  boundaries alone (snapped, old->walk): {a("walked bounds + snap")-a("old rule + snap"):+.3f}')
print(f'  width alone (snapped, walk->fixed):    {a("fixed +-2 + snap")-a("walked bounds + snap"):+.3f}')
