#!/usr/bin/env python3
"""Codex: "A useful ablation is simply fixed +-3 cycles. If fixed seven performs
as well as the full floor/rebound algorithm, the more complicated boundary logic
has not earned its complexity."

Agreed, and it is the right question to ask of any rule with four constants in
it. Same paired separation probe, four boundary rules, stratified by how far the
wrong candidate sits.
"""
import sys, re, numpy as np, pyarrow.parquet as pq
D='/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'; SP=1.385
def smooth(v,h): return np.convolve(v,np.ones(2*h+1)/(2*h+1),mode='same')
def snap(sm,a):
    lo=max(a-2,0); hi=min(a+2,len(sm)-1); return int(lo+np.argmax(sm[lo:hi+1]))
def adaptive(sm,apex,minc=5,frac=0.10,sigmas=1.0,maxh=20):
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
def fixed(sm,apex,half):
    a=snap(sm,apex); return max(a-half,0),min(a+half,len(sm)-1)
def old(sm,apex,frac=0.10):
    n=len(sm); f=frac*sm[apex]; l=apex; r=apex
    while l>0 and sm[l-1]>f: l-=1
    while r+1<n and sm[r+1]>f: r+=1
    return l,r
RULES=[('old (shipped before)',old),('fixed +-2 (w 5)',lambda s,a:fixed(s,a,2)),
       ('fixed +-3 (w 7)',lambda s,a:fixed(s,a,3)),('adaptive min 5',adaptive)]
def sc(tr,msk,rel,lo,hi):
    lo,hi=int(lo),int(hi); seg=tr[:,lo:hi+1]
    fl=np.concatenate([tr[:,max(lo-20,0):lo],tr[:,hi+1:hi+21]],axis=1)
    base=np.median(fl,axis=1) if fl.shape[1] else np.zeros(tr.shape[0])
    area=np.clip(seg-base[:,None],0,None).sum(1); f=msk&(rel>0)
    if f.sum()<3: return np.nan,np.nan
    a,l_=area[f],rel[f]
    lc=np.corrcoef(a,l_)[0,1] if a.std()>0 and l_.std()>0 else np.nan
    s=seg[f]; s=s-s.mean(1,keepdims=True); nn=np.linalg.norm(s,axis=1); ok=nn>0
    if ok.sum()<3: return lc,np.nan
    s=s[ok]/nn[ok][:,None]; C=s@s.T; iu=np.triu_indices(len(s),1)
    return lc,float(C[iu].mean())
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
A={'UniMod:4':'Carbamidomethyl'}
nrm=lambda p: re.sub(r'\((.*?)\)',lambda m:'('+A.get(m.group(1),m.group(1))+')',str(p))
dn=pq.read_table('/scratch/kohlbach/xic/dn_xic.parquet',columns=['Precursor.Id','RT','Q.Value'])
q=dn.column('Q.Value').to_numpy()
TRT={nrm(p):r*60.0 for p,r,k in zip(dn.column('Precursor.Id').to_pylist(),dn.column('RT').to_numpy(),q<=0.01) if k}
BANDS=[('3-6',3,6),('7-19',7,19),('>=20',20,999)]
acc={(rn,bn,s):{'p':[],'n':[]} for rn,_ in RULES for bn,_,_ in BANDS for s in ('lc','co')}
wid={rn:[] for rn,_ in RULES}
sel=[i for i in np.flatnonzero(label=='pos')[:12000] if ids[i] in TRT]
for a0 in range(0,len(sel),400):
    s=np.array(sel[a0:a0+400]); T=np.asarray(X[s],dtype=np.float32)
    for j in range(len(s)):
        i=s[j]; tr=T[j]; msk=M[i]; rel=REL[i].astype(float)
        c=int(round((TRT[ids[i]]-RT0[i])/SP))
        if not (6<=c<tr.shape[1]-6): continue
        tot=(tr*msk[:,None]).sum(0).astype(float); sm=smooth(tot,1)
        idx=np.arange(len(sm)); ok=(idx>=6)&(idx<len(sm)-6)
        for rn,fn in RULES:
            lo,hi=fn(sm,c); wid[rn].append(hi-lo+1)
            lc,co=sc(tr,msk,rel,lo,hi)
            for bn,dlo,dhi in BANDS:
                band=ok&(np.abs(idx-c)>=dlo)&(np.abs(idx-c)<=dhi)
                if not band.any(): continue
                w=int(idx[band][np.argmax(sm[band])])
                lo2,hi2=fn(sm,w); lc2,co2=sc(tr,msk,rel,lo2,hi2)
                acc[(rn,bn,'lc')]['p'].append(lc); acc[(rn,bn,'lc')]['n'].append(lc2)
                acc[(rn,bn,'co')]['p'].append(co); acc[(rn,bn,'co')]['n'].append(co2)
print(f'{len(sel):,} positives. Separation AUC, wrong apex by distance in cycles.\n')
print(f'{"rule":<22}{"med w":>7}' + ''.join(f'{b+" lc":>10}{b+" co":>10}' for b,_,_ in BANDS))
for rn,_ in RULES:
    row=f'{rn:<22}{np.median(wid[rn]):>7.0f}'
    for bn,_,_ in BANDS:
        row+=f'{auc(acc[(rn,bn,"lc")]["p"],acc[(rn,bn,"lc")]["n"]):>10.3f}'
        row+=f'{auc(acc[(rn,bn,"co")]["p"],acc[(rn,bn,"co")]["n"]):>10.3f}'
    print(row)
