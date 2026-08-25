#!/usr/bin/env python3
"""Phase-4 candidate, measured before it is written: per-fragment boundaries
with a cross-fragment consensus, the way OpenSWATH's recalculatePeakBorders_
does it -- and border dispersion as a sub-score in its own right.

ODIA (and ODIA's OpenSWATH picker) determine boundaries once, on the SUMMED
trace. OpenSWATH determines them per transition and takes a consensus: if the
best transition's border is a z-score outlier against the others, it is
replaced by the median. The fragments of one eluting molecule agree on where
the peak starts and stops; an interference's do not. So the agreement is both
a better estimator and, separately, a discriminant.

Three questions, all on the same paired correct/wrong candidates:
  1. do consensus bounds separate better than summed-trace bounds?
  2. is border dispersion itself discriminative, and how strongly?
  3. does it add anything over co-elution, or is it the same signal renamed?
"""
import sys, re, numpy as np, pyarrow.parquet as pq
D = '/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'; SP = 1.385
def smooth(v, h): return np.convolve(v, np.ones(2*h+1)/(2*h+1), mode='same')

def walk(sm, apex, frac=0.10, sigmas=1.0, minc=7, maxh=20, snap=True):
    n = len(sm)
    if snap:
        lo = max(apex-2,0); hi = min(apex+2,n-1); apex = int(lo+np.argmax(sm[lo:hi+1]))
    a = sm[apex]
    if not np.isfinite(a) or a <= 0: return apex, apex
    f = frac*a; core = 6
    far = np.concatenate([sm[:max(apex-core,0)], sm[min(apex+core+1,n):]])
    if far.size >= 8:
        base = np.median(far); s = 1.4826*np.median(np.abs(far-base))
        f = max(f, min(base+sigmas*s, a-1e-9*max(1.0,abs(a))))
    reb = 0.25*a
    l = apex; run = a; g = 0
    while l > 0 and g < maxh:
        v = sm[l-1]
        if v <= f or v > run+reb: break
        run = min(run,v); l -= 1; g += 1
    r = apex; run = a; g = 0
    while r+1 < n and g < maxh:
        v = sm[r+1]
        if v <= f or v > run+reb: break
        run = min(run,v); r += 1; g += 1
    while r-l+1 < minc and (l > 0 or r+1 < n):
        if l > 0 and (r+1 >= n or sm[l-1] >= sm[r+1]): l -= 1
        elif r+1 < n: r += 1
        else: break
    return l, r

def frag_borders(tr, msk, apex):
    """Per-fragment boundaries, each fragment walked on its OWN trace."""
    L, R = [], []
    for k in range(tr.shape[0]):
        if not msk[k]: continue
        v = tr[k].astype(float)
        if v.max() <= 0: continue
        sm = smooth(v, 1)
        # no snap: every fragment must be judged at the candidate's position,
        # otherwise each one drifts to its own maximum and the agreement
        # measured below is agreement about nothing.
        l, r = walk(sm, apex, snap=False)
        L.append(l); R.append(r)
    return np.array(L, float), np.array(R, float)

def sc(tr, msk, rel, lo, hi):
    lo, hi = int(lo), int(hi)
    seg = tr[:, lo:hi+1]
    fl = np.concatenate([tr[:, max(lo-20,0):lo], tr[:, hi+1:hi+21]], axis=1)
    base = np.median(fl, axis=1) if fl.shape[1] else np.zeros(tr.shape[0])
    area = np.clip(seg-base[:,None], 0, None).sum(1)
    f = msk & (rel > 0)
    if f.sum() < 3: return np.nan, np.nan
    a, l_ = area[f], rel[f]
    lc = np.corrcoef(a, l_)[0,1] if a.std()>0 and l_.std()>0 else np.nan
    s = seg[f]; s = s - s.mean(1, keepdims=True); nn = np.linalg.norm(s, axis=1)
    ok = nn > 0
    if ok.sum() < 3: return lc, np.nan
    s = s[ok]/nn[ok][:,None]; C = s @ s.T; iu = np.triu_indices(len(s),1)
    return lc, float(C[iu].mean())

def auc(p, n):
    p = np.asarray(p,float); n = np.asarray(n,float)
    m = np.isfinite(p)&np.isfinite(n); p,n = p[m],n[m]
    if len(p) < 10: return np.nan
    a = np.concatenate([p,n]); r = a.argsort().argsort().astype(float)+1
    return (r[:len(p)].sum()-len(p)*(len(p)+1)/2)/(len(p)*len(n))

lab = pq.read_table(f'{D}/tensor_s08_labels.parquet')
ids = np.array(lab.column('Precursor.Id').to_pylist()); label = np.array(lab.column('Label').to_pylist())
RT0 = pq.read_table(f'{D}/tensor_s08_meta.parquet').column('RT0').to_numpy()
X = np.load(f'{D}/tensor_s08_traces.npy', mmap_mode='r'); M = np.load(f'{D}/tensor_s08_mask.npy')
REL = np.load(f'{D}/desc_relint.npy')
ALIAS = {'UniMod:4':'Carbamidomethyl'}
nrm = lambda p: re.sub(r'\((.*?)\)', lambda m:'('+ALIAS.get(m.group(1),m.group(1))+')', str(p))
dn = pq.read_table('/scratch/kohlbach/xic/dn_xic.parquet', columns=['Precursor.Id','RT','Q.Value'])
q = dn.column('Q.Value').to_numpy()
TRT = {nrm(p): r*60.0 for p,r,k in zip(dn.column('Precursor.Id').to_pylist(), dn.column('RT').to_numpy(), q<=0.01) if k}

sel = [i for i in np.flatnonzero(label=='pos')[:12000] if ids[i] in TRT]
K = ['sum_lc','sum_co','con_lc','con_co','disp','width_sum','width_con']
P = {k: [] for k in K}; N = {k: [] for k in K}
for a in range(0, len(sel), 400):
    s = np.array(sel[a:a+400]); T = np.asarray(X[s], dtype=np.float32)
    for j in range(len(s)):
        i = s[j]; tr = T[j]; msk = M[i]; rel = REL[i].astype(float)
        c = int(round((TRT[ids[i]]-RT0[i])/SP))
        if not (6 <= c < tr.shape[1]-6): continue
        tot = (tr*msk[:,None]).sum(0).astype(float); sm = smooth(tot,1)
        far = np.arange(len(sm)); far = far[(np.abs(far-c)>=20)&(far>=6)&(far<len(sm)-6)]
        if far.size == 0: continue
        w = int(far[np.argmax(sm[far])])
        for apex, acc in ((c,P),(w,N)):
            lo, hi = walk(sm, apex)
            lc, co = sc(tr,msk,rel,lo,hi)
            acc['sum_lc'].append(lc); acc['sum_co'].append(co); acc['width_sum'].append(hi-lo+1)
            L, R = frag_borders(tr, msk, apex)
            if len(L) < 3:
                for k in ('con_lc','con_co','disp','width_con'): acc[k].append(np.nan)
                continue
            cl, cr = np.median(L), np.median(R)
            # dispersion in cycles, robust, on both borders; negated so that
            # larger is better and the AUC reads the same way as the others
            d = 1.4826*(np.median(np.abs(L-cl)) + np.median(np.abs(R-cr)))
            acc['disp'].append(-d)
            lc, co = sc(tr,msk,rel,cl,cr)
            acc['con_lc'].append(lc); acc['con_co'].append(co); acc['width_con'].append(cr-cl+1)

print(f'{len(P["sum_lc"]):,} paired candidates\n')
print(f'{"score":<26}{"AUC correct-vs-wrong":>22}')
for k, nm in (('sum_lc','lib_corr  @ summed bounds'), ('con_lc','lib_corr  @ consensus'),
              ('sum_co','coelution @ summed bounds'), ('con_co','coelution @ consensus'),
              ('disp','border dispersion (new)')):
    print(f'  {nm:<24}{auc(P[k],N[k]):>22.3f}')
for k, nm in (('width_sum','summed'), ('width_con','consensus')):
    v = np.array(P[k],float); v = v[np.isfinite(v)]
    print(f'  median width {nm:<12}{np.median(v):>10.0f} cycles')
d1 = np.array(P['disp'],float); d0 = np.array(N['disp'],float)
m = np.isfinite(d1)&np.isfinite(d0)
print(f'\n  border dispersion: correct {-np.median(d1[m]):.2f} cycles, '
      f'wrong {-np.median(d0[m]):.2f} cycles')
co1 = np.array(P['sum_co'],float)
mm = m & np.isfinite(co1)
print(f'  corr(dispersion, coelution) on correct candidates: '
      f'{np.corrcoef(d1[mm], co1[mm])[0,1]:+.3f}  '
      f'(near 0 => genuinely new information, near 1 => the same signal renamed)')
