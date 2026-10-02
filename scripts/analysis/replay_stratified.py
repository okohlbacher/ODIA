#!/usr/bin/env python3
"""Two reviewer challenges, one measurement.

Kimi, adversarially:
  (a) the >=20-cycle wrong candidate is the EASY mode -- most of the real
      confusion mass sits at |dapex| <= 3, where a width-7 group shares 5 of its
      7 points with the correct one, and the probe says nothing about it;
  (b) `peak_min_cycles = 7` may be too NARROW, not too wide: DIA-NN computes its
      discriminating correlations over a fixed W = 2S+1 with S = 2.2*PeakWidth,
      which at our FWHM of ~3 cycles is about 12.

So: separation AUC stratified by how far the wrong candidate sits, swept over
min-width, and with DIA-NN's valley-relative rebound guard against our absolute
one. Judged on SEPARATION, which is what the boundary is for -- not on apex
containment, which the previous sweep showed is flat in every knob.
"""
import sys, re, numpy as np, pyarrow.parquet as pq
D = '/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'; SP = 1.385
def smooth(v, h): return np.convolve(v, np.ones(2*h+1)/(2*h+1), mode='same')

def walk(sm, apex, minc, frac=0.10, sigmas=1.0, maxh=20, valley_guard=False):
    n = len(sm)
    lo = max(apex-2,0); hi = min(apex+2,n-1); apex = int(lo+np.argmax(sm[lo:hi+1]))
    a = sm[apex]
    if not np.isfinite(a) or a <= 0: return apex, apex
    f = frac*a; core = 6
    far = np.concatenate([sm[:max(apex-core,0)], sm[min(apex+core+1,n):]])
    if far.size >= 8:
        base = np.median(far); s = 1.4826*np.median(np.abs(far-base))
        f = max(f, min(base+sigmas*s, a-1e-9*max(1.0,abs(a))))
    reb = 0.25*a
    def side(step, start, limit):
        p = start; run = a; g = 0
        while 0 <= p+step < n and g < maxh:
            v = sm[p+step]
            if v <= f: break
            if valley_guard:
                # DIA-NN: a valley only counts if it is DEEP (below apex/3) and
                # the signal has since doubled off it. Valley-relative, so it
                # fires on a shallow saddle between co-eluting isomers that an
                # absolute 0.25*apex threshold walks straight through.
                if run < a/3.0 and v >= 2.0*run: break
            else:
                if v > run + reb: break
            run = min(run, v); p += step; g += 1
        return p
    l = side(-1, apex, None); r = side(1, apex, None)
    while r-l+1 < minc and (l > 0 or r+1 < n):
        if l > 0 and (r+1 >= n or sm[l-1] >= sm[r+1]): l -= 1
        elif r+1 < n: r += 1
        else: break
    return l, r

def sc(tr, msk, rel, lo, hi):
    lo, hi = int(lo), int(hi); seg = tr[:, lo:hi+1]
    fl = np.concatenate([tr[:, max(lo-20,0):lo], tr[:, hi+1:hi+21]], axis=1)
    base = np.median(fl, axis=1) if fl.shape[1] else np.zeros(tr.shape[0])
    area = np.clip(seg-base[:,None], 0, None).sum(1)
    f = msk & (rel > 0)
    if f.sum() < 3: return np.nan, np.nan
    a, l_ = area[f], rel[f]
    lc = np.corrcoef(a,l_)[0,1] if a.std()>0 and l_.std()>0 else np.nan
    s = seg[f]; s = s - s.mean(1,keepdims=True); nn = np.linalg.norm(s,axis=1); ok = nn>0
    if ok.sum() < 3: return lc, np.nan
    s = s[ok]/nn[ok][:,None]; C = s@s.T; iu = np.triu_indices(len(s),1)
    return lc, float(C[iu].mean())

def auc(p, n):
    p=np.asarray(p,float); n=np.asarray(n,float); m=np.isfinite(p)&np.isfinite(n)
    p,n=p[m],n[m]
    if len(p)<30: return np.nan
    a=np.concatenate([p,n]); r=a.argsort().argsort().astype(float)+1
    return (r[:len(p)].sum()-len(p)*(len(p)+1)/2)/(len(p)*len(n))

lab = pq.read_table(f'{D}/tensor_ih1_labels.parquet')
ids = np.array(lab.column('Precursor.Id').to_pylist()); label = np.array(lab.column('Label').to_pylist())
RT0 = pq.read_table(f'{D}/tensor_ih1_meta.parquet').column('RT0').to_numpy()
X = np.load(f'{D}/tensor_ih1_traces.npy', mmap_mode='r'); M = np.load(f'{D}/tensor_ih1_mask.npy')
REL = np.load(f'{D}/desc_relint.npy')
ALIAS={'UniMod:4':'Carbamidomethyl'}
nrm=lambda p: re.sub(r'\((.*?)\)', lambda m:'('+ALIAS.get(m.group(1),m.group(1))+')', str(p))
dn = pq.read_table('/scratch/kohlbach/xic/dn_xic.parquet', columns=['Precursor.Id','RT','Q.Value'])
q = dn.column('Q.Value').to_numpy()
TRT = {nrm(p): r*60.0 for p,r,k in zip(dn.column('Precursor.Id').to_pylist(), dn.column('RT').to_numpy(), q<=0.01) if k}

BANDS = [('close  3-6', 3, 6), ('mid   7-19', 7, 19), ('far    >=20', 20, 999)]
MINS = [5, 7, 9, 12]
acc = {(m,g,b,s): {'p':[], 'n':[]} for m in MINS for g in (0,1) for b,_,_ in BANDS for s in ('lc','co')}
wid = {(m,g): [] for m in MINS for g in (0,1)}

sel = [i for i in np.flatnonzero(label=='pos')[:12000] if ids[i] in TRT]
for a0 in range(0, len(sel), 400):
    s = np.array(sel[a0:a0+400]); T = np.asarray(X[s], dtype=np.float32)
    for j in range(len(s)):
        i = s[j]; tr = T[j]; msk = M[i]; rel = REL[i].astype(float)
        c = int(round((TRT[ids[i]]-RT0[i])/SP))
        if not (6 <= c < tr.shape[1]-6): continue
        tot = (tr*msk[:,None]).sum(0).astype(float); sm = smooth(tot,1)
        idx = np.arange(len(sm)); ok = (idx>=6)&(idx<len(sm)-6)
        for m in MINS:
            for g in (0,1):
                lo, hi = walk(sm, c, m, valley_guard=bool(g))
                wid[(m,g)].append(hi-lo+1)
                lc, co = sc(tr,msk,rel,lo,hi)
                for bn, dlo, dhi in BANDS:
                    band = ok & (np.abs(idx-c)>=dlo) & (np.abs(idx-c)<=dhi)
                    if not band.any(): continue
                    w = int(idx[band][np.argmax(sm[band])])
                    lo2, hi2 = walk(sm, w, m, valley_guard=bool(g))
                    lc2, co2 = sc(tr,msk,rel,lo2,hi2)
                    acc[(m,g,bn,'lc')]['p'].append(lc); acc[(m,g,bn,'lc')]['n'].append(lc2)
                    acc[(m,g,bn,'co')]['p'].append(co); acc[(m,g,bn,'co')]['n'].append(co2)

print(f'{len(sel):,} positives; separation AUC (correct vs wrong) by how far the wrong apex sits\n')
for g, gn in ((0,'absolute rebound (current, 0.25*apex)'), (1,'DIA-NN valley-relative guard')):
    print(f'{gn}')
    print(f'{"min":>5}{"med w":>7}' + ''.join(f'{b+" lc":>15}{b+" co":>15}' for b,_,_ in BANDS))
    for m in MINS:
        row = f'{m:>5}{np.median(wid[(m,g)]):>7.0f}'
        for bn,_,_ in BANDS:
            row += f'{auc(acc[(m,g,bn,"lc")]["p"], acc[(m,g,bn,"lc")]["n"]):>15.3f}'
            row += f'{auc(acc[(m,g,bn,"co")]["p"], acc[(m,g,bn,"co")]["n"]):>15.3f}'
        print(row)
    print()
