#!/usr/bin/env python3
"""Do the narrowed bounds SEPARATE a correct candidate from a wrong one?

Picking is a different problem, measured separately. This asks the only
question phases 1-3 can be held responsible for: given a candidate at the true
apex and a candidate at a wrong position in the same window, do the sub-scores
computed inside the group tell them apart better under the new rule?

Two scores, both computed strictly inside the group, both real sub-scores:
  lib_corr   Pearson of per-fragment integrated area against library intensity
  coel       mean pairwise Pearson of the fragment traces

Reported as AUC over paired (correct, wrong) candidates. An AUC of 0.5 means
the score cannot tell the true position from a false one.

The old rule is EXPECTED to sit near 0.5 -- its group is the whole window, so
both candidates integrate the same interval and get the same number almost by
construction. That is not a subtle finding; it is the defect stated in the form
that shows what it costs.
"""
import sys, re, numpy as np, pyarrow.parquet as pq
D = '/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'; SP = 1.385

def smooth(v, h): return np.convolve(v, np.ones(2*h+1)/(2*h+1), mode='same')

def bounds_old(sm, apex, frac=0.10):
    n = len(sm); f = frac*sm[apex]; l = apex; r = apex
    while l > 0 and sm[l-1] > f: l -= 1
    while r+1 < n and sm[r+1] > f: r += 1
    return l, r

def bounds_new(sm, apex, frac=0.10, sigmas=1.0, minc=7, maxh=20):
    n = len(sm); lo = max(apex-2,0); hi = min(apex+2,n-1)
    apex = int(lo+np.argmax(sm[lo:hi+1])); a = sm[apex]; f = frac*a; core = 6
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

def scores(tr, msk, rel, lo, hi):
    """tr (F,T) fragment traces; integrate inside [lo,hi] over a flank baseline."""
    seg = tr[:, lo:hi+1]
    fl = np.concatenate([tr[:, max(lo-20,0):lo], tr[:, hi+1:hi+21]], axis=1)
    base = np.median(fl, axis=1) if fl.shape[1] else np.zeros(tr.shape[0])
    area = np.clip(seg - base[:, None], 0, None).sum(1)
    f = msk & (rel > 0)
    if f.sum() < 3: return np.nan, np.nan
    a, l_ = area[f], rel[f]
    lc = np.corrcoef(a, l_)[0,1] if a.std() > 0 and l_.std() > 0 else np.nan
    s = seg[f]
    s = s - s.mean(1, keepdims=True)
    nn = np.linalg.norm(s, axis=1)
    ok = nn > 0
    if ok.sum() < 3: return lc, np.nan
    s = s[ok]/nn[ok][:,None]; C = s @ s.T
    iu = np.triu_indices(len(s), 1)
    return lc, float(C[iu].mean())

def auc(pos, neg):
    pos = np.asarray(pos); neg = np.asarray(neg)
    m = np.isfinite(pos) & np.isfinite(neg); pos, neg = pos[m], neg[m]
    if len(pos) < 10: return np.nan, 0
    allv = np.concatenate([pos, neg]); r = allv.argsort().argsort().astype(float)+1
    n1 = len(pos)
    return (r[:n1].sum() - n1*(n1+1)/2) / (n1*len(neg)), n1

lab = pq.read_table(f'{D}/tensor_s08_labels.parquet')
ids = np.array(lab.column('Precursor.Id').to_pylist()); label = np.array(lab.column('Label').to_pylist())
RT0 = pq.read_table(f'{D}/tensor_s08_meta.parquet').column('RT0').to_numpy()
X = np.load(f'{D}/tensor_s08_traces.npy', mmap_mode='r'); M = np.load(f'{D}/tensor_s08_mask.npy')
REL = np.load(f'{D}/desc_relint.npy')
ALIAS = {'UniMod:4': 'Carbamidomethyl'}
nrm = lambda p: re.sub(r'\((.*?)\)', lambda m: '('+ALIAS.get(m.group(1), m.group(1))+')', str(p))
dn = pq.read_table('/scratch/kohlbach/xic/dn_xic.parquet', columns=['Precursor.Id','RT','Q.Value'])
q = dn.column('Q.Value').to_numpy()
TRT = {nrm(p): r*60.0 for p,r,k in zip(dn.column('Precursor.Id').to_pylist(), dn.column('RT').to_numpy(), q<=0.01) if k}

sel = [i for i in np.flatnonzero(label == 'pos')[:12000] if ids[i] in TRT]
out = {('old','lc'):[], ('old','co'):[], ('new','lc'):[], ('new','co'):[]}
outw = {k: [] for k in out}
wid = {'old':[], 'new':[]}
for a in range(0, len(sel), 400):
    s = np.array(sel[a:a+400])
    T = np.asarray(X[s], dtype=np.float32)
    for j in range(len(s)):
        i = s[j]; tr = T[j]; msk = M[i]; rel = REL[i].astype(float)
        c = int(round((TRT[ids[i]] - RT0[i]) / SP))
        if not (6 <= c < tr.shape[1]-6): continue
        tot = (tr * msk[:,None]).sum(0).astype(float); sm = smooth(tot, 1)
        # a wrong candidate: the tallest local maximum at least 20 cycles away
        far = np.arange(len(sm)); far = far[np.abs(far - c) >= 20]
        far = far[(far >= 6) & (far < len(sm)-6)]
        if far.size == 0: continue
        w = int(far[np.argmax(sm[far])])
        for tag, fn in (('old', bounds_old), ('new', bounds_new)):
            lo, hi = fn(sm, c); wid[tag].append(hi-lo+1)
            lc, co = scores(tr, msk, rel, lo, hi); out[(tag,'lc')].append(lc); out[(tag,'co')].append(co)
            lo, hi = fn(sm, w)
            lc, co = scores(tr, msk, rel, lo, hi); outw[(tag,'lc')].append(lc); outw[(tag,'co')].append(co)

print(f'{len(out[("new","lc")]):,} paired (correct apex, wrong apex) candidates\n')
print(f'{"":14}{"median width":>14}{"AUC lib_corr":>14}{"AUC coelution":>15}')
for tag in ('old','new'):
    al,_ = auc(out[(tag,'lc')], outw[(tag,'lc')])
    ac,n = auc(out[(tag,'co')], outw[(tag,'co')])
    print(f'  {tag:<12}{np.median(wid[tag]):>14.0f}{al:>14.3f}{ac:>15.3f}')
for tag in ('old','new'):
    p = np.array(out[(tag,'co')]); q_ = np.array(outw[(tag,'co')])
    m = np.isfinite(p)&np.isfinite(q_)
    print(f'  {tag}: coelution  correct {np.median(p[m]):+.3f}  wrong {np.median(q_[m]):+.3f}'
          f'  gap {np.median(p[m]-q_[m]):+.3f}')
