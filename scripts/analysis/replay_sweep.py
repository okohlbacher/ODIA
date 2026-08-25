#!/usr/bin/env python3
"""Why does the narrow rule miss the true apex? And is k=1.0 too aggressive?

Splits the 36.8% miss rate into its two possible causes:
  (a) the candidate apex is simply far from the true one -- a PICKING error the
      boundary rule cannot fix and must not be blamed for;
  (b) the candidate apex is close but the group is cut too narrow -- a BOUNDARY
      error, which is ours.
Then sweeps the floor to see which k buys apex coverage at what width.
"""
import sys, re, numpy as np, pyarrow.parquet as pq
D = '/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'
SP = 1.385

def smooth(v, half):
    k = 2 * half + 1
    return np.convolve(v, np.ones(k) / k, mode='same')

def bounds(sm, apex, frac, sigmas, minc, maxh):
    n = len(sm)
    lo = max(apex - 2, 0); hi = min(apex + 2, n - 1)
    apex = int(lo + np.argmax(sm[lo:hi + 1])); a = sm[apex]
    f = frac * a; core = 6
    far = np.concatenate([sm[:max(apex - core, 0)], sm[min(apex + core + 1, n):]])
    if far.size >= 8:
        base = np.median(far); s = 1.4826 * np.median(np.abs(far - base))
        f = max(f, min(base + sigmas * s, a - 1e-9 * max(1.0, abs(a))))
    reb = 0.25 * a
    l = apex; run = a; g = 0
    while l > 0 and g < maxh:
        v = sm[l - 1]
        if v <= f or v > run + reb: break
        run = min(run, v); l -= 1; g += 1
    r = apex; run = a; g = 0
    while r + 1 < n and g < maxh:
        v = sm[r + 1]
        if v <= f or v > run + reb: break
        run = min(run, v); r += 1; g += 1
    walked = r - l + 1
    while r - l + 1 < minc and (l > 0 or r + 1 < n):
        if l > 0 and (r + 1 >= n or sm[l - 1] >= sm[r + 1]): l -= 1
        elif r + 1 < n: r += 1
        else: break
    return l, r, walked, apex

lab = pq.read_table(f'{D}/tensor_s08_labels.parquet')
ids = np.array(lab.column('Precursor.Id').to_pylist()); label = np.array(lab.column('Label').to_pylist())
RT0 = pq.read_table(f'{D}/tensor_s08_meta.parquet').column('RT0').to_numpy()
A0 = np.load(f'{D}/apex_s08.npy'); X = np.load(f'{D}/tensor_s08_traces.npy', mmap_mode='r'); M = np.load(f'{D}/tensor_s08_mask.npy')
ALIAS = {'UniMod:4': 'Carbamidomethyl'}
nrm = lambda p: re.sub(r'\((.*?)\)', lambda m: '(' + ALIAS.get(m.group(1), m.group(1)) + ')', str(p))
dn = pq.read_table('/scratch/kohlbach/xic/dn_xic.parquet', columns=['Precursor.Id','RT','Q.Value'])
q = dn.column('Q.Value').to_numpy()
true_rt = {nrm(p): r*60.0 for p, r, k in zip(dn.column('Precursor.Id').to_pylist(), dn.column('RT').to_numpy(), q <= 0.01) if k}

sel = np.flatnonzero((label == 'pos') & (A0 >= 0))[:4000]
KS = [0.0, 0.25, 0.5, 1.0, 2.0]
res = {k: {'w': [], 'in': [], 'walk': []} for k in KS}
dist, fwhm = [], []
for a in range(0, len(sel), 500):
    s = sel[a:a+500]
    t = np.asarray(X[s], dtype=np.float32) * M[s][:, :, None]; tot = t.sum(1)
    for j in range(len(s)):
        i = int(A0[s[j]]); v = tot[j].astype(float); sm = smooth(v, 1)
        r = true_rt.get(ids[s[j]])
        c = int(round((r - RT0[s[j]]) / SP)) if r is not None else None
        if c is not None and not (0 <= c < len(v)): c = None
        # empirical FWHM around the candidate apex, over the local baseline
        ap = sm[i]; base = np.median(sm); h = base + 0.5*(ap-base)
        l_ = i
        while l_ > 0 and sm[l_-1] > h: l_ -= 1
        r_ = i
        while r_+1 < len(sm) and sm[r_+1] > h: r_ += 1
        fwhm.append(r_-l_+1)
        if c is not None: dist.append(abs(i - c))
        for k in KS:
            lo, hi, walked, ap2 = bounds(sm, i, 0.10, k, 7, 20)
            res[k]['w'].append(hi-lo+1); res[k]['walk'].append(walked)
            if c is not None: res[k]['in'].append(lo <= c <= hi)

dist = np.array(dist); fwhm = np.array(fwhm)
print(f'{len(dist):,} positives with a known true apex\n')
print('CANDIDATE APEX ERROR |A0 - true| in cycles  (a boundary rule cannot fix this)')
for p in (50, 75, 90, 95): print(f'    p{p:<3}{np.percentile(dist,p):>6.0f}')
for d in (2, 3, 5, 10, 20):
    print(f'    within {d:>2} cycles: {100*(dist<=d).mean():5.1f}%')
print(f'\nEMPIRICAL FWHM at the candidate apex (cycles): p25 {np.percentile(fwhm,25):.0f} '
      f'p50 {np.percentile(fwhm,50):.0f} p75 {np.percentile(fwhm,75):.0f} p90 {np.percentile(fwhm,90):.0f}')
print(f'\n{"k":>5}{"med width":>11}{"p90 w":>8}{"walked=1":>10}{"apex in":>9}'
      f'{"in|d<=3":>9}')
for k in KS:
    w = np.array(res[k]['w']); wk = np.array(res[k]['walk']); ii = np.array(res[k]['in'])
    near = ii[dist <= 3]
    print(f'{k:>5.2f}{np.median(w):>11.0f}{np.percentile(w,90):>8.0f}'
          f'{100*(wk<=1).mean():>9.1f}%{100*ii.mean():>8.1f}%{100*near.mean():>8.1f}%')
