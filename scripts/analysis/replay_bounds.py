#!/usr/bin/env python3
"""Paired replay of the OLD and NEW boundary rules over real traces.

The smallest measurement that can say whether three code changes are
net-positive, short of a full run -- which is the shape codex asked for: no
classifier, no FDR, no identifications, just the geometry and the quantities
that feed off it, on the same precursors, old rule against new.

Both rules are reimplemented here rather than called through the binary. That
is a real limitation and it is the reason this is a PROBE and not an
acceptance test: the C++ is what ships, and a Python restatement of it can
drift. It is checked against the C++ fixture values where those overlap.

Reports, for confident positives whose DIA-NN apex is known:
  * candidate width
  * whether the true apex stays inside the group
  * how often the walk ends on the span bound rather than the floor
  * background estimate: whole-trace against local flanks
"""
import sys, re, numpy as np, pyarrow.parquet as pq

D = sys.argv[1] if len(sys.argv) > 1 else '/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'
N = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
SP, MINC, MAXH, FRAC, SIG = 1.385, 7, 20, 0.10, 1.0

def smooth(v, half):
    if half <= 0: return v.copy()
    k = 2 * half + 1
    return np.convolve(v, np.ones(k) / k, mode='same')

def bounds_old(sm, apex, frac):
    """Descend while above frac*apex. No rebound guard, no span bound."""
    n = len(sm); f = frac * sm[apex]
    l = apex
    while l > 0 and sm[l - 1] > f: l -= 1
    r = apex
    while r + 1 < n and sm[r + 1] > f: r += 1
    return l, r

def bounds_new(sm, apex, frac, sigmas, minc, maxh):
    n = len(sm)
    lo = max(apex - 2, 0); hi = min(apex + 2, n - 1)
    apex = int(lo + np.argmax(sm[lo:hi + 1]))
    a = sm[apex]
    f = frac * a
    core = 6
    far = np.concatenate([sm[:max(apex - core, 0)], sm[min(apex + core + 1, n):]])
    if far.size >= 8:
        base = np.median(far)
        s = 1.4826 * np.median(np.abs(far - base))
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
    while r - l + 1 < minc and (l > 0 or r + 1 < n):
        if l > 0 and (r + 1 >= n or sm[l - 1] >= sm[r + 1]): l -= 1
        elif r + 1 < n: r += 1
        else: break
    return l, r, (g >= maxh)

lab = pq.read_table(f'{D}/tensor_ih1_labels.parquet')
ids = np.array(lab.column('Precursor.Id').to_pylist())
label = np.array(lab.column('Label').to_pylist())
meta = pq.read_table(f'{D}/tensor_ih1_meta.parquet'); RT0 = meta.column('RT0').to_numpy()
A0 = np.load(f'{D}/apex_ih1.npy')
X = np.load(f'{D}/tensor_ih1_traces.npy', mmap_mode='r')
M = np.load(f'{D}/tensor_ih1_mask.npy')

ALIAS = {'UniMod:4': 'Carbamidomethyl'}
nrm = lambda p: re.sub(r'\((.*?)\)', lambda m: '(' + ALIAS.get(m.group(1), m.group(1)) + ')', str(p))
dn = pq.read_table('/scratch/kohlbach/xic/dn_xic.parquet', columns=['Precursor.Id','RT','Q.Value'])
q = dn.column('Q.Value').to_numpy()
true_rt = {nrm(p): r * 60.0 for p, r, k in zip(dn.column('Precursor.Id').to_pylist(),
                                               dn.column('RT').to_numpy(), q <= 0.01) if k}

sel = np.flatnonzero((label == 'pos') & (A0 >= 0))[:N]
wo, wn, ino, inn, bound_hits, bg_delta = [], [], [], [], 0, []
for a in range(0, len(sel), 500):
    s = sel[a:a + 500]
    t = np.asarray(X[s], dtype=np.float32) * M[s][:, :, None]
    tot = t.sum(1)
    for k in range(len(s)):
        i = int(A0[s[k]])
        sm_old = smooth(tot[k].astype(float), 2)     # smooth_half_width
        sm_new = smooth(tot[k].astype(float), 1)     # boundary_smooth_half
        lo, hi = bounds_old(sm_old, i, FRAC)
        ln, hn, hit = bounds_new(sm_new, i, FRAC, SIG, MINC, MAXH)
        wo.append(hi - lo + 1); wn.append(hn - ln + 1)
        bound_hits += int(hit)
        r = true_rt.get(ids[s[k]])
        if r is not None:
            c = int(round((r - RT0[s[k]]) / SP))
            if 0 <= c < tot.shape[1]:
                ino.append(lo <= c <= hi); inn.append(ln <= c <= hn)
        v = tot[k].astype(float)
        out_all = np.concatenate([v[:ln], v[hn + 1:]])
        fl = np.concatenate([v[max(ln - 20, 0):ln], v[hn + 1:hn + 21]])
        if out_all.size and fl.size:
            bg_delta.append(abs(np.median(fl) - np.median(out_all)))

wo, wn = np.array(wo), np.array(wn)
print(f'{len(wo):,} confident positives replayed\n')
print(f'{"":22}{"OLD":>10}{"NEW":>10}')
for q_ in (25, 50, 75, 90):
    print(f'  width p{q_:<3} (cycles){np.percentile(wo,q_):>10.0f}{np.percentile(wn,q_):>10.0f}')
print(f'  mean width          {wo.mean():>10.1f}{wn.mean():>10.1f}')
print(f'  spans >80% window   {100*(wo>104).mean():>9.1f}%{100*(wn>104).mean():>9.1f}%')
if ino:
    print(f'\n  TRUE apex inside the group ({len(ino):,} with a known apex):')
    print(f'    old {100*np.mean(ino):5.1f}%      new {100*np.mean(inn):5.1f}%')
print(f'\n  new rule ending on the SPAN BOUND rather than the floor: '
      f'{100*bound_hits/len(wo):.1f}%')
if bg_delta:
    print(f'  |flank - whole-trace| background, new bounds: median '
          f'{np.median(bg_delta):.1f}')
