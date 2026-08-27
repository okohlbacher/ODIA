#!/usr/bin/env python3
"""Characterise extracted chromatograms, per trace, for cross-engine comparison.

    scripts/trace_stats.py diann  <dn_xic.parquet>   --out stats.tsv [--sample N]
    scripts/trace_stats.py odia   <out_chrom.tsv>    --out stats.tsv [--sample N]

One row per (precursor, fragment). The metrics are chosen so that a difference
points at a CAUSE rather than merely existing:

  n_pts, n_nz, occupancy   how much of the window carries signal at all. An
                           engine that extracts the same peak into a wider
                           window has lower occupancy without being wrong; an
                           engine that misses the peak has near-zero n_nz.
  apex_rt, apex_val        where the maximum sits. Disagreement here is peak
                           CHOICE, which no intensity ratio can separate from
                           extraction quality.
  fwhm_s                   width at half maximum, in seconds. Sensitive to
                           smoothing, to window aperture, and to co-eluting
                           interference merging two peaks into one.
  area                     summed intensity over the window.
  frac_in_fwhm             area within the FWHM divided by total area. A clean
                           single peak concentrates its area; an interfered or
                           noisy trace spreads it. This is the closest thing
                           here to a purity measure that needs no reference.
  n_maxima                 local maxima above 20% of the apex. >1 means the
                           window contains more than one candidate peak, which
                           is what makes a peak-picking disagreement possible.
  baseline_frac            median divided by apex -- a high value is a trace
                           with no peak, only floor.

Nothing is smoothed here. Smoothing is a controlled variable in the comparison,
not a property of the measurement.
"""
import argparse, sys, math, random

SMOOTH = False


def diann_smooth(v):
    """DIA-NN's smooth(), verbatim: diann_1.7.12.cpp:1035.

    A single-pass binomial [0.25, 0.5, 0.25] with mirrored 2/3-1/3 ends -- NOT a
    moving average, which would be [1/3,1/3,1/3] and different at the edges.
    Applied here as a CONTROL: if ODIA's traces reproduce DIA-NN's statistics
    once this is applied, smoothing was the difference and nothing else was.

    Read from 1.7.12 source while the compared file came from the closed 2.0
    binary, so this is the best available control, not a proven equivalence.
    """
    n = len(v)
    if n < 2:
        return list(v)
    o = [0.0] * n
    o[0] = (2.0 / 3.0) * v[0] + (1.0 / 3.0) * v[1]
    o[n - 1] = (2.0 / 3.0) * v[n - 1] + (1.0 / 3.0) * v[n - 2]
    for i in range(1, n - 1):
        o[i] = 0.5 * v[i] + 0.25 * (v[i - 1] + v[i + 1])
    return o


def trace_metrics(rts, vals):
    if SMOOTH:
        vals = diann_smooth(vals)
    n = len(vals)
    if n == 0:
        return None
    nz = [v for v in vals if v > 0]
    apex_i = max(range(n), key=lambda i: vals[i])
    apex = vals[apex_i]
    total = sum(vals)
    out = {
        'n_pts': n, 'n_nz': len(nz),
        'occupancy': len(nz) / n,
        'apex_rt': rts[apex_i], 'apex_val': apex,
        'area': total,
        'rt_lo': rts[0], 'rt_hi': rts[-1],
    }
    if apex <= 0:
        out.update({'fwhm_s': 0.0, 'frac_in_fwhm': 0.0, 'n_maxima': 0,
                    'baseline_frac': 0.0})
        return out
    half = apex / 2.0
    lo = apex_i
    while lo > 0 and vals[lo - 1] >= half:
        lo -= 1
    hi = apex_i
    while hi < n - 1 and vals[hi + 1] >= half:
        hi += 1
    out['fwhm_s'] = abs(rts[hi] - rts[lo])
    seg = sum(vals[lo:hi + 1])
    out['frac_in_fwhm'] = seg / total if total > 0 else 0.0
    thr = 0.2 * apex
    m = 0
    for i in range(1, n - 1):
        if vals[i] > thr and vals[i] >= vals[i - 1] and vals[i] > vals[i + 1]:
            m += 1
    out['n_maxima'] = max(m, 1)
    s = sorted(vals)
    out['baseline_frac'] = s[n // 2] / apex
    return out


COLS = ['Precursor.Id', 'fragment', 'n_pts', 'n_nz', 'occupancy', 'apex_rt',
        'apex_val', 'area', 'fwhm_s', 'frac_in_fwhm', 'n_maxima',
        'baseline_frac', 'rt_lo', 'rt_hi']


def _norm(pid):
    return pid.replace('(Carbamidomethyl)', '(UniMod:4)')


def emit(out, pid, frag, m):
    pid = _norm(pid)
    out.write('\t'.join([pid, frag] + [
        (f'{m[c]:.6g}' if isinstance(m[c], float) else str(m[c])) for c in COLS[2:]]) + '\n')


def run_diann(path, outp, sample, seed, explicit=None):
    import pyarrow.parquet as pq
    f = pq.ParquetFile(path)
    keep = explicit
    if sample and keep is None:
        ids = set()
        for b in f.iter_batches(batch_size=2_000_000, columns=['pr']):
            ids |= set(b.column('pr').to_pylist())
        rnd = random.Random(seed)
        keep = set(rnd.sample(sorted(ids), min(sample, len(ids))))
        print(f'sampling {len(keep):,} of {len(ids):,} precursors', file=sys.stderr)
    cur = None
    rts, vals = [], []
    n = 0
    with open(outp, 'w') as out:
        out.write('\t'.join(COLS) + '\n')
        for b in f.iter_batches(batch_size=1_000_000,
                                columns=['pr', 'feature', 'rt', 'value']):
            d = b.to_pydict()
            for pr, fe, rt, va in zip(d['pr'], d['feature'], d['rt'], d['value']):
                if keep is not None and pr not in keep:
                    continue
                k = (pr, fe)
                if k != cur:
                    if cur is not None:
                        m = trace_metrics(rts, vals)
                        if m:
                            emit(out, cur[0], cur[1], m); n += 1
                    cur, rts, vals = k, [], []
                rts.append(rt * 60.0)   # DIA-NN rt is MINUTES; everything here is seconds
                vals.append(va)
        if cur is not None:
            m = trace_metrics(rts, vals)
            if m:
                emit(out, cur[0], cur[1], m); n += 1
    print(f'wrote {n:,} traces to {outp}', file=sys.stderr)


def run_odia(path, outp, sample, seed, explicit=None):
    """ODIA -out_chrom: one row per (precursor, transition, cycle)."""
    hdr = None
    cur = None
    rts, vals = [], []
    n = 0
    keep = explicit
    if sample and keep is None:
        ids = set()
        with open(path) as fh:
            h = fh.readline().rstrip('\n').split('\t')
            ip = h.index('Precursor.Id')
            for line in fh:
                ids.add(line.split('\t', ip + 1)[ip])
        rnd = random.Random(seed)
        keep = set(rnd.sample(sorted(ids), min(sample, len(ids))))
        print(f'sampling {len(keep):,} of {len(ids):,} precursors', file=sys.stderr)
    with open(path) as fh, open(outp, 'w') as out:
        out.write('\t'.join(COLS) + '\n')
        hdr = fh.readline().rstrip('\n').split('\t')
        ix = {c: i for i, c in enumerate(hdr)}
        need = ['Precursor.Id', 'RT', 'Intensity']
        # Product m/z, not Transition.Index: the index is the GLOBAL library
        # index (transition_begin[i]+k), not a per-precursor slot, so it cannot
        # be joined to another engine's fragment without the offsets. The m/z
        # maps to a canonical (type, ordinal, charge) through the library.
        fragcol = 'Product.Mz' if 'Product.Mz' in ix else (
            'Transition.Index' if 'Transition.Index' in ix else None)
        for c in need:
            if c not in ix:
                raise SystemExit(f'{path}: missing column {c}; has {hdr}')
        dec = ix.get('Decoy')
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if dec is not None and p[dec] not in ('0', 'false', 'False'):
                continue
            pid = p[ix['Precursor.Id']]
            if keep is not None and _norm(pid) not in keep:
                continue
            frag = p[ix[fragcol]] if fragcol else '?'
            if fragcol == 'Product.Mz':
                frag = f'{float(frag):.3f}'
            k = (pid, frag)
            if k != cur:
                if cur is not None:
                    m = trace_metrics(rts, vals)
                    if m:
                        emit(out, cur[0], cur[1], m); n += 1
                cur, rts, vals = k, [], []
            rts.append(float(p[ix['RT']]))
            vals.append(float(p[ix['Intensity']]))
        if cur is not None:
            m = trace_metrics(rts, vals)
            if m:
                emit(out, cur[0], cur[1], m); n += 1
    print(f'wrote {n:,} traces to {outp}', file=sys.stderr)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('engine', choices=['diann', 'odia'])
    ap.add_argument('path')
    ap.add_argument('--out', required=True)
    ap.add_argument('--sample', type=int, default=0)
    ap.add_argument('--precursors', help='file of Precursor.Id, one per line -- the SAME '
                                         'list for both engines. Overrides --sample, which '
                                         'draws from each engine\'s own pool and so cannot '
                                         'produce a comparable pair.')
    ap.add_argument('--seed', type=int, default=20260827)
    ap.add_argument('--smooth', action='store_true',
                    help="apply DIA-NN's [0.25,0.5,0.25] kernel before measuring")
    a = ap.parse_args()
    explicit = None
    if a.precursors:
        with open(a.precursors) as fh:
            explicit = {l.strip() for l in fh if l.strip()}
        print(f'restricted to {len(explicit):,} listed precursors', file=sys.stderr)
    SMOOTH = a.smooth
    globals()['SMOOTH'] = SMOOTH
    (run_diann if a.engine == 'diann' else run_odia)(a.path, a.out, a.sample, a.seed, explicit)
