#!/usr/bin/env python3
"""Compare two engines' extracted traces on a common precursor set.

    scripts/trace_compare.py --odia stats_odia.tsv --diann stats_diann.tsv \
                             --lib matched_extract_lib.parquet [--report out.md]

Joins ODIA's product-m/z fragment key to DIA-NN's ion label through the library,
which is the only stable correspondence: DIA-NN names fragments ('y7^1'), ODIA
emits a product m/z, and the raw m/z rounded to 3 dp collides across series and
moves under calibration.

Reports, in the order a difference can be ATTRIBUTED rather than merely seen:

  1. SELECTION   which fragments each engine extracted at all -- an outer join,
                 so a fragment missing from one side is counted, not skipped.
  2. OCCUPANCY   fraction of window cycles carrying non-zero signal.
  3. SHAPE       FWHM, local maxima above 20% of apex, fraction of area inside
                 the FWHM.
  4. APEX        where each engine puts the maximum. Disagreement here is peak
                 CHOICE and contaminates every intensity comparison downstream,
                 so intensity is reported BOTH overall and restricted to
                 fragments where the apexes agree.
  5. INTENSITY   apex and area ratios, reported per precursor rather than per
                 cell: 46,000 cells from 4,000 precursors are not 46,000
                 independent observations.

Every ratio is a median of per-precursor medians for that reason.
"""
import argparse, math, statistics as st, sys, os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from odia_join import canonical_transition


def read_stats(path):
    out = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        ix = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            p = line.rstrip('\n').split('\t')
            rec = {}
            for c in ('n_pts', 'n_nz', 'n_maxima'):
                rec[c] = int(float(p[ix[c]]))
            for c in ('occupancy', 'apex_rt', 'apex_val', 'area', 'fwhm_s',
                      'frac_in_fwhm', 'baseline_frac', 'rt_lo', 'rt_hi'):
                rec[c] = float(p[ix[c]])
            out.setdefault(p[0], {})[p[1]] = rec
    return out


def lib_labels(path):
    """precursor -> [(product_mz, canonical label)], for NEAREST-m/z matching.

    Not a string key. ODIA's -out_chrom writes Product.Mz at 6 significant
    figures, so a fragment at 1086.5478 Th is written "1086.55" while the
    library holds 1086.5478 -- and an exact-string join silently drops every
    fragment above 1000 Th. That was 4,093 of 46,648 traces here, 959 of them
    peak-like, i.e. a 9% hole in the selection comparison that looked like a
    real difference between the engines.
    """
    import pyarrow.parquet as pq
    t = pq.read_table(path, columns=['Precursor.Id', 'Product.Mz', 'Fragment.Type',
                                     'Fragment.Charge', 'Fragment.Series.Number'])
    d = t.to_pydict()
    m = {}
    for pid, mzs, ty, ch, sn in zip(d['Precursor.Id'], d['Product.Mz'],
                                    d['Fragment.Type'], d['Fragment.Charge'],
                                    d['Fragment.Series.Number']):
        key = pid.replace('(Carbamidomethyl)', '(UniMod:4)')
        lst = m.setdefault(key, [])
        for mz, t_, c_, s_ in zip(mzs, ty, ch, sn):
            try:
                lst.append((mz, canonical_transition(t_, s_, c_)))
            except ValueError:
                continue
    return m


def label_for(lib, pid, mz_str, tol=0.02):
    """Nearest library fragment to an emitted product m/z, within tol Th.

    Library fragments within one precursor are separated by far more than
    0.02 Th, so this cannot pick the wrong one; and 6-significant-figure
    rounding at 1685 Th moves a value by at most 0.005 Th, so it cannot miss
    the right one.
    """
    lst = lib.get(pid)
    if not lst:
        return None
    v = float(mz_str)
    best = min(lst, key=lambda r: abs(r[0] - v))
    return best[1] if abs(best[0] - v) <= tol else None


def q(v, p):
    if not v:
        return float('nan')
    v = sorted(v)
    return v[min(len(v) - 1, max(0, int(p * len(v))))]


def describe(name, vals, fmt='{:.4g}'):
    if not vals:
        return f'{name:26s}  (none)'
    return (f'{name:26s} n={len(vals):7,}  p05 {fmt.format(q(vals,.05)):>9}  '
            f'median {fmt.format(q(vals,.50)):>9}  p95 {fmt.format(q(vals,.95)):>9}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--odia', required=True)
    ap.add_argument('--diann', required=True)
    ap.add_argument('--lib', required=True)
    ap.add_argument('--peaklike', action='store_true',
                    help='restrict to traces that look like a chromatographic peak in BOTH '
                         'engines: <=3 local maxima above 20%% of apex, >=30%% of area inside '
                         'the FWHM, and apex at least 5x the median. Without this the median '
                         'trace is noise -- 10 maxima and 15%% of area in the FWHM -- and the '
                         'comparison measures agreement between two noise fields.')
    ap.add_argument('--apex-tol', type=float, default=5.0,
                    help='seconds within which two apexes count as agreeing')
    a = ap.parse_args()

    od = read_stats(a.odia)
    dn = read_stats(a.diann)
    lab = lib_labels(a.lib)

    def peaklike(r):
        return (r['n_maxima'] <= 3 and r['frac_in_fwhm'] >= 0.30
                and r['apex_val'] > 0 and r['baseline_frac'] <= 0.20)

    common = sorted(set(od) & set(dn))
    print(f'precursors: odia {len(od):,}  diann {len(dn):,}  common {len(common):,}\n')

    # ---- 1. selection -------------------------------------------------------
    both = o_only = d_only = 0
    unmapped = 0
    pairs = []            # (pid, label, odia_rec, diann_rec)
    for pid in common:
        omap = {}
        for mz, rec in od[pid].items():
            L = label_for(lab, pid, mz)
            if L is None:
                unmapped += 1
                continue
            omap[L] = rec
        dmap = {k: v for k, v in dn[pid].items() if k not in ('ms1', 'index')}
        for L in set(omap) | set(dmap):
            if L in omap and L in dmap:
                if a.peaklike and not (peaklike(omap[L]) and peaklike(dmap[L])):
                    continue
                both += 1
                pairs.append((pid, L, omap[L], dmap[L]))
            elif L in omap:
                o_only += 1
            else:
                d_only += 1
    tot = both + o_only + d_only
    print('== 1. SELECTION (fragment level, outer join) ==')
    print(f'  extracted by BOTH        {both:8,}  ({100*both/tot:5.2f}%)')
    print(f'  ODIA only                {o_only:8,}  ({100*o_only/tot:5.2f}%)')
    print(f'  DIA-NN only              {d_only:8,}  ({100*d_only/tot:5.2f}%)')
    print(f'  ODIA fragments unmapped to a library label: {unmapped:,}')

    # ---- 2-4 ----------------------------------------------------------------
    print('\n== 2. OCCUPANCY (fraction of window cycles non-zero) ==')
    print(' ', describe('ODIA',   [p[2]['occupancy'] for p in pairs]))
    print(' ', describe('DIA-NN', [p[3]['occupancy'] for p in pairs]))
    print('\n  window length, cycles')
    print(' ', describe('ODIA',   [p[2]['n_pts'] for p in pairs], '{:.0f}'))
    print(' ', describe('DIA-NN', [p[3]['n_pts'] for p in pairs], '{:.0f}'))

    print('\n== 3. SHAPE ==')
    for key, nm in (('fwhm_s', 'FWHM (s)'), ('frac_in_fwhm', 'area inside FWHM'),
                    ('n_maxima', 'local maxima >20% apex'), ('baseline_frac', 'median/apex')):
        print(f'  {nm}')
        print('   ', describe('ODIA',   [p[2][key] for p in pairs]))
        print('   ', describe('DIA-NN', [p[3][key] for p in pairs]))

    print('\n== 4. APEX AGREEMENT ==')
    dd = [abs(p[2]['apex_rt'] - p[3]['apex_rt']) for p in pairs
          if p[2]['apex_val'] > 0 and p[3]['apex_val'] > 0]
    print(' ', describe('|apex_odia - apex_diann| s', dd))
    for tol in (1.0, 2.0, 5.0, 10.0, 20.0):
        k = sum(1 for x in dd if x <= tol)
        print(f'    within {tol:5.1f} s : {k:7,} / {len(dd):,} = {100*k/max(1,len(dd)):5.2f}%')

    # ---- 5. intensity, per precursor ---------------------------------------
    print('\n== 5. INTENSITY (median of per-precursor medians) ==')
    for restrict in (False, True):
        per = []
        for pid in common:
            r = [p for p in pairs if p[0] == pid]
            if restrict:
                r = [p for p in r if abs(p[2]['apex_rt'] - p[3]['apex_rt']) <= a.apex_tol]
            vals = [p[2]['apex_val'] / p[3]['apex_val'] for p in r
                    if p[3]['apex_val'] > 0 and p[2]['apex_val'] > 0]
            if vals:
                per.append(st.median(vals))
        tag = f'apex ratio (|dRT|<={a.apex_tol}s)' if restrict else 'apex ratio (all)'
        print(' ', describe(tag, per))


if __name__ == '__main__':
    main()
