#!/usr/bin/env python3
"""Summarise several ODIA extraction arms against one DIA-NN reference.

    scripts/trace_sweep.py --diann stats_diann.tsv --lib lib.parquet \
                           --arm 0.020=stats_odia_p12_i020.tsv --arm ...

One row per arm, so the response of each statistic to the swept parameter is
visible rather than inferred from separate runs. Restricted throughout to
fragments BOTH engines extracted, and reported both over all of those and over
the subset that looks like a chromatographic peak in both -- because ~79% of
co-extracted fragment traces are noise in both engines, and a statistic
dominated by them describes the noise field, not the extraction.
"""
import argparse, statistics as st, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trace_compare import read_stats, lib_labels, label_for


def peaklike(r):
    return (r['n_maxima'] <= 3 and r['frac_in_fwhm'] >= 0.30
            and r['apex_val'] > 0 and r['baseline_frac'] <= 0.20)


def med(v):
    return st.median(v) if v else float('nan')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--diann', required=True)
    ap.add_argument('--lib', required=True)
    ap.add_argument('--arm', action='append', required=True,
                    help='LABEL=path/to/stats_odia.tsv')
    a = ap.parse_args()
    dn = read_stats(a.diann)
    lab = lib_labels(a.lib)

    print(f'{"arm":>10} {"co-extr":>8} {"pk both":>8} {"pk O-only":>9} {"pk D-only":>9} '
          f'{"pk agree":>8} {"occ O":>6} {"occ D":>6} {"apex ratio":>10} {"|dRT| p95":>9} '
          f'{"fwhm O":>7} {"fwhm D":>7}')
    for spec in a.arm:
        label, path = spec.split('=', 1)
        od = read_stats(path)
        bb = bo = bd = nn = 0
        ratios, apexd, occo, occd, fwo, fwd = [], [], [], [], [], []
        for pid in sorted(set(od) & set(dn)):
            omap = {}
            for mz, rec in od[pid].items():
                L = label_for(lab, pid, mz)
                if L:
                    omap[L] = rec
            dmap = {k: v for k, v in dn[pid].items() if k not in ('ms1', 'index')}
            pr = []
            for L in set(omap) & set(dmap):
                o, d = omap[L], dmap[L]
                po, pd = peaklike(o), peaklike(d)
                if po and pd:
                    bb += 1
                    apexd.append(abs(o['apex_rt'] - d['apex_rt']))
                    occo.append(o['occupancy']); occd.append(d['occupancy'])
                    fwo.append(o['fwhm_s']);     fwd.append(d['fwhm_s'])
                    if d['apex_val'] > 0:
                        pr.append(o['apex_val'] / d['apex_val'])
                elif po:
                    bo += 1
                elif pd:
                    bd += 1
                else:
                    nn += 1
            if pr:
                ratios.append(st.median(pr))
        tot = bb + bo + bd + nn
        agree = 100 * bb / max(1, bb + bo + bd)
        ad = sorted(apexd)
        p95 = ad[int(0.95 * len(ad))] if ad else float('nan')
        print(f'{label:>10} {tot:8,} {bb:8,} {bo:9,} {bd:9,} {agree:7.2f}% '
              f'{med(occo):6.3f} {med(occd):6.3f} {med(ratios):10.3f} {p95:9.3f} '
              f'{med(fwo):7.3f} {med(fwd):7.3f}')


if __name__ == '__main__':
    main()
