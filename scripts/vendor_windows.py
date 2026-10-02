#!/usr/bin/env python3
"""The diaPASEF acquisition geometry, straight from the vendor .d.

    scripts/vendor_windows.py <run.d> [--tsv out.tsv]

Why this exists. Both ODIA and DIA-NN *infer* the isolation scheme -- ODIA reads
it out of the converted spectra rather than the stated method, and neither side
wrote it anywhere, so a misread scheme (and on diaPASEF, a collapsed m/z x 1/K0
geometry) was invisible in both directions. Bruker states the scheme exactly, in
plain SQLite inside the .d, so there is no reason for either engine's guess to
be the reference.

`DiaFrameMsMsWindows` is authoritative: one row per (window group, mobility
scan range), with the isolation m/z and width. On IH1 that is 12 groups x 2
windows = 24 windows, isolation widths spanning 23.0 to 266.6 Th -- so a
comparison that assumes uniform windows is checking something the instrument
never did.

Mobility is reported in SCAN NUMBERS, which decrease in 1/K0 as they increase.
Converting to 1/K0 needs the TimsCalibration, which this script does not do;
scan bounds are compared against scan bounds, and any engine reporting 1/K0
must be converted to scans, not the other way round, because the scan grid is
the thing the instrument actually sampled.

Read-only: the .d lives under data/ and is never written.
"""
import sqlite3, sys, os, json


def read(dpath):
    tdf = os.path.join(dpath, 'analysis.tdf')
    if not os.path.exists(tdf):
        raise SystemExit(f'no analysis.tdf in {dpath} -- is this a Bruker .d?')
    c = sqlite3.connect(f'file:{tdf}?mode=ro', uri=True)
    win = [dict(zip(('group', 'scan_begin', 'scan_end', 'iso_mz', 'iso_width', 'ce'), r))
           for r in c.execute(
               'select WindowGroup,ScanNumBegin,ScanNumEnd,IsolationMz,IsolationWidth,'
               'CollisionEnergy from DiaFrameMsMsWindows order by WindowGroup,ScanNumBegin')]
    for w in win:
        w['mz_low'] = w['iso_mz'] - w['iso_width'] / 2.0
        w['mz_high'] = w['iso_mz'] + w['iso_width'] / 2.0
    q = lambda s: c.execute(s).fetchone()
    meta = {
        'window_groups': q('select count(*) from DiaFrameMsMsWindowGroups')[0],
        'windows': len(win),
        'ms1_frames': q('select count(*) from Frames where MsMsType=0')[0],
        'ms2_frames': q('select count(*) from Frames where MsMsType=9')[0],
        'num_scans_per_frame': [r[0] for r in c.execute('select distinct NumScans from Frames')],
        'rt_min_s': q('select min(Time) from Frames')[0],
        'rt_max_s': q('select max(Time) from Frames')[0],
        'mz_low': min(w['mz_low'] for w in win),
        'mz_high': max(w['mz_high'] for w in win),
        'scan_min': min(w['scan_begin'] for w in win),
        'scan_max': max(w['scan_end'] for w in win),
    }
    return win, meta


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    dpath = sys.argv[1]
    win, meta = read(dpath)
    print(json.dumps(meta, indent=2))
    print()
    print(' grp scanBeg scanEnd   isolMz  width     mz_low    mz_high     CE')
    for w in win:
        print('%4d %7d %7d %8.2f %6.1f %10.2f %10.2f %6.1f'
              % (w['group'], w['scan_begin'], w['scan_end'], w['iso_mz'],
                 w['iso_width'], w['mz_low'], w['mz_high'], w['ce']))
    if '--tsv' in sys.argv:
        out = sys.argv[sys.argv.index('--tsv') + 1]
        with open(out, 'w') as fh:
            fh.write('group\tscan_begin\tscan_end\tiso_mz\tiso_width\tmz_low\tmz_high\tce\n')
            for w in win:
                fh.write('%d\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.2f\n'
                         % (w['group'], w['scan_begin'], w['scan_end'], w['iso_mz'],
                            w['iso_width'], w['mz_low'], w['mz_high'], w['ce']))
        print(f'\nwrote {len(win)} windows to {out}')


if __name__ == '__main__':
    main()
