#!/usr/bin/env python3
"""Pull OpenSWATH's per-candidate peak groups out of an .osw, for comparison.

    scripts/osw_candidates.py <run.osw> --out cand.tsv [--precursors ids.txt]

Why OSW is the third engine. DIA-NN cannot adjudicate whether ODIA picked the
right peak: ODIA's features were ported FROM DIA-NN, its XIC is centred on
DIA-NN's own chosen apex, and its calibration is the surface that maximised its
own yield. OSW is independent of both, and it read the SAME vendor .d that
DIA-NN's XIC run read -- while ODIA reads a converted .mzpeak. So the pattern
of agreement localises a defect that no pairwise comparison can:

    ODIA differs from BOTH, and DIA-NN agrees with OSW  ->  ODIA (or its reader)
    ODIA agrees with OSW, both differ from DIA-NN       ->  DIA-NN-specific
    all three differ                                     ->  escalate to synthetic

SCORING IS NOT REQUIRED. The IH1 .osw has no SCORE_MS2 table and was set aside
as "unscored, therefore deferred". That was wrong: FEATURE and
FEATURE_TRANSITION are fully populated -- 11,453,918 features over 2,291,374
precursors, up to 5 candidates each -- and every question here is
pre-classifier. PyProphet would only add a ranking we deliberately do not use.

The join is by (modified sequence, charge), normalised through odia_join:
OSW writes ODIA's (Carbamidomethyl) dialect while DIA-NN writes (UniMod:4), and
a verbatim join silently drops 10%.
"""
import argparse, os, sqlite3, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from odia_join import to_diann_dialect

FEATURE_COLS = ['EXP_RT', 'EXP_IM', 'DELTA_RT', 'LEFT_WIDTH', 'RIGHT_WIDTH', 'NORM_RT']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('osw')
    ap.add_argument('--out', required=True)
    ap.add_argument('--precursors', help='file of Precursor.Id, one per line; '
                                         'restricts output to these (any dialect)')
    ap.add_argument('--targets-only', action='store_true', default=True)
    a = ap.parse_args()

    keep = None
    if a.precursors:
        with open(a.precursors) as fh:
            keep = {to_diann_dialect(l.strip()) for l in fh if l.strip()}
        print(f'restricting to {len(keep):,} precursors', file=sys.stderr)

    c = sqlite3.connect(f'file:{a.osw}?mode=ro', uri=True)
    have = {r[1] for r in c.execute('pragma table_info(FEATURE)')}
    cols = [x for x in FEATURE_COLS if x in have]
    print(f'FEATURE columns used: {cols}', file=sys.stderr)

    sel = ', '.join(f'f.{x}' for x in cols)
    q = f'''select pep.MODIFIED_SEQUENCE, p.CHARGE, p.DECOY, p.PRECURSOR_MZ,
                   p.LIBRARY_RT, p.LIBRARY_DRIFT_TIME, f.ID, {sel}
            from FEATURE f
            join PRECURSOR p                 on p.ID = f.PRECURSOR_ID
            join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID = p.ID
            join PEPTIDE pep                 on pep.ID = m.PEPTIDE_ID
            {'where p.DECOY = 0' if a.targets_only else ''}'''

    n_in = n_out = 0
    with open(a.out, 'w') as out:
        out.write('\t'.join(['Precursor.Id', 'Charge', 'Decoy', 'Precursor.Mz',
                             'Library.RT', 'Library.IM', 'Feature.Id'] + cols) + '\n')
        for row in c.execute(q):
            n_in += 1
            seq, charge = row[0], row[1]
            pid = to_diann_dialect(f'{seq}{charge}')
            if keep is not None and pid not in keep:
                continue
            out.write(pid + '\t' + '\t'.join('' if v is None else str(v) for v in row[1:]) + '\n')
            n_out += 1
            if n_out % 1_000_000 == 0:
                print(f'  {n_out:,} written', file=sys.stderr, flush=True)
    print(f'read {n_in:,} candidate features, wrote {n_out:,} to {a.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
