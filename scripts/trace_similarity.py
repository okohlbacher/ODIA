#!/usr/bin/env python3
"""Pointwise similarity and lag between corresponding ODIA and DIA-NN traces.

    scripts/trace_similarity.py --odia out_chrom.tsv --diann dn_xic.parquet \
        --lib lib.parquet --precursors ids.txt [--maxlag 4] [--out sim.tsv]

Summary statistics can agree while the traces disagree: two engines can produce
the same apex, width and area from different signal. This asks the question that
cannot be answered that way -- cycle by cycle, is it the SAME trace?

Correspondence is established in three steps, each of which can fail loudly:

  1. precursor   the shared list, one dialect (UniMod:4).
  2. transition  ODIA's product m/z -> the library's canonical (type, ordinal,
                 charge) -> DIA-NN's ion label. Nearest m/z within 0.02 Th,
                 because ODIA writes Product.Mz at 6 significant figures and an
                 exact key drops everything above 1000 Th.
  3. cycle       nearest retention time within half a cycle. The grids are
                 checked for identity first and the check is reported, because
                 if they ARE identical then every lag below is a real shift and
                 not a resampling artefact.

Reported per trace pair:
  pearson        on the overlapping cycles, raw intensities.
  cosine         scale-free, so a pure intensity offset does not read as
                 disagreement.
  best_lag       the cycle shift maximising Pearson, searched over +/-maxlag.
                 A systematic non-zero lag is an RT offset between the engines;
                 a symmetric spread around zero is noise.
  r_at_best      Pearson at that lag, so the gain from shifting is visible.

THE NULL IS NOT ZERO. Two traces of the same precursor in the same window share
a baseline and an elution envelope, so a WRONG fragment of the SAME precursor
already correlates highly -- measured here, around 0.97 on peak-like traces.
Reading r against zero therefore overstates agreement enormously. Every
correlation below is reported beside that null, and the margin between them is
the statistic that carries information.

AND PEARSON CANNOT SEE AMPLITUDE. Pearson is affine-invariant and cosine is
scale-invariant, so neither can detect that one engine returns systematically
more signal than the other. The area ratio is reported alongside for that
reason; it is the only number here that moves with aperture.

Traces flat in either engine are excluded and counted: Pearson is undefined
there, and including them as zero would report a disagreement that is really an
absence.
"""
import argparse, math, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from odia_join import canonical_transition


def pearson(a, b):
    n = len(a)
    if n < 3:
        return float('nan')
    ma, mb = sum(a) / n, sum(b) / n
    va = sum((x - ma) ** 2 for x in a)
    vb = sum((x - mb) ** 2 for x in b)
    if va <= 0 or vb <= 0:
        return float('nan')
    cov = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    return cov / math.sqrt(va * vb)


def cosine(a, b):
    na = math.sqrt(sum(x * x for x in a))
    nb = math.sqrt(sum(y * y for y in b))
    if na <= 0 or nb <= 0:
        return float('nan')
    return sum(x * y for x, y in zip(a, b)) / (na * nb)


def load_lib(path):
    import pyarrow.parquet as pq
    d = pq.read_table(path, columns=['Precursor.Id', 'Product.Mz', 'Fragment.Type',
                                     'Fragment.Charge', 'Fragment.Series.Number']).to_pydict()
    m = {}
    for pid, mzs, ty, ch, sn in zip(d['Precursor.Id'], d['Product.Mz'], d['Fragment.Type'],
                                    d['Fragment.Charge'], d['Fragment.Series.Number']):
        k = pid.replace('(Carbamidomethyl)', '(UniMod:4)')
        lst = m.setdefault(k, [])
        for mz, a, b, c in zip(mzs, ty, ch, sn):
            try:
                lst.append((mz, canonical_transition(a, c, b)))
            except ValueError:
                pass
    return m


def label_for(lib, pid, mz, tol=0.02):
    lst = lib.get(pid)
    if not lst:
        return None
    best = min(lst, key=lambda r: abs(r[0] - mz))
    return best[1] if abs(best[0] - mz) <= tol else None


def load_odia(path, keep, lib):
    out = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        ix = {c: i for i, c in enumerate(hdr)}
        dec = ix.get('Decoy')
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if dec is not None and p[dec] not in ('0', 'false', 'False'):
                continue
            pid = p[ix['Precursor.Id']].replace('(Carbamidomethyl)', '(UniMod:4)')
            if pid not in keep:
                continue
            lab = label_for(lib, pid, float(p[ix['Product.Mz']]))
            if lab is None:
                continue
            out.setdefault((pid, lab), ([], []))
            out[(pid, lab)][0].append(float(p[ix['RT']]))
            out[(pid, lab)][1].append(float(p[ix['Intensity']]))
    return out


def load_diann(path, keep):
    import pyarrow.parquet as pq
    f = pq.ParquetFile(path)
    out = {}
    for b in f.iter_batches(batch_size=1_000_000, columns=['pr', 'feature', 'rt', 'value']):
        d = b.to_pydict()
        for pr, fe, rt, va in zip(d['pr'], d['feature'], d['rt'], d['value']):
            if pr not in keep or fe in ('ms1', 'index'):
                continue
            out.setdefault((pr, fe), ([], []))
            out[(pr, fe)][0].append(rt * 60.0)   # minutes -> seconds
            out[(pr, fe)][1].append(va)
    return out


def align(ort, oval, drt, dval, tol):
    """Nearest-cycle join. One-to-one: each DIA-NN cycle consumes one ODIA cycle."""
    a, b = [], []
    i = 0
    used = -1
    for t, v in zip(drt, dval):
        while i + 1 < len(ort) and abs(ort[i + 1] - t) <= abs(ort[i] - t):
            i += 1
        if i != used and abs(ort[i] - t) <= tol:
            a.append(oval[i]); b.append(v); used = i
    return a, b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--odia', required=True)
    ap.add_argument('--diann', required=True)
    ap.add_argument('--lib', required=True)
    ap.add_argument('--precursors', required=True)
    ap.add_argument('--maxlag', type=int, default=4)
    ap.add_argument('--cycle', type=float, default=1.385)
    ap.add_argument('--out')
    a = ap.parse_args()

    keep = {l.strip() for l in open(a.precursors) if l.strip()}
    lib = load_lib(a.lib)
    print(f'loading ODIA traces for {len(keep):,} precursors...', file=sys.stderr)
    od = load_odia(a.odia, keep, lib)
    print(f'  {len(od):,} ODIA traces', file=sys.stderr)
    dn = load_diann(a.diann, keep)
    print(f'  {len(dn):,} DIA-NN traces', file=sys.stderr)

    common = sorted(set(od) & set(dn))
    print(f'  {len(common):,} corresponding traces\n', file=sys.stderr)

    # Are the two cycle grids the same grid?
    gd = []
    for k in common[:400]:
        o, d = od[k][0], dn[k][0]
        for t in d:
            gd.append(min(abs(t - x) for x in o))
    gd.sort()
    print(f'GRID: |odia_rt - diann_rt| over {len(gd):,} cycles of 400 traces: '
          f'median {gd[len(gd)//2]:.4f}s  p95 {gd[int(.95*len(gd))]:.4f}s  max {gd[-1]:.4f}s')
    print(f'      cycle is {a.cycle:.3f}s, so the grids are '
          f'{"IDENTICAL" if gd[int(.95*len(gd))] < 0.01 else "DIFFERENT"}\n')

    # Wrong-fragment null: the same ODIA trace against a DIFFERENT fragment of the
    # SAME precursor, in the same window. This is the floor any true correlation
    # must clear, and it is nowhere near zero.
    by_prec = {}
    for (p_, f_) in common:
        by_prec.setdefault(p_, []).append(f_)

    rows = []
    flat = 0
    for k in common:
        ort, oval = od[k]
        drt, dval = dn[k]
        o, d = align(ort, oval, drt, dval, a.cycle / 2.0)
        if len(o) < 10:
            continue
        r = pearson(o, d)
        if math.isnan(r):
            flat += 1
            continue
        cs = cosine(o, d)
        best_r, best_lag = r, 0
        for L in range(-a.maxlag, a.maxlag + 1):
            if L == 0:
                continue
            if L > 0:
                x, y = o[L:], d[:len(d) - L]
            else:
                x, y = o[:len(o) + L], d[-L:]
            if len(x) < 10:
                continue
            rr = pearson(x, y)
            if not math.isnan(rr) and rr > best_r:
                best_r, best_lag = rr, L
        # area ratio -- the quantity the two correlation metrics are blind to
        so, sd = sum(o), sum(d)
        area = (so / sd) if sd > 0 else float('nan')

        sib = [f_ for f_ in by_prec[k[0]] if f_ != k[1]]
        null_r = float('nan')
        if sib:
            # best-of-siblings, the hardest null: if a wrong fragment can match
            # ODIA better than the right one, the metric is not identifying.
            best_null = -2.0
            for f_ in sib:
                ort2, oval2 = od[(k[0], f_)]
                o2, d2 = align(ort2, oval2, drt, dval, a.cycle / 2.0)
                if len(o2) < 10:
                    continue
                rr = pearson(o2, d2)
                if not math.isnan(rr) and rr > best_null:
                    best_null = rr
            if best_null > -2.0:
                null_r = best_null
        rows.append((k[0], k[1], len(o), r, cs, best_lag, best_r, area, null_r))

    print(f'compared {len(rows):,} traces ({flat:,} flat in one engine, excluded)\n')

    def dist(nm, vals, fmt='{:+.4f}'):
        v = sorted(vals)
        n = len(v)
        print(f'{nm:22s} p05 {fmt.format(v[int(.05*n)])}  p25 {fmt.format(v[int(.25*n)])}  '
              f'median {fmt.format(v[n//2])}  p75 {fmt.format(v[int(.75*n)])}  '
              f'p95 {fmt.format(v[int(.95*n)])}')

    dist('pearson (lag 0)', [r[3] for r in rows])
    dist('cosine  (lag 0)', [r[4] for r in rows])
    nulls = [r[8] for r in rows if not math.isnan(r[8])]
    dist('WRONG-FRAGMENT null', nulls)
    marg = [r[3] - r[8] for r in rows if not math.isnan(r[8])]
    dist('margin over null', marg)
    beat = sum(1 for r in rows if not math.isnan(r[8]) and r[3] > r[8])
    print(f'{"":22s} true beats null in {beat:,}/{len(nulls):,} = {100*beat/max(1,len(nulls)):.2f}%')
    dist('AREA ratio odia/diann', [r[7] for r in rows if not math.isnan(r[7])], '{:.4f}')
    print()
    for thr in (0.9, 0.8, 0.5, 0.2, 0.0):
        k = sum(1 for r in rows if r[3] >= thr)
        print(f'  pearson >= {thr:4.1f} : {k:7,} / {len(rows):,} = {100*k/len(rows):5.2f}%')
    print('\nLAG (cycles; +1 means ODIA must shift LATER to match DIA-NN)')
    from collections import Counter
    c = Counter(r[5] for r in rows)
    for L in sorted(c):
        print(f'  {L:+d} : {c[L]:7,} ({100*c[L]/len(rows):5.2f}%)  '
              f'= {L*a.cycle:+.3f}s')
    gain = [r[6] - r[3] for r in rows]
    dist('pearson gain from lag', gain)

    if a.out:
        with open(a.out, 'w') as fh:
            fh.write('Precursor.Id\tfragment\tn\tpearson\tcosine\tbest_lag\tr_at_best'
                     '\tarea_ratio\tnull_r\n')
            for r in rows:
                fh.write(f'{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.6f}\t{r[4]:.6f}\t{r[5]}\t{r[6]:.6f}'
                         f'\t{r[7]:.6f}\t{r[8]:.6f}\n')
        print(f'\nwrote {len(rows):,} rows to {a.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
