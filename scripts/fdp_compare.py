#!/usr/bin/env python3
"""Compare runs at MATCHED EMPIRICAL FDP. The acceptance instrument.

    scripts/fdp_compare.py base=rank_a.tsv arm=rank_b.tsv [...]

Why this exists. Nominal q on this data delivers 5-7x its claim (doc/46), and
the error is 109x at depth 2,000 against 1.4x at depth 20,000 -- so it is not
even a constant. Two arms at "nominal 1%" therefore sit at DIFFERENT empirical
FDPs and comparing them there compares two different thresholds. That mistake
reversed the Gate C verdict once (doc/46: -17% at nominal, +12% at matched).

What it reports, per doc/49's control 2 and 3:
  * identifications at a grid of MATCHED empirical FDP targets;
  * the Poisson error on every FDP, because ~40-130 entrapment hits is the whole
    story and differences inside it are not differences;
  * recovery of DIA-NN's confident set at each depth, as cross-engine
    concordance -- NOT as truth, since that set is itself ~7% false.

Entrapment FDP is (e/r)/t over ordinary targets. r is computed on the library
the arm ACTUALLY SEARCHED: the fragment floor changes it, and using the wrong
one shifts the number ~0.4%. The rule used is printed so it can be checked.
"""
import re, sys, os
import pyarrow.parquet as pq, pyarrow.compute as pc

L = os.environ.get('ODIA_LIBV2', '/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/libv2')
args = [a for a in sys.argv[1:] if '=' in a]
if not args:
    sys.exit(__doc__)

arab = set(l.strip() for l in open(f'{L}/arab_acc.txt') if l.strip())
t = pq.read_table(f'{L}/human_v2.parquet',
                  columns=['Precursor.Id', 'Protein.Group', 'Decoy', 'Product.Mz'])
nfr = pc.list_value_length(t.column('Product.Mz')).to_numpy()
ent, nf = {}, {}
for i, (p, g, d) in enumerate(zip(t.column('Precursor.Id').to_pylist(),
                                  t.column('Protein.Group').to_pylist(),
                                  t.column('Decoy').to_numpy())):
    if d:
        continue
    ent[p] = any(a in arab for a in str(g).replace(';', '|').split('|') if a)
    nf[p] = int(nfr[i])

# DIA-NN writes C(UniMod:4); our library writes C(Carbamidomethyl). It is the
# only modification on either side, so matching Precursor.Id verbatim drops
# EVERY cysteine peptide -- 3,902 of DIA-NN's 39,149 confident precursors on
# IH1, 10.0% -- and understates concordance by that much. Alias, do not strip:
# stripping the parentheses would also merge modified and unmodified forms.
DN_ALIAS = {'UniMod:4': 'Carbamidomethyl'}


def dn_norm(pid):
    return re.sub(r'\((.*?)\)',
                  lambda m: '(' + DN_ALIAS.get(m.group(1), m.group(1)) + ')',
                  str(pid))


DN = set()
try:
    rr = pq.read_table(f'{L}/dn_xic_report.parquet',
                       columns=['Precursor.Id', 'Q.Value']).to_pydict()
    DN = {dn_norm(p) for p, q in zip(rr['Precursor.Id'], rr['Q.Value'])
          if q is None or q <= 0.01}
except Exception:
    pass


def ratio(floored):
    keep = (lambda p: nf[p] >= 3) if floored else (lambda p: True)
    E = sum(1 for p, e in ent.items() if e and keep(p))
    T = sum(1 for p, e in ent.items() if not e and keep(p))
    return E / T


def curve(path):
    rows, acc_short = [], False
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if f[0] not in ent:
            continue
        try:
            s, q = float(f[1]), float(f[2])
        except (ValueError, IndexError):
            continue
        if q <= 0.01 and nf[f[0]] < 3:
            acc_short = True
        rows.append((s, ent[f[0]], f[0]))
    r = ratio(not acc_short)
    rows.sort(key=lambda x: -x[0])
    he = ht = 0
    out = []
    for _, e, p in rows:
        he += e
        ht += (not e)
        out.append((ht, (he / r) / max(ht, 1), he, p))
    return out, r, (not acc_short)


C = {}
for a in args:
    tag, path = a.split('=', 1)
    C[tag] = curve(path)
    _, r, floored = C[tag]
    print(f"  {tag:>14}: r = {r:.5f}  (library {'floored at 3 fragments' if floored else 'unfiltered'})")

TARGETS = [0.02, 0.03, 0.05, 0.0572, 0.0742, 0.10, 0.15]
tags = list(C)
print(f"\n{'FDP target':>11}" + ''.join(f'{t:>22}' for t in tags))
print(f"{'':>11}" + ''.join(f"{'IDs   (+-)   DIA-NN':>22}" for _ in tags))
for tgt in TARGETS:
    cells = []
    for tag in tags:
        rows, _, _ = C[tag]
        best = idx = 0
        for i, (ht, f, he, _) in enumerate(rows):
            if f <= tgt and ht > best:
                best, idx = ht, i
        if best == 0:
            cells.append(f"{'--':>22}")
            continue
        he = rows[idx][2]
        poi = 100 * tgt / max(he, 1) ** 0.5
        rec = sum(1 for _, _, _, p in rows[:idx + 1] if p in DN)
        cells.append(f'{best:>9,}{poi:>7.2f}{rec:>6,}')
    star = '  <- DIA-NN operating point' if abs(tgt - 0.0742) < 1e-9 else ''
    print(f'{100*tgt:>10.2f}%' + ''.join(cells) + star)

if len(tags) == 2:
    a, b = tags
    print(f'\n  {b} vs {a}, at each matched FDP:')
    for tgt in TARGETS:
        v = []
        for tag in (a, b):
            rows, _, _ = C[tag]
            best = idx = 0
            for i, (ht, f, he, _) in enumerate(rows):
                if f <= tgt and ht > best:
                    best, idx = ht, i
            v.append((best, rows[idx][2] if best else 0))
        if not v[0][0] or not v[1][0]:
            continue
        d = v[1][0] - v[0][0]
        sig = (100 * tgt) * ((1 / max(v[0][1], 1) + 1 / max(v[1][1], 1)) ** 0.5)
        print(f'{100*tgt:>10.2f}%  {d:>+8,} IDs   ({100*d/v[0][0]:+.1f}%)   '
              f'entrapment {v[0][1]}/{v[1][1]}, FDP sigma ~{sig:.2f} pp')
    print('\n  Differences smaller than the FDP sigma are NOT differences.')
