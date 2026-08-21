#!/usr/bin/env python3
"""Report one bench arm against the stored reference arms.

Entrapment FDP is (e/r)/t over ordinary targets, with r computed on the library
the arm actually searched -- the fragment floor changes it, and using the wrong
one shifts the number by ~0.4%.

Poisson error is printed on every FDP because it is the whole story on this
fixture: ~42 entrapment hits give +-1.2 pp, so differences under ~2.5 pp are
not differences.
"""
import sys, os, datetime
import pyarrow.parquet as pq, pyarrow.compute as pc

L = '/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/libv2'
arm, path, wall, mem = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

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

acc = []
for line in open(path):
    f = line.rstrip('\n').split('\t')
    if f[0] not in ent:
        continue
    try:
        if float(f[2]) <= 0.01:
            acc.append(f[0])
    except (ValueError, IndexError):
        pass

# r on the searched library: if nothing below 3 fragments was accepted the arm
# almost certainly ran with the floor on.
floored = not any(nf[p] < 3 for p in acc)
keep = (lambda p: nf[p] >= 3) if floored else (lambda p: True)
E = sum(1 for p, e in ent.items() if e and keep(p))
T = sum(1 for p, e in ent.items() if not e and keep(p))
r = E / T

e = sum(1 for p in acc if ent[p])
tt = len(acc) - e
fdp = (e / r) / max(tt, 1) * 100
poi = fdp / max(e, 1) ** 0.5
true = tt - e / r

# est.true is carried explicitly rather than recomputed: each arm searched a
# library with its own entrapment ratio, and dividing them all by one constant
# made the same arm print two different numbers.
#
# Keyed by fixture, because an Astral arm compared against the S08 baseline is
# a comparison between two different files. The arm tag carries the fixture as
# a suffix (bench.sh writes <arm>_<fixture>); anything unsuffixed is S08, which
# is what every stored row was measured on.
REFS = {
  's08': [('full v5  (3h00, 124GB)', 13268, 130, 5.72, 0.50, 12387),
          ('fx baseline', 3282, 42, 7.49, 1.16, 2997),
          ('fx no-floor', 2986, 44, 8.64, 1.30, 2688),
          ('fx DIA-NN', 4606, 35, 4.44, 0.75, 4404)],
  # DIA-NN on the FULL library (1,533). The old 1,066 searched a 10,891-precursor
  # pre-selected library and is not comparable. The ODIA baseline row is filled
  # in from the first arm that runs here.
  'astral': [('fx DIA-NN (full lib)', 1533, None, None, None, None)],
}
fixture = 'astral' if arm.endswith('_astral') else 's08'
REF = REFS[fixture]
# What the new arm is judged against: the stored ODIA baseline for this fixture,
# or nothing when there is not one yet.
BASE = {'s08': (7.49, 1.16)}.get(fixture)
print(f"\n{'arm':>24}{'IDs':>8}{'entrap':>8}{'FDP':>9}{'+-':>7}{'est.true':>10}")
for n, i, en, f, p, tr in REF:
    if f is None:
        print(f'{n:>24}{i:>8,}{"-":>8}{"-":>9}{"-":>7}{"-":>10}')
    else:
        print(f'{n:>24}{i:>8,}{en:>8}{f:>8.2f}%{p:>6.2f}%{tr:>10,}')
print(f"{'-'*66}")
print(f'{arm:>24}{len(acc):>8,}{e:>8}{fdp:>8.2f}%{poi:>6.2f}%{true:>10,.0f}'
      f'   wall {wall}  {int(mem)/1000:.0f}GB')

if BASE is None:
    print(f'\n  no stored ODIA baseline for the {fixture} fixture yet -- '
          f'this arm becomes one. Nothing to compare against, and a comparison '
          f'against the other fixture would be a comparison between two files.')
    raise SystemExit(0)
base_fdp, base_sig = BASE
d = fdp - base_fdp
sig = (poi ** 2 + base_sig ** 2) ** 0.5
verdict = 'INSIDE the noise' if abs(d) < sig else f'{abs(d)/sig:.1f} sigma'
base_ids = next(i for n, i, *_ in REF if n == 'fx baseline')
print(f'\n  vs fx baseline ({fixture}): FDP {d:+.2f} pp ({verdict}), '
      f'IDs {len(acc)-base_ids:+,}')

with open(f'{L}/bench_results.tsv', 'a') as fh:
    fh.write(f'{datetime.datetime.now():%Y-%m-%d %H:%M}\t{arm}\t{len(acc)}\t{e}'
             f'\t{fdp:.2f}\t{poi:.2f}\t{true:.0f}\t{wall}\t{mem}\n')
