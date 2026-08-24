#!/usr/bin/env python3
"""Attach labels to the tensor corpus, with the exclusions doc/56 requires.

    pos          DIA-NN confident, and its apex is INSIDE the window
    unlabelled   DIA-NN confident, apex OUTSIDE the window (A11): a true
                 positive our window does not contain, so it teaches
                 "no peak implies present" if left in the positive class
    ent          entrapment DIA-NN does not call
    ent_hard     entrapment DIA-NN DOES call -- verified false positives, real
                 peptide masses, confidently misidentified. Held out, never
                 trained on: training only on uncalled entrapment measures the
                 EASY absent population and leaves the tail unmeasured (A3)
    unid         ordinary human targets DIA-NN does not call. Absent is
                 UNVERIFIABLE for these, so they are unlabelled too
    decoy        held out entirely -- the check that a presence model ranks
                 them low, not a training class (A2)

Also reports the bias the A11 exclusion induces, because the surviving
positives are those whose iRT prediction is good and that correlates with
sequence properties. Stated as a measurement rather than as a caveat.
"""
import sys, re, numpy as np, pyarrow.parquet as pq, pyarrow as pa

meta_pq, ids_tsv, dn_pq, lib_pq, out_pq = sys.argv[1:6]
SLOPE, INTERCEPT, HALF = 1086.50, 473.77, 90.0

ALIAS = {'UniMod:4': 'Carbamidomethyl'}
def norm(p):
    return re.sub(r'\((.*?)\)', lambda m: '(' + ALIAS.get(m.group(1), m.group(1)) + ')', str(p))

dn = pq.read_table(dn_pq, columns=['Precursor.Id', 'RT', 'Q.Value', 'Precursor.Quantity'])
q = dn.column('Q.Value').to_numpy()
nm = [norm(p) for p in dn.column('Precursor.Id').to_pylist()]
rt = dn.column('RT').to_numpy() * 60.0
qt = dn.column('Precursor.Quantity').to_numpy()
called = {nm[i]: (rt[i], qt[i]) for i in np.where(q <= 0.01)[0]}

grp = {}
with open(ids_tsv) as f:
    f.readline()
    for line in f:
        pid, dcy, g = line.split('\t')[:3]
        grp[(pid, int(dcy))] = g

lib = pq.read_table(lib_pq, columns=['Precursor.Id', 'Decoy', 'RT'])
librt = {(p, int(d)): r for p, d, r in zip(lib.column('Precursor.Id').to_pylist(),
                                           lib.column('Decoy').to_numpy(),
                                           lib.column('RT').to_numpy())}

m = pq.read_table(meta_pq)
ids = m.column('Precursor.Id').to_pylist()
dec = m.column('Decoy').to_numpy()

label, resid = [], []
for pid, d in zip(ids, dec):
    d = int(d)
    if d:
        label.append('decoy'); resid.append(np.nan); continue
    g = grp.get((pid, 0), 'unknown')
    lr = librt.get((pid, 0))
    if g == 'ent':
        label.append('ent_hard' if pid in called else 'ent')
        resid.append(np.nan); continue
    if g == 'unid':
        label.append('unid'); resid.append(np.nan); continue
    if pid in called and lr is not None:
        dt = called[pid][0] - (SLOPE * lr + INTERCEPT)
        resid.append(dt)
        label.append('pos' if abs(dt) <= HALF else 'unlabelled')
    else:
        label.append('unknown'); resid.append(np.nan)

label = np.array(label); resid = np.array(resid, dtype=np.float32)
import collections
for k, v in collections.Counter(label).most_common():
    print(f'  {k:<12} {v:>8,}')

# The A11 exclusion bias, measured rather than asserted.
inw  = np.array([called[p][1] for p, l in zip(ids, label) if l == 'pos'])
outw = np.array([called[p][1] for p, l in zip(ids, label) if l == 'unlabelled'])
if len(outw):
    print(f'\nA11 exclusion bias -- DIA-NN abundance of the positives kept vs dropped:')
    print(f'  kept    n={len(inw):>6,}  median {np.median(inw):>12,.0f}')
    print(f'  dropped n={len(outw):>6,}  median {np.median(outw):>12,.0f}'
          f'   ratio {np.median(outw)/max(np.median(inw),1):.2f}x')

pq.write_table(pa.table({'Precursor.Id': pa.array(ids),
                         'Decoy': pa.array(dec),
                         'Label': pa.array(label),
                         'RT.Residual': pa.array(resid)}), out_pq)
print(f'\n-> {out_pq}')
