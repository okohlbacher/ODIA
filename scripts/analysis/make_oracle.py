#!/usr/bin/env python3
"""Write an -oracle_rt file from a DIA-NN report: Precursor.Id and run seconds.

The alias is applied here, once, so the tool never has to know about it. If the
run log then reports a low match rate, that is a real coverage gap and not a
spelling mismatch.
"""
import sys, re, numpy as np, pyarrow.parquet as pq

dn_pq, out = sys.argv[1], sys.argv[2]
ALIAS = {'UniMod:4': 'Carbamidomethyl'}
def norm(p):
    return re.sub(r'\((.*?)\)', lambda m: '(' + ALIAS.get(m.group(1), m.group(1)) + ')', str(p))

t = pq.read_table(dn_pq, columns=['Precursor.Id', 'RT', 'Q.Value'])
q = t.column('Q.Value').to_numpy()
ids = [norm(p) for p in t.column('Precursor.Id').to_pylist()]
rt = t.column('RT').to_numpy() * 60.0        # DIA-NN reports minutes
keep = np.where(q <= 0.01)[0]
with open(out, 'w') as f:
    f.write('Precursor.Id\tRT\n')
    for i in keep:
        f.write(f'{ids[i]}\t{rt[i]:.3f}\n')
print(f'{len(keep)} precursors -> {out}')
