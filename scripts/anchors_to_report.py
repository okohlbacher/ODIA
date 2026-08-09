#!/usr/bin/env python
"""Turn -out_anchors TSV into the report.parquet shape finetune_rt.py reads.

The fine-tuner was written against a DIA-NN report and wants Modified.Sequence,
Precursor.Charge, RT and Q.Value, with RT in MINUTES. Our anchors carry the
same facts under different names and in run SECONDS, so this is a rename and a
divide -- but the divide is the whole point: feeding seconds silently trains
the model on a 60x stretched gradient.
"""
import csv, re, sys
import pyarrow as pa, pyarrow.parquet as pq

# Named modifications -> UniMod, because finetune_rt.py accepts only UniMod and
# REFUSES to guess at names rather than corrupt the training data silently.
# Extend this table explicitly; an unmapped name drops the peptide, and the
# count is printed so a silent loss is impossible.
UNIMOD = {
    "Carbamidomethyl": 4,
    "Oxidation": 35,
    "Acetyl": 1,
    "Phospho": 21,
    "GG": 121,
}
NAMED = re.compile(r"\(([A-Za-z][A-Za-z0-9 -]*)\)")

def to_unimod(seq):
    """Return the sequence with named mods rewritten, or None if any is unknown."""
    unknown = []
    def sub(m):
        name = m.group(1)
        if name.startswith("UniMod:"):
            return m.group(0)
        if name in UNIMOD:
            return f"(UniMod:{UNIMOD[name]})"
        unknown.append(name)
        return m.group(0)
    out = NAMED.sub(sub, seq)
    return None if unknown else out
src, dst = sys.argv[1], sys.argv[2]
seq, ch, rt, q = [], [], [], []
with open(src) as f:
    dropped = {}
    for r in csv.DictReader(f, delimiter='\t'):
        m = to_unimod(r['Modified.Sequence'])
        if m is None:
            for n in NAMED.findall(r['Modified.Sequence']):
                if n not in UNIMOD and not n.startswith("UniMod:"):
                    dropped[n] = dropped.get(n, 0) + 1
            continue
        seq.append(m); ch.append(int(r['Precursor.Charge']))
        rt.append(float(r['Observed.RT']) / 60.0); q.append(float(r['QValue']))
if dropped:
    print("DROPPED unmapped modifications:", dict(sorted(dropped.items(),
          key=lambda kv: -kv[1])), file=sys.stderr)
pq.write_table(pa.table({"Modified.Sequence": seq, "Precursor.Charge": ch,
                         "RT": rt, "Q.Value": q}), dst)
print(f"{src} -> {dst}: {len(seq)} rows, RT {min(rt):.2f}-{max(rt):.2f} min")
