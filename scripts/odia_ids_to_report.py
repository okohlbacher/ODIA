#!/usr/bin/env python
"""ODIA's own confident identifications, in the shape the RT fine-tuner reads.

Round 2 of the bootstrap. Pass 1's anchors are a few thousand; a completed
search identifies several times more, and doc/06 measured the fine-tune's gain
scaling with the training set. So the second round trains on what the first
round's search found.

The circularity is real and is the point: this is a per-run refinement learning
THIS chromatography, exactly as doc/06 section R-A specifies. What it must not
become is a model persisted and reused on another run -- that cost 2,027
precursors when it happened.

Only targets at q <= the threshold, best row per precursor, RT in MINUTES.
"""
import csv, re, sys
import pyarrow as pa, pyarrow.parquet as pq

src, dst = sys.argv[1], sys.argv[2]
qmax = float(sys.argv[3]) if len(sys.argv) > 3 else 0.01
UNIMOD = {"Carbamidomethyl": 4, "Oxidation": 35, "Acetyl": 1, "Phospho": 21, "GG": 121}
NAMED = re.compile(r"\(([A-Za-z][A-Za-z0-9 -]*)\)")

def to_unimod(seq):
    bad = []
    def sub(m):
        n = m.group(1)
        if n.startswith("UniMod:"): return m.group(0)
        if n in UNIMOD: return f"(UniMod:{UNIMOD[n]})"
        bad.append(n); return m.group(0)
    out = NAMED.sub(sub, seq)
    return None if bad else out

best = {}
csv.field_size_limit(10**9)
with open(src) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r['Decoy'] != '0': continue
        q = float(r['QValue'])
        if q > qmax: continue
        pid = r['Precursor.Id']
        d = float(r['DScore'])
        if pid not in best or d > best[pid][0]:
            best[pid] = (d, float(r['RT']), q)

seq, ch, rt, qs, dropped = [], [], [], [], {}
# Two formats, because the tool writes two. ChromatogramTsv emits
# "<modified sequence><charge>" with no separator (AAAATGTIFTFR2) while
# -out_anchors emits "<modified sequence>_<charge>". Parsing only the
# underscore form silently produced ZERO rows from a 45,638-row file.
ID = re.compile(r"^(.*[A-Za-z\)])_?(\d+)$")
for pid, (d, t, q) in best.items():
    hit = ID.match(pid)
    if not hit: continue
    m, c = hit.group(1), hit.group(2)
    u = to_unimod(m)
    if u is None:
        for n in NAMED.findall(m):
            if n not in UNIMOD and not n.startswith("UniMod:"):
                dropped[n] = dropped.get(n, 0) + 1
        continue
    try: charge = int(c)
    except ValueError: continue
    seq.append(u); ch.append(charge); rt.append(t / 60.0); qs.append(q)

pq.write_table(pa.table({"Modified.Sequence": seq, "Precursor.Charge": ch,
                         "RT": rt, "Q.Value": qs}), dst)
print(f"{src} -> {dst}: {len(seq)} precursors at q<={qmax}, "
      f"RT {min(rt):.2f}-{max(rt):.2f} min")
if dropped: print("DROPPED unmapped modifications:", dropped, file=sys.stderr)
