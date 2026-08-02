#!/usr/bin/env python
"""Self-checks for the PeptDeep reference encoder.

The reference is the oracle ODIA's C++ encoder will be validated against, so it
needs its own guard. These assertions come from doc/04-peptdeep-encoding.md and
from values obtained by running the shipped models directly.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from peptdeep_reference import MOD_ELEMENTS, encode, mod_vector, parse

failures = []


def check(cond, msg):
    if not cond:
        failures.append(msg)


check(len(MOD_ELEMENTS) == 109, f"element list is {len(MOD_ELEMENTS)}, must be 109")
check(MOD_ELEMENTS[0] == "C", "index 0 must be C")
check(MOD_ELEMENTS[-1] == "?", "last element must be the '?' catch-all")
check(MOD_ELEMENTS[:6] == ["C", "H", "N", "O", "P", "S"],
      "the first six elements are fixed by the yaml's own comment")

# Composition -> vector, with signed counts.
v = mod_vector("Carbamidomethyl")          # H(3)C(2)N(1)O(1)
check(v[MOD_ELEMENTS.index("C")] == 2 and v[MOD_ELEMENTS.index("H")] == 3
      and v[MOD_ELEMENTS.index("N")] == 1 and v[MOD_ELEMENTS.index("O")] == 1,
      "Carbamidomethyl composition is wrong")
check(float(v.sum()) == 7.0, "Carbamidomethyl should contribute exactly 7 atoms")

d = mod_vector("Deamidated")               # H(-1)N(-1)O(1) -- counts are signed
check(d[MOD_ELEMENTS.index("H")] == -1 and d[MOD_ELEMENTS.index("N")] == -1,
      "negative composition counts must survive; ~21% of UniMod rows have them")

# Encoding: terminal tokens, 1-based residue placement, A->1.
aa, mod_x = encode("PEPC(Carbamidomethyl)TIDEK")
check(aa[0] == 0 and aa[-1] == 0, "both terminal tokens must be 0")
# PEPCTIDEK is nine residues, so the encoded length is n+2 = 11.
check(len(aa) == 11 and mod_x.shape == (11, 109),
      f"a 9-mer must encode to n+2=11, got {len(aa)} and {mod_x.shape}")
check(aa[1] == ord("P") - ord("A") + 1, "residue encoding must be A->1")
check(mod_x[4].any() and not mod_x[3].any(),
      "a modification on residue 4 belongs at index 4, not 3")

# An N-terminal modification sits at index 0, not on the first residue.
aa2, mod_x2 = encode("(UniMod:1)PEPTIDEK")
check(mod_x2[0].any(), "an N-terminal modification belongs at index 0")
check(not mod_x2[1].any(), "an N-terminal modification must not land on residue 1")

# Modifications accumulate at a site rather than overwriting.
_, mod_x3 = encode("PEPTIDEK")
check(not mod_x3.any(), "an unmodified peptide must encode an all-zero mod_x")

# The name must never be scanned for residues: "UniMod" contains an M.
residues, _ = parse("(UniMod:1)ADAWEEIRR")
check("".join(residues) == "ADAWEEIRR",
      f"modification names must not contribute residues, got {''.join(residues)}")

for f in failures:
    print(f"  FAIL {f}")
print(f"peptdeep encoding: {len(failures)} failures")
sys.exit(1 if failures else 0)
