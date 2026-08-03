#!/usr/bin/env python
"""Compare ODIA's C++ PeptDeep encoder against the independent Python reference.

The two were written from doc/04-peptdeep-encoding.md separately -- the C++ from
OpenMS's AASequence, the Python from the raw string -- so agreement is evidence
that the spec was implemented rather than that one copied the other.

Usage: compare_encoders.py <odia_encode_dump> [sequence...]
"""
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import peptdeep_reference as ref

DEFAULT = [
    "PEPTIDEK",                          # unmodified
    "PEPC(Carbamidomethyl)TIDEK",        # residue modification
    ".(Acetyl)PEPTIDEK",                 # N-terminal, row 0
    "PEPTIDER.(Amidated)",               # C-terminal, row nAA+1, negative count
    "PEPTIDEM(Oxidation)K",              # a different residue modification
    ".(Acetyl)PEPC(Carbamidomethyl)TIDEM(Oxidation)K",   # several at once
    "AVVPASLSGQDVGSFAYLTIK",             # long, unmodified
    "C(Carbamidomethyl)C(Carbamidomethyl)PEPTIDEK",      # adjacent modifications
]


def main(dump, sequences):
    failures = 0
    for seq in sequences:
        proc = subprocess.run([dump, seq], capture_output=True, text=True)
        if proc.returncode != 0:
            print(f"  FAIL {seq}: C++ encoder failed: {proc.stderr.strip()}")
            failures += 1
            continue
        cpp = json.loads(proc.stdout)

        aa, mod_x = ref.encode(seq)
        if cpp["aa_indices"] != aa.tolist():
            print(f"  FAIL {seq}: aa_indices differ\n"
                  f"        C++    {cpp['aa_indices']}\n        python {aa.tolist()}")
            failures += 1
            continue

        py_rows = {str(r): {str(c): float(mod_x[r][c])
                            for c in np.nonzero(mod_x[r])[0]}
                   for r in range(mod_x.shape[0]) if mod_x[r].any()}
        cpp_rows = {r: {c: float(v) for c, v in cols.items()}
                    for r, cols in cpp["mod_x"].items()}
        if cpp_rows != py_rows:
            print(f"  FAIL {seq}: mod_x differs\n"
                  f"        C++    {cpp_rows}\n        python {py_rows}")
            failures += 1
            continue
        print(f"  ok   {seq}")

    print(f"encoder comparison: {failures} failures over {len(sequences)} sequences")
    return 1 if failures else 0


if __name__ == "__main__":
    seqs = sys.argv[2:] or DEFAULT
    sys.exit(main(sys.argv[1], seqs))
