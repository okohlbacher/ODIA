#!/usr/bin/env python
"""A reference PeptDeep encoder and RT predictor, written from the spec.

This exists to cross-check ODIA's C++ implementation. It is written from
`doc/04-peptdeep-encoding.md` only -- not from the C++ -- so agreement between
the two is evidence that the spec was implemented, rather than evidence that one
implementation copied the other.

That distinction has mattered repeatedly here: OpenMS's own PeptDeep binding
takes its constants from the same upstream files, so agreeing with it would only
show a shared reading. It is also unusable from outside OpenMS -- the headers are
compiled into libOpenMS.so but never installed.

Usage:
  peptdeep_reference.py encode <modified-sequence>      # dump the tensors
  peptdeep_reference.py rt <model.onnx> <seq> [seq...]  # predict iRT
"""

import json
import os
import sys

import numpy as np

# --- the 109 mod_elements, in yaml order; index IS the feature position -------
# Read from the extracted data file rather than typed: a hand-typed copy came
# out with 112 entries and diverged at index 38.
_ELEMENTS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "..", "data", "peptdeep_mod_elements.txt")
MOD_ELEMENTS = [ln.strip() for ln in open(_ELEMENTS_FILE)
                if ln.strip() and not ln.startswith("#")]
assert len(MOD_ELEMENTS) == 109, f"expected 109 elements, got {len(MOD_ELEMENTS)}"

# Compositions for the modifications the fixtures use, as UniMod gives them.
COMPOSITION = {
    "Acetyl": {"H": 2, "C": 2, "O": 1},
    "Carbamidomethyl": {"H": 3, "C": 2, "N": 1, "O": 1},
    "Oxidation": {"O": 1},
    "GG": {"H": 6, "C": 4, "N": 2, "O": 2},
    "Deamidated": {"H": -1, "N": -1, "O": 1},          # counts are signed
    "Phospho": {"H": 1, "O": 3, "P": 1},
    "1": {"H": 2, "C": 2, "O": 1},
    "4": {"H": 3, "C": 2, "N": 1, "O": 1},
    "35": {"O": 1},
    "21": {"H": 1, "O": 3, "P": 1},
    "121": {"H": 6, "C": 4, "N": 2, "O": 2},
}

CHARGE_SCALE = 0.1
NCE_SCALE = 0.01


def mod_vector(name):
    """Composition -> 109-vector. Unknown elements accumulate in the '?' slot."""
    v = np.zeros(len(MOD_ELEMENTS), dtype=np.float32)
    for element, count in COMPOSITION[name].items():
        if element in MOD_ELEMENTS:
            v[MOD_ELEMENTS.index(element)] = count      # assigned, not added
        else:
            v[-1] += count                              # accumulated
    return v


def parse(seq):
    """(residues, [(site, mod_name)]). Site 0 is the N-terminus, 1..n residues."""
    residues, mods, i = [], [], 0
    while i < len(seq):
        c = seq[i]
        if c in "([":
            close = ")" if c == "(" else "]"
            j = seq.index(close, i)
            name = seq[i + 1:j]
            if name.startswith("UniMod:"):
                name = name.split(":")[1]
            mods.append((len(residues), name))          # attaches to the last residue
            i = j + 1
            continue
        if c.isalpha():
            residues.append(c)
        i += 1
    return residues, mods


def encode(seq):
    """aa_indices [n+2] and mod_x [n+2, 109]."""
    residues, mods = parse(seq)
    n = len(residues)
    aa = np.zeros(n + 2, dtype=np.int64)
    for k, r in enumerate(residues):
        aa[k + 1] = ord(r) - ord("A") + 1                # A->1 .. Z->26, 0 pads
    mod_x = np.zeros((n + 2, len(MOD_ELEMENTS)), dtype=np.float32)
    for site, name in mods:
        mod_x[site] += mod_vector(name)                  # accumulates
    return aa, mod_x


def predict_rt(model_path, sequences):
    import onnxruntime as ort

    session = ort.InferenceSession(model_path, providers=["CPUExecutionProvider"])
    names = [i.name for i in session.get_inputs()]
    out = {}
    # One batch per encoded length: padding is not inert, because index 0 is
    # one-hot encoded and no model applies a padding mask.
    by_length = {}
    for s in sequences:
        by_length.setdefault(len(parse(s)[0]), []).append(s)
    for _, group in sorted(by_length.items()):
        aas, mods = zip(*(encode(s) for s in group))
        feed = {names[0]: np.stack(aas), names[1]: np.stack(mods)}
        values = session.run(None, feed)[0].reshape(-1)
        for s, v in zip(group, values):
            out[s] = float(v)
    return out


if __name__ == "__main__":
    if sys.argv[1] == "encode":
        aa, mod_x = encode(sys.argv[2])
        print(json.dumps({
            "aa_indices": aa.tolist(),
            "mod_x_nonzero": {str(i): mod_x[i][mod_x[i] != 0].tolist()
                              for i in range(mod_x.shape[0]) if mod_x[i].any()},
            "mod_x_nonzero_idx": {str(i): np.nonzero(mod_x[i])[0].tolist()
                                  for i in range(mod_x.shape[0]) if mod_x[i].any()},
        }))
    elif sys.argv[1] == "rt":
        for seq, rt in predict_rt(sys.argv[2], sys.argv[3:]).items():
            print(f"{seq}\t{rt:.6f}")
