#!/usr/bin/env python
"""Compare ODIA's predicted iRT against the independent Python reference.

Tensor agreement is not prediction agreement: dtype, layout, stride or ordering
assumptions can differ while the encoded values match. This compares what the
models actually return.

The sequence list deliberately mixes lengths, so a broken length-grouping or a
lost reordering shows up as wrong values rather than as an error.

Usage: compare_predictions.py <odia_predict_rt> <model.onnx>
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import peptdeep_reference as ref

SEQUENCES = [
    "ELVISLIVESK",                       # 11
    "PEPTIDEK",                          # 8  -- different length, tests grouping
    "PEPTIDER",                          # 8
    "PEPC(Carbamidomethyl)TIDEK",        # 9, residue modification
    ".(Acetyl)PEPTIDEK",                 # 8, N-terminal
    "PEPTIDER.(Amidated)",               # 8, C-terminal
    "AVVPASLSGQDVGSFAYLTIK",             # 21
    "PEPTIDEM(Oxidation)K",              # 9
]
TOLERANCE = 1e-5


def main(tool, model):
    proc = subprocess.run([tool, model, *SEQUENCES], capture_output=True, text=True)
    if proc.returncode != 0:
        print(f"  FAIL ODIA predictor failed: {proc.stderr.strip()}")
        return 1

    cpp = {}
    for line in proc.stdout.strip().split("\n"):
        seq, value = line.rsplit("\t", 1)
        cpp[seq] = float(value)

    python = dict(zip(SEQUENCES, ref.predict_rt(model, SEQUENCES)))

    failures = 0
    for seq in SEQUENCES:
        if seq not in cpp:
            print(f"  FAIL {seq}: no prediction from ODIA")
            failures += 1
            continue
        d = abs(cpp[seq] - python[seq])
        if d > TOLERANCE:
            print(f"  FAIL {seq}: ODIA {cpp[seq]:.6f} vs reference {python[seq]:.6f} "
                  f"(delta {d:.2e})")
            failures += 1
        else:
            print(f"  ok   {seq}: {cpp[seq]:.6f}")

    print(f"prediction comparison: {failures} failures over {len(SEQUENCES)} sequences")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
