#!/usr/bin/env python
"""Compare ODIA's predicted fragment spectra against the Python reference.

Covers the inputs the RT model does not have -- charge, NCE and the rank-1
instrument index -- each of which fails silently rather than loudly when wrong:
passing raw NCE gives a spectrum with cosine 0.0028 against the correct one,
and raw charge 0.6377.

Usage: compare_ms2.py <odia_predict_ms2> <model.onnx>
"""
import json
import math
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import peptdeep_reference as ref

# Sequence, charge, NCE, instrument. Charge, NCE and instrument are varied
# deliberately: holding them fixed would let a defect that ignores one of them
# pass, which is how the scale factors were missed in the first place.
CASES = [
    ("ELVISLIVESK", 2, 30.0, "QE"),
    ("ELVISLIVESK", 3, 30.0, "QE"),
    ("PEPTIDEK", 2, 27.0, "Lumos"),
    ("PEPC(Carbamidomethyl)TIDEK", 2, 35.0, "timsTOF"),
    (".(Acetyl)PEPTIDEK", 2, 30.0, "QE"),
    ("PEPTIDER.(Amidated)", 2, 30.0, "SciexTOF"),
    ("AVVPASLSGQDVGSFAYLTIK", 3, 28.0, "NotAnInstrument"),
]
TOLERANCE = 1e-5


def main(tool, model):
    failures = 0
    for seq, charge, nce, instrument in CASES:
        proc = subprocess.run(
            [tool, model, seq, str(charge), str(nce), instrument],
            capture_output=True, text=True,
            env=dict(os.environ, ODIA_ORT_THREADS="1"))
        if proc.returncode != 0:
            print(f"  FAIL {seq} z={charge}: {proc.stderr.strip()}")
            failures += 1
            continue
        cpp = json.loads(proc.stdout)

        os.environ["ODIA_ORT_THREADS"] = "1"
        want = ref.predict_ms2(model, [seq], [charge], nce, instrument)[0]

        if list(cpp["shape"]) != list(want.shape):
            print(f"  FAIL {seq} z={charge}: shape {cpp['shape']} vs {list(want.shape)}")
            failures += 1
            continue

        worst = 0.0
        nonfinite = 0
        for p, row in enumerate(cpp["values"]):
            for c, value in enumerate(row):
                if not math.isfinite(value):
                    nonfinite += 1
                    continue
                worst = max(worst, abs(value - float(want[p][c])))
        if nonfinite:
            print(f"  FAIL {seq} z={charge}: {nonfinite} non-finite intensities")
            failures += 1
        elif worst > TOLERANCE:
            print(f"  FAIL {seq} z={charge} nce={nce} {instrument}: worst {worst:.2e}")
            failures += 1
        else:
            print(f"  ok   {seq} z={charge} nce={nce} {instrument}")

    print(f"MS2 comparison: {failures} failures over {len(CASES)} cases")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
