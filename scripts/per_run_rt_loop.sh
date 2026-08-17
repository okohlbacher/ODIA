#!/usr/bin/env bash
#
# Per-run RT fine-tuning loop: pass 1 -> anchors -> fine-tune -> re-predict ->
# full search. The model is built from THIS run and discarded with it.
#
#   per_run_rt_loop.sh <run.mzpeak> <library.parquet> <workdir> [threads]
#
# WHY THIS SHAPE
#
# Both the calibration and the fine-tuned model are PER-RUN objects. Training on
# the run's own confident identifications and predicting for that same run is
# the intended use, not leakage. The failure mode is the opposite one: applying
# a model fine-tuned on one run to ANOTHER cost 2,027 confident precursors
# (library v5 at 35,556 against v4's 37,583). See doc/28 and the -rt_model help.
#
# NO PRE-WIRED SLOPE. -irt_slope/-irt_intercept are deliberately NOT passed.
# -repredict_irt writes raw model output into the iRT column with no standards
# rescale, so that axis carries unlabelled run-fraction units. That is only safe
# because ODIA fits its own map afterwards and the scale cancels under a
# consistently refitted calibration. Freezing an affine on top of an axis whose
# units just changed is precisely how this breaks.
#
# WHAT THIS DOES NOT DO
#
# It never writes fine-tuned predictions into a shipped library. A library is a
# cross-run artefact and must carry stock (OpenMS-bundled) values; the -rt_model
# default is correct as it stands and is not touched here.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/env.sh"

RUN="${1:?usage: per_run_rt_loop.sh <run.mzpeak> <library.parquet> <workdir> [threads]}"
LIB="${2:?library.parquet required}"
W="${3:?workdir required}"
THREADS="${4:-48}"
ODIA="${ODIA_SCRATCH}/build/odia/OpenDIAlyzer"
PY="${ODIA_SCRATCH}/pylibs"
mkdir -p "$W"

# Shared by both passes. Note the absence of -irt_slope/-irt_intercept.
COMMON=(-in "$RUN" -decoys none -entrapment_prefix ENTRAP_
        -rt_window_pass1 60 -rt_window 60
        -live_memory_gb 400 -gate_alpha 0.05 -threads "$THREADS")

echo "### 1/5 pass 1: extract and calibrate, harvest anchors ###"
"$ODIA" "${COMMON[@]}" -tr "$LIB" \
        -stop_after calib -out_anchors "$W/anchors.tsv" \
        -out "$W/pass1.tsv" 2>&1 | tee "$W/pass1.log"
[ -s "$W/anchors.tsv" ] || { echo "NO ANCHORS -- aborting"; exit 1; }
n_anchors=$(( $(wc -l < "$W/anchors.tsv") - 1 ))
echo "    anchors: $n_anchors"
# Measured: 500 anchors improved the RT residual but bought ZERO identifications.
[ "$n_anchors" -ge 2000 ] || { echo "FEWER THAN 2000 ANCHORS ($n_anchors) -- aborting; below this the fine-tune does not pay"; exit 1; }

echo "### 2/5 anchors -> report shape (seconds -> minutes) ###"
PYTHONPATH="$PY" python3 "${here}/anchors_to_report.py" \
        "$W/anchors.tsv" "$W/anchors_report.parquet" 2>&1 | tee "$W/convert.log"

echo "### 3/5 fine-tune on THIS run, protein-level held-out split ###"
# --rt-max-minutes: the RUN's gradient end, so the rt_norm denominator is a
# per-run constant and not a property of whichever peptides were sampled.
RT_MAX=$(PYTHONPATH="$PY" python3 - "$W/anchors_report.parquet" <<'PY'
import sys; sys.path.insert(0,'/scratch/kohlbach/odia/pylibs')
import pyarrow.parquet as pq
print(f"{max(pq.read_table(sys.argv[1], columns=['RT']).column('RT').to_pylist()):.6f}")
PY
)
echo "    rt_norm denominator (run gradient end): ${RT_MAX} min"
PYTHONPATH="$PY" python3 "${here}/finetune_rt.py" \
        --identifications "$W/anchors_report.parquet" \
        --out-dir "$W/model" \
        --rt-max-minutes "$RT_MAX" \
        --method direct --evaluate 2>&1 | tee "$W/finetune.log"

echo "### 4/5 export to ONNX ###"
"${here}/export_finetuned_rt.sh" "$W/model" 2>&1 | tee "$W/export.log"
TUNED="$W/model/peptdeep_rt_dynamic.onnx"
[ -s "$TUNED" ] || { echo "EXPORT PRODUCED NO ONNX -- aborting"; exit 1; }

echo "### 5/5 full search with the per-run model, ODIA fits its own map ###"
"$ODIA" "${COMMON[@]}" -tr "$LIB" \
        -rt_model "$TUNED" -repredict_irt \
        -out "$W/tuned.tsv" 2>&1 | tee "$W/pass2.log"

echo "### validation: identifications, NOT residuals ###"
PYTHONPATH="$PY" python3 "${ODIA_SCRATCH}/massanchor/bench/validate_run.py" \
        "$W/tuned.tsv" "per-run fine-tuned"
echo "PER-RUN RT LOOP DONE"
