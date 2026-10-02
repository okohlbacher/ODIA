#!/usr/bin/env bash
#
# Fine-tune the retention-time model on a run's own identifications, and export
# it as an ONNX file ODIA reads with -rt_model.
#
#   finetune_rt.sh <report.parquet for THIS run> <output dir> [max peptides]
#
# What this is for, and what it is not: the stock model's output saturates at
# the end of the gradient, and no calibration can undo that. Fine-tuning does.
# It buys a narrower extraction window, not identifications -- a perfect RT
# column was measured to be worth -108 precursors on IH1. See
# doc/06-rt-refinement-plan.md before deciding it is worth running.
#
# The model produced is specific to the run it was tuned on. Reusing it across
# runs or gradients is wrong by construction; the sidecar says so.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/env.sh"

ids="${1:?usage: finetune_rt.sh <report.parquet> <outdir> [max peptides]}"
out="${2:?usage: finetune_rt.sh <report.parquet> <outdir> [max peptides]}"
n="${3:-2000}"
env_py="${ODIA_FINETUNE_PYTHON:-/scratch/kohlbach/odia/rtfinetune/env/bin/python}"

[[ -x "${env_py}" ]] || {
  echo "no python with peptdeep and torch at ${env_py}." >&2
  echo "Create one -- NOT in ${ODIA_ENV}, which is pinned and shared:" >&2
  echo "  micromamba create -p <prefix> python=3.11 && <prefix>/bin/pip install peptdeep onnx" >&2
  exit 1
}

# --evaluate is ON by default. The docstring claims every run reports whether
# fine-tuning helped on its own data; without this flag the shipped path
# reported nothing, and the one defensive feature that stops a user inheriting
# our verdict from a different instrument was off for everybody.
echo "==> fine-tuning on ${ids}"
# Neither --epochs nor --q-value was passed, so both looked like knobs and were
# silently ignored. ODIA_RT_EPOCHS=80 ran the default 40 and duplicated its own
# control; and an arm fed a q<=0.05 selection of 7,173 precursors was filtered
# back to 5,867 by finetune_rt.py's own --q-value default, so it duplicated the
# q<=0.01 arm and reported it as a data-size result. Both plumbed through.
"${env_py}" "${here}/finetune_rt.py" "${ids}" "${out}" --max-peptides "${n}" \
  --epochs "${ODIA_RT_EPOCHS:-40}" \
  --q-value "${ODIA_RT_QVALUE:-0.01}" \
  ${ODIA_RT_CURVE:+--curve "${ODIA_RT_CURVE}"} \
  --device "${ODIA_RT_DEVICE:-gpu}" \
  ${ODIA_RT_THREADS:+} \
  --method "${ODIA_RT_METHOD:-direct}" --evaluate

echo "==> exporting to ONNX"
"${here}/export_finetuned_rt.sh" "${out}/rt.pth" "${out}"

echo
echo "Use it with:  OpenDIAlyzer -fasta <fasta> -rt_model ${out}/peptdeep_rt_dynamic.onnx"
echo "Provenance:   ${out}/rt_provenance.json"
