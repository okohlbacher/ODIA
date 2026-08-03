#!/usr/bin/env bash
#
# Turn a fine-tuned AlphaPeptDeep RT checkpoint into an ONNX file ODIA can use.
#
# The export itself is OpenMS's own script, unmodified -- it already produces
# exactly the input names (input_sequences, mod_x) and output name (rt_pred)
# that ODIA's PeptDeepPredictor expects, so a fine-tuned model reaches the tool
# with no code change on either side.
#
# This is why torch stays out of ODIA's runtime: it is needed here, in an
# offline library-preparation step, and nowhere else.
#
#   export_finetuned_rt.sh <fine-tuned rt.pth> <output dir>
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/env.sh"

pth="${1:?usage: export_finetuned_rt.sh <rt.pth> <outdir>}"
out="${2:?usage: export_finetuned_rt.sh <rt.pth> <outdir>}"
env_py="${ODIA_FINETUNE_PYTHON:-/scratch/kohlbach/odia/rtfinetune/env/bin/python}"
exporter="${ODIA_OPENMS_SOURCE}/tools/scripts/export_peptdeep_models_to_onnx.py"

[[ -x "${env_py}" ]] || { echo "no fine-tuning python at ${env_py}" >&2; exit 1; }
[[ -f "${exporter}" ]] || { echo "no exporter at ${exporter}" >&2; exit 1; }

# The exporter takes a DIRECTORY and insists on all three checkpoints, so the
# stock ms2.pth and ccs.pth are staged beside the fine-tuned rt.pth. Only the
# RT model is replaced; the other two are exported unchanged and are identical
# to the ones already installed, which is why they can simply be overwritten.
stock="${ODIA_PEPTDEEP_MODELS:-/tmp/ptm/generic}"
[[ -f "${stock}/ms2.pth" ]] || {
  echo "no stock checkpoints at ${stock}; unzip pretrained_models_v3.zip there" >&2
  exit 1
}
stage="$(mktemp -d)"
trap 'rm -rf "${stage}"' EXIT
cp "${stock}/ms2.pth" "${stock}/ccs.pth" "${stage}/"
cp "${pth}" "${stage}/rt.pth"

mkdir -p "${out}"
echo "==> exporting $(basename "${pth}") through OpenMS's exporter"
"${env_py}" "${exporter}" --pretrained-dir "${stage}" --out-dir "${out}"

echo "==> ODIA reads it, and refits the iRT line from the standards:"
"${ODIA_BUILD:-${ODIA_SCRATCH}/build/odia}/odia_irt_calibration" \
  "${out}/peptdeep_rt_dynamic.onnx" "${here}/../data/irt_standards.tsv" | head -4
