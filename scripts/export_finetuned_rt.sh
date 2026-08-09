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

# Record the ONNX ODIA actually consumes. The sidecar written by finetune_rt.py
# hashes rt.pth, which ODIA never opens, so on its own it could not tell you
# whether the model in use is the model that was trained.
onnx="${out}/peptdeep_rt_dynamic.onnx"
if [[ -f "${out}/rt_provenance.json" ]]; then
  "${env_py}" - "${out}/rt_provenance.json" "${onnx}" <<'PYEOF'
import hashlib, json, sys
path, onnx = sys.argv[1], sys.argv[2]
d = json.load(open(path))
d["onnx_model"] = onnx
d["onnx_sha256"] = hashlib.sha256(open(onnx, "rb").read()).hexdigest()
json.dump(d, open(path, "w"), indent=2)
print(f"recorded onnx sha256 {d['onnx_sha256'][:16]}... in {path}")
PYEOF
fi

# ${ODIA_SCRATCH}/build, not ${ODIA_SCRATCH}/build/odia -- the extra path
# component was stale and made this exit 127 AFTER a successful fine-tune and
# export, so the run looked failed when the model was already on disk.
# Non-fatal on purpose: this step is a sanity print, and losing it must not
# discard a model that took real compute to produce.
irt_probe="${ODIA_BUILD:-${ODIA_SCRATCH}/build}/odia_irt_calibration"
if [[ -x "${irt_probe}" ]]; then
  echo "==> ODIA reads it, and refits the iRT line from the standards:"
  "${irt_probe}" "${out}/peptdeep_rt_dynamic.onnx" \
    "${here}/../data/irt_standards.tsv" | head -4
else
  echo "note: no odia_irt_calibration at ${irt_probe}; skipping the read-back check." >&2
fi
