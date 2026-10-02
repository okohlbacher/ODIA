#!/usr/bin/env bash
#
# Copy the example DIA runs to node-local scratch for benchmarking.
#
# The canonical copies live on ceph (${ODIA_DATA}) and are read-only. Ceph is a
# shared network filesystem, so timing anything against it measures the network
# as much as the code. /scratch is node-local NVMe; benchmarks must read from
# there, and every node used for benchmarking needs its own copy.
#
# Usage:
#   scripts/stage_data.sh              # stage the small run only (default)
#   scripts/stage_data.sh 12_80 astral # stage named runs
#   scripts/stage_data.sh --all        # stage everything (~24 GB)
#
# Re-running is cheap: files already staged with a matching size are skipped.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${here}/env.sh"

ODIA_DATA="${ODIA_DATA:-${ODIA_ROOT}/data}"
STAGE_DIR="${STAGE_DIR:-${ODIA_SCRATCH}/data}"

ALL_RUNS=(12_80 astral IH1_diaPASEF)

if [[ $# -eq 0 ]]; then
  runs=(12_80)
elif [[ "${1:-}" == "--all" ]]; then
  runs=("${ALL_RUNS[@]}")
else
  runs=("$@")
fi

mkdir -p "${STAGE_DIR}"
echo "==> staging to ${STAGE_DIR} on $(hostname)"

for run in "${runs[@]}"; do
  for ext in mzML mzpeak; do
    srcf="${ODIA_DATA}/${run}.${ext}"
    dstf="${STAGE_DIR}/${run}.${ext}"

    [[ -e "${srcf}" ]] || continue

    if [[ -e "${dstf}" ]] &&
       [[ "$(stat -c%s "${srcf}")" == "$(stat -c%s "${dstf}")" ]]; then
      echo "  skip  ${run}.${ext} (already staged)"
      continue
    fi

    echo "  copy  ${run}.${ext} ($(du -h "${srcf}" | cut -f1))"
    cp "${srcf}" "${dstf}.partial"
    mv "${dstf}.partial" "${dstf}"
  done
done

echo "==> staged:"
du -sh "${STAGE_DIR}"/* 2>/dev/null || echo "  (nothing)"
df -h "${STAGE_DIR}" | tail -1
