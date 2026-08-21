#!/usr/bin/env bash
#
# The default ODIA benchmark: one command, ~25 minutes, on the RT-sliced S08
# fixture. Use it to iterate on performance and to try ideas.
#
#   scripts/bench.sh <arm-name> [extra OpenDIAlyzer args...]
#
#   scripts/bench.sh baseline
#   scripts/bench.sh no-floor -min_library_fragments 0
#   scripts/bench.sh wide-im  -precursor_im_window 0.05
#
# Results append to shared/libv2/bench_results.tsv and are printed against the
# stored reference arms. doc/47 has the construction and the evidence.
#
# WHAT THIS CANNOT ANSWER -- printed again at the end of every run, because the
# fixture is 7.5x faster precisely by throwing away 81% of the retention time:
#
#   * FDR / FDP acceptance. Effects attenuate ~5x. The fragment floor is worth
#     6.91 pp on the full run and 1.15 pp here, inside a Poisson sigma of 1.74.
#   * "Where we stand" against DIA-NN. The fixture compresses the ratio from
#     2.96x to 1.40x and moves the two tools' FDP in OPPOSITE directions.
#   * Retention-time calibration. The fixture cannot seed its own map (the fit
#     collapses to a slope of 2.35 against the full run's 1086.50), so the map
#     is supplied and that code path is untested here.
#
# Anything in those three categories needs a full run: shared/libv2/run_full_v5.sh.
set -uo pipefail
[[ $# -ge 1 ]] || { sed -n '2,12p' "$0"; exit 2; }
arm=$1; shift
R=/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer
S=/scratch/kohlbach/odia2x2
L=$R/shared/libv2
FX=/scratch/kohlbach/fixtures/s08_6x60/s08_6x60.mzpeak
[[ -f "$FX" ]] || { echo "fixture missing: $FX (rebuild: shared/libv2/build_fixture.sh)" >&2; exit 1; }
source $R/ODIA/scripts/env.sh

echo "== bench arm '$arm' on s08_6x60 (6 x 60 s, 19.0% of spectra), extra args: $*"
/usr/bin/time -v $R/build-gpu/OpenDIAlyzer \
  -in "$FX" -tr $S/human_v2.parquet \
  -irt_slope 1086.50 -irt_intercept 473.77 -rt_window_pass1 75.4069 \
  -threads 96 -live_memory_gb 400 "$@" \
  -out $S/bench_${arm}.tsv > $L/bench_${arm}.log 2>&1
rc=$?

wall=$(grep -oE "took [0-9:]+ [mh]" $L/bench_${arm}.log | head -1 | awk '{print $2}')
mem=$(grep -oE "Peak Memory Usage: [0-9]+ MB" $L/bench_${arm}.log | grep -oE "[0-9]+")
ids=$(grep -oE "identified [0-9]+ precursors at q <= 0.01" $L/bench_${arm}.log | tail -1 | grep -oE "[0-9]+" | head -1)
echo "   exit $rc  wall $wall  peak ${mem}MB  ids ${ids:-0}"

awk -F'\t' 'FNR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
 ($h["Decoy"]==0||$h["Decoy"]=="false"){k=$h["Precursor.Id"]; s=$h["DScore"]+0;
  if(!(k in b)||s>b[k]){b[k]=s; q[k]=$h["QValue"]}}
 END{for(k in b) print k"\t"b[k]"\t"q[k]}' $S/bench_${arm}.tsv > $L/rank_bench_${arm}.tsv

/ceph/ibmi/abi/oliver/envs/pyprophet/bin/python3 $R/ODIA/scripts/bench_report.py \
  "$arm" "$L/rank_bench_${arm}.tsv" "${wall:-?}" "${mem:-0}"

cat <<'GUARD'

-- guardrails ------------------------------------------------------------
This fixture answers PERFORMANCE and MECHANISM, not FDR and not standing.
  FDP effects attenuate ~5x   (floor: 6.91 pp full -> 1.15 pp here, sigma 1.74)
  DIA-NN ratio compresses     (2.96x full -> 1.40x here)
  RT calibration is untested  (map supplied; the fixture cannot fit one)
For any of those, run shared/libv2/run_full_v5.sh (3h) and compare there.
--------------------------------------------------------------------------
GUARD
