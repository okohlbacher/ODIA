#!/usr/bin/env bash
#
# Build the assay-library test fixtures for the Phase 1 reader.
#
# Source: a DIA-NN 2.0 empirical library written during earlier benchmarking.
# DIA-NN's Parquet and TSV libraries share column names, so one file gives us
# both dialects and no DIA-NN run is needed.
#
# Fixtures land on node-local scratch, not in the repo -- they are derived data
# and the Parquet one is several MB. Re-run on any node that needs them.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${here}/env.sh"

SRC="${DIANN_LIB_PARQUET:-/scratch/kohlbach/bench/s23_diann/report-lib.parquet}"
DEST="${ODIA_SCRATCH}/fixtures"

if [[ ! -f "${SRC}" ]]; then
  echo "source library not found: ${SRC}" >&2
  exit 1
fi

mkdir -p "${DEST}"
cp -n "${SRC}" "${DEST}/diann_library.parquet"

echo "==> writing TSV and a small slice"
"${ODIA_ENV}/bin/python" - "$DEST" <<'PY'
import sys, csv
import pyarrow.parquet as pq

dest = sys.argv[1]
table = pq.read_table(f"{dest}/diann_library.parquet")
cols = table.column_names

def write_tsv(path, tbl):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for batch in tbl.to_batches(max_chunksize=50_000):
            d = batch.to_pydict()
            w.writerows(zip(*(d[c] for c in cols)))

write_tsv(f"{dest}/diann_library.tsv", table)

# A small slice keeps the unit tests fast; take whole precursors, not rows, so
# the transition grouping stays intact.
pid = table.column("Precursor.Id").to_pylist()
keep, seen = set(), []
for p in pid:
    if p not in keep:
        if len(keep) >= 500:
            break
        keep.add(p)
mask = [p in keep for p in pid]
import pyarrow as pa
small = table.filter(pa.array(mask))
write_tsv(f"{dest}/diann_library_small.tsv", small)
pq.write_table(small, f"{dest}/diann_library_small.parquet")

print(f"  full : {table.num_rows:,} rows, {len(set(pid)):,} precursors")
print(f"  small: {small.num_rows:,} rows, {len(keep):,} precursors")
PY

ls -la "${DEST}"
