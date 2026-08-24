#!/usr/bin/env python3
"""Turn a -out_chrom dump into a fixed (precursors, fragments, cycles) tensor.

The dump is long format and contiguous per (precursor, transition), so the whole
thing can be scattered into place vectorised rather than parsed row by row.

Traces are NOT re-centred on anything. The dump's own window -- centred on the
library-predicted retention time by one procedure for every class -- is kept
whole. Re-centring positives on an observed apex is the leak both reviewers
named as the most likely way this experiment produces an impressive and
meaningless number (doc/56 A1).
"""
import sys, numpy as np, pyarrow as pa, pyarrow.csv as pv, pyarrow.parquet as pq

src, out_prefix = sys.argv[1], sys.argv[2]
F = int(sys.argv[3]) if len(sys.argv) > 3 else 12     # fragment slots
C = int(sys.argv[4]) if len(sys.argv) > 4 else 128    # cycle slots

t = pv.read_csv(src,
    read_options=pv.ReadOptions(block_size=1 << 27),
    parse_options=pv.ParseOptions(delimiter='\t'),
    convert_options=pv.ConvertOptions(column_types={
        'Precursor.Id': pa.string(), 'Decoy': pa.int8(),
        'Transition.Index': pa.int32(), 'Product.Mz': pa.float32(),
        'RT': pa.float32(), 'Intensity': pa.float32()}))
print(f'{t.num_rows:,} rows read')

pid = t.column('Precursor.Id').dictionary_encode()
code = np.asarray(pid.combine_chunks().indices, dtype=np.int64)
names = [str(x) for x in pid.combine_chunks().dictionary]
dec  = t.column('Decoy').to_numpy(zero_copy_only=False).astype(np.int8)
tix  = t.column('Transition.Index').to_numpy(zero_copy_only=False).astype(np.int32)
mz   = t.column('Product.Mz').to_numpy(zero_copy_only=False)
rt   = t.column('RT').to_numpy(zero_copy_only=False)
inten= t.column('Intensity').to_numpy(zero_copy_only=False)
del t

# A precursor is (name, decoy): a decoy reconstructs its target's Precursor.Id,
# so keying on the name alone would sum the two into one trace. That exact bug
# cost a 1,000-precursor dump 24 fragments per precursor against 12 in the
# library, and every statistic taken from it.
key = code * 2 + dec
# Runs are contiguous, so a change in (key, transition) starts a new trace.
newrun = np.empty(len(key), dtype=bool)
newrun[0] = True
np.not_equal(key[1:], key[:-1], out=newrun[1:])
newrun[1:] |= (tix[1:] != tix[:-1])
# Transition.Index is the GLOBAL library index (transition_begin[i] + k), not a
# per-precursor slot -- it runs past 32767 on this corpus. Transitions are
# contiguous per precursor, so the slot is the offset from the precursor's own
# first transition.
newprec = np.empty(len(key), dtype=bool)
newprec[0] = True
np.not_equal(key[1:], key[:-1], out=newprec[1:])
prec_start = np.flatnonzero(newprec)
prec_id = np.cumsum(newprec) - 1
tix = tix - tix[prec_start][prec_id]
run_id = np.cumsum(newrun) - 1
run_start = np.flatnonzero(newrun)
cyc = np.arange(len(key), dtype=np.int64) - run_start[run_id]

uk, inv = np.unique(key, return_inverse=True)
P = len(uk)
print(f'{P:,} precursors, {run_start.size:,} traces')

keep = (tix < F) & (cyc < C)
traces = np.zeros((P, F, C), dtype=np.float32)
traces[inv[keep], tix[keep], cyc[keep]] = inten[keep]
mask = np.zeros((P, F), dtype=bool)
mask[inv[keep], tix[keep]] = True

# Per-precursor and per-fragment metadata, taken at each trace's first row.
first = run_start
frag_mz = np.zeros((P, F), dtype=np.float32)
frag_mz[inv[first][tix[first] < F], tix[first][tix[first] < F]] = mz[first][tix[first] < F]
rt0 = np.zeros(P, dtype=np.float32)
rt0[inv[first]] = rt[first]
np.save(f'{out_prefix}_traces.npy', traces)
np.save(f'{out_prefix}_mask.npy', mask)
np.save(f'{out_prefix}_fragmz.npy', frag_mz)

name = [names[k // 2] for k in uk]
pq.write_table(pa.table({'Precursor.Id': pa.array(name),
                         'Decoy': pa.array((uk % 2).astype(np.int8)),
                         'RT0': pa.array(rt0),
                         'N.Fragments': pa.array(mask.sum(1).astype(np.int16))}),
               f'{out_prefix}_meta.parquet')
nz = (traces > 0).mean()
print(f'traces {traces.shape} {traces.nbytes/1e9:.2f} GB   nonzero {100*nz:.1f}%')
print(f'fragments per precursor: median {int(np.median(mask.sum(1)))}  min {mask.sum(1).min()}')
print(f'-> {out_prefix}_{{traces,mask,fragmz,meta}}')
