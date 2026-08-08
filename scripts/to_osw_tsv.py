#!/usr/bin/env python3
"""Convert ODIA/DIA-NN library TSV to OpenSWATH's transition TSV convention.

OpenSwathWorkflow's TSV reader wants specific column NAMES; ours are DIA-NN's.
It failed with "Could not convert string 'AAAATGTIFTFR2' to a double value",
i.e. it read our Precursor.Id column where it expected a numeric field.
"""
import csv, sys
src, dst = sys.argv[1], sys.argv[2]
COLS = ["PrecursorMz","ProductMz","LibraryIntensity","NormalizedRetentionTime",
        "ProteinId","PeptideSequence","ModifiedPeptideSequence","PrecursorCharge",
        "ProductCharge","FragmentType","FragmentSeriesNumber","transition_group_id",
        "transition_name","Decoy"]
n=0
with open(src, newline='') as fi, open(dst,'w',newline='') as fo:
    r=csv.DictReader(fi, delimiter='\t'); w=csv.writer(fo, delimiter='\t')
    w.writerow(COLS)
    for i,row in enumerate(r):
        pid=row['Precursor.Id']
        seq=row.get('Modified.Sequence', pid)
        # OpenSWATH wants a bare sequence in PeptideSequence; ours may carry mods.
        bare=''.join(c for c in seq if c.isalpha())
        w.writerow([row['Precursor.Mz'], row['Product.Mz'], row['Relative.Intensity'],
                    row.get('RT','0'), row.get('Protein.Group','NA'), bare, seq,
                    row['Precursor.Charge'], row.get('Fragment.Charge','1'),
                    row.get('Fragment.Type','y'), row.get('Fragment.Series.Number','1'),
                    pid, f"{pid}_{i}", row.get('Decoy','0')])
        n+=1
print(f"wrote {n} transitions -> {dst}")
