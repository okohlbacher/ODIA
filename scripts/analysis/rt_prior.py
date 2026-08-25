#!/usr/bin/env python3
"""How much does the RT prior already know, and what would using it buy?

Both comparators use predicted retention time at candidate-SELECTION time:
OpenSWATH constrains picking with the normalised-RT window, DIA-NN restricts
its first pass around the predicted RT. ODIA uses it only to CENTRE the
extraction window and then treats every position inside it alike.

This does not propose an RT-delta SCORE -- that was removed deliberately, on
the grounds that a feature's fragments must co-elute and its RT is fixed, so
the quantity can only mark interference. Restricting WHERE a candidate may be
found is a different operation: it is the same thing the extraction window
already does, done less bluntly, and it is the approved plan's
`-rt_window_p95_factor` item.

Measured here: where does the true apex actually sit relative to the window,
and what would tightening the window cost and buy?
"""
import re, numpy as np, pyarrow.parquet as pq
D='/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/corpus'; SP=1.385
lab=pq.read_table(f'{D}/tensor_s08_labels.parquet')
ids=np.array(lab.column('Precursor.Id').to_pylist()); label=np.array(lab.column('Label').to_pylist())
RT0=pq.read_table(f'{D}/tensor_s08_meta.parquet').column('RT0').to_numpy()
AL={'UniMod:4':'Carbamidomethyl'}
nrm=lambda p: re.sub(r'\((.*?)\)',lambda m:'('+AL.get(m.group(1),m.group(1))+')',str(p))
dn=pq.read_table('/scratch/kohlbach/xic/dn_xic.parquet',columns=['Precursor.Id','RT','Q.Value'])
q=dn.column('Q.Value').to_numpy()
TRT={nrm(p):r*60.0 for p,r,k in zip(dn.column('Precursor.Id').to_pylist(),dn.column('RT').to_numpy(),q<=0.01) if k}
sel=[i for i in np.flatnonzero(label=='pos') if ids[i] in TRT]
cy=np.array([ (TRT[ids[i]]-RT0[i])/SP for i in sel ])
inside=cy[(cy>=0)&(cy<128)]
print(f'{len(sel):,} confident positives; {len(inside):,} with the true apex inside the '
      f'128-cycle dump ({100*len(inside)/len(sel):.1f}%)\n')
d=inside-64.0
print(f'  true apex relative to the window CENTRE, in cycles:')
print(f'    median {np.median(d):+.1f}   mean {d.mean():+.1f}   sd {d.std():.1f}')
for p in (5,25,50,75,95):
    print(f'    p{p:<3}{np.percentile(d,p):+8.1f}')
print(f'\n  a window of +-W cycles about the centre would contain:')
for W in (8,16,24,32,48,64):
    keep=100*(np.abs(d)<=W).mean()
    print(f'    +-{W:<3} ({W*SP:5.1f} s)   {keep:5.1f}% of true apexes   '
          f'and {100*(2*W+1)/128:5.1f}% of the positions a candidate could be found at')
print(f'\n  The second column is what the prior BUYS -- every position removed is a '
      f'position\n  that cannot outrank the right one. The first is what it COSTS.')
