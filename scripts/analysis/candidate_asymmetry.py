import csv, sys, collections, statistics as st
cnt=collections.Counter(); best={}
with open(sys.argv[1]) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        k=(row['Precursor.Id'], row['Decoy']); cnt[k]+=1
        d=float(row['DScore'])
        if k not in best or d>best[k]: best[k]=d
for dec in ('0','1'):
    ns=[v for (p,dd),v in cnt.items() if dd==dec]
    bs=sorted(v for (p,dd),v in best.items() if dd==dec)
    lab='decoy' if dec=='1' else 'target'
    print(f"{lab:7s} precursors {len(ns):6d}  cands/prec mean {sum(ns)/len(ns):5.2f} median {st.median(ns):4.0f}  "
          f"best-dscore p50 {bs[len(bs)//2]:6.3f} p99 {bs[int(.99*len(bs))]:6.3f} max {bs[-1]:6.3f}")
