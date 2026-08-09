import csv, sys
TRUTH='/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/cmp/diann_v6_anchors.tsv'
tset=set()
with open(TRUTH) as f:
    r=csv.DictReader(f, delimiter='\t')
    for row in r: tset.add(row['Precursor.Id'])
print("truth precursors:", len(tset))

SRC=sys.argv[1]
best={}      # pid -> dscore
qof={}
dec=[]
with open(SRC) as f:
    r=csv.DictReader(f, delimiter='\t')
    for row in r:
        d=float(row['DScore']); pid=row['Precursor.Id']
        if row['Decoy']=='1':
            dec.append(d); continue
        if pid not in best or d>best[pid]:
            best[pid]=d; qof[pid]=float(row['QValue'])
print("target precursors ranked:", len(best), " decoy rows:", len(dec))

order=sorted(best.items(), key=lambda kv:-kv[1])
dec.sort(reverse=True)
cum=0; rows=[]
for i,(pid,d) in enumerate(order,1):
    if pid in tset: cum+=1
    rows.append((i,cum,d,qof[pid]))
print("\n  rank  cum_true   emp_FDR  reported_q   dscore")
for r_ in (50,100,200,300,500,738,1000,2000,5000,10000):
    if r_<=len(rows):
        i,c,d,q=rows[r_-1]
        print(f"{i:6d} {c:9d} {1-c/i:9.3f} {q:11.4f} {d:8.3f}")
for tgt in (0.01,0.05,0.10,0.20,0.50,0.70):
    bestp=0; bestt=0
    for i,c,d,q in rows:
        if 1-c/i<=tgt: bestp=i; bestt=c
    print(f"empirical FDR <= {tgt:.2f}: longest prefix {bestp}, truths {bestt}")
# Where do decoys sit relative to FALSE targets?
false_d=[d for (pid,d) in order if pid not in tset]
import statistics as st
print(f"\ndecoy dscore   p50 {st.median(dec):.3f}  p99 {sorted(dec)[int(.99*len(dec))]:.3f}  max {max(dec):.3f}")
print(f"false-target   p50 {st.median(false_d):.3f}  p99 {sorted(false_d)[int(.99*len(false_d))]:.3f}  max {max(false_d):.3f}")
true_d=[d for (pid,d) in order if pid in tset]
print(f"true-target    p50 {st.median(true_d):.3f}  p99 {sorted(true_d)[int(.99*len(true_d))]:.3f}  max {max(true_d):.3f}")
