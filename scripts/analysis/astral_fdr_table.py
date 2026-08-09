import csv, sys, statistics as st
tset=set(l.strip() for l in open('/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/cmp/astral_truth_ids.txt') if l.strip())
print("truth precursors:", len(tset))
best={}; qof={}; dec=[]
with open(sys.argv[1]) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        d=float(row['DScore']); pid=row['Precursor.Id']
        if row['Decoy']=='1': dec.append(d); continue
        if pid not in best or d>best[pid]: best[pid]=d; qof[pid]=float(row['QValue'])
order=sorted(best.items(), key=lambda kv:-kv[1])
print("target precursors ranked:", len(order))
cum=0; rows=[]
for i,(pid,d) in enumerate(order,1):
    if pid in tset: cum+=1
    rows.append((i,cum,d,qof[pid]))
print("\n  rank  cum_true   emp_FDR  reported_q")
for r_ in (1000,3000,5000,5729,7000,8765,10000):
    if r_<=len(rows):
        i,c,d,q=rows[r_-1]; print(f"{i:6d} {c:9d} {1-c/i:9.4f} {q:11.4f}")
for tgt in (0.01,0.02,0.05):
    bp=bt=0
    for i,c,d,q in rows:
        if 1-c/i<=tgt: bp=i; bt=c
    print(f"empirical FDR <= {tgt:.2f}: longest prefix {bp}, truths {bt}")
# reported q<=0.01 set: what is its true empirical FDR?
sel=[(i,c) for i,c,d,q in rows if q<=0.01]
if sel:
    i,c=sel[-1]
    print(f"\nreported q<=0.01 -> {i} precursors, {c} true, EMPIRICAL FDR {1-c/i:.4f}")
false_d=[d for (pid,d) in order if pid not in tset]
true_d=[d for (pid,d) in order if pid in tset]
q=lambda v,p: sorted(v)[int(p*len(v))-1]
print(f"\ndecoy        p50 {st.median(dec):.3f} p99 {q(dec,.99):.3f} max {max(dec):.3f}  n={len(dec)}")
print(f"false-target p50 {st.median(false_d):.3f} p99 {q(false_d,.99):.3f} max {max(false_d):.3f}  n={len(false_d)}")
print(f"true-target  p50 {st.median(true_d):.3f} p99 {q(true_d,.99):.3f} max {max(true_d):.3f}  n={len(true_d)}")
