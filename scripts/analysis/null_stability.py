#!/usr/bin/env python3
"""A column of pure noise moves identifications by up to 11.9% at the operating
point. WHERE does that instability live -- in the discriminant, or in the
q-value threshold placed on it?

The distinction decides what to fix and how expensive it is. If the RANKING is
stable and only the cut moves, the discriminant is fine and the FDR estimator is
the unstable part, which is a far smaller and better-understood object. If the
ranking itself churns, the semi-supervised loop is not converging and the fix is
upstream of everything.

Uses the rank files four arms already produced -- no new run.
"""
import sys, numpy as np
from scipy.stats import spearmanr

L='/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/shared/libv2'
ARMS=['inertcol','nullfeat','null2','null3']
def load(tag):
    d={}
    for line in open(f'{L}/rank_bench_{tag}_ih1.tsv'):
        f=line.rstrip('\n').split('\t')
        if len(f)<3: continue
        try: d[f[0]]=(float(f[1]), float(f[2]))
        except ValueError: continue
    return d
A={t:load(t) for t in ARMS}
for t in ARMS: print(f'{t:>10}: {len(A[t]):,} precursors')
common=set(A[ARMS[0]])
for t in ARMS[1:]: common &= set(A[t])
common=sorted(common)
print(f'\ncommon to all four: {len(common):,}\n')

base=ARMS[0]
b_s=np.array([A[base][k][0] for k in common])
print(f'{"arm":>10}{"rho(DScore)":>14}{"rank churn":>13}{"q<=0.01 base":>14}{"q<=0.01 arm":>13}{"jaccard":>10}')
for t in ARMS[1:]:
    s=np.array([A[t][k][0] for k in common])
    rho=spearmanr(b_s,s).statistic
    rb=np.argsort(np.argsort(-b_s)); ra=np.argsort(np.argsort(-s))
    churn=float(np.mean(np.abs(rb-ra)))/len(common)
    qb=set(k for k in common if A[base][k][1]<=0.01)
    qa=set(k for k in common if A[t][k][1]<=0.01)
    jac=len(qb&qa)/max(len(qb|qa),1)
    print(f'{t:>10}{rho:>14.4f}{churn:>12.2%}{len(qb):>14,}{len(qa):>13,}{jac:>10.3f}')

print('\n  rho near 1 with a moving accepted set => the DISCRIMINANT is stable and')
print('  the THRESHOLD is what moves. rho well below 1 => the semi-supervised loop')
print('  itself is not converging to the same solution.')

# Where in the ranked list does the accepted set diverge?
print('\n  top-N agreement with the baseline ranking:')
print(f'{"arm":>10}' + ''.join(f'{"top-"+str(n):>10}' for n in (500,1000,2000,3000,5000)))
ob=[k for _,k in sorted(zip(-b_s,common))]
for t in ARMS[1:]:
    s=np.array([A[t][k][0] for k in common])
    oa=[k for _,k in sorted(zip(-s,common))]
    row=f'{t:>10}'
    for n in (500,1000,2000,3000,5000):
        row+=f'{len(set(ob[:n])&set(oa[:n]))/n:>10.3f}'
    print(row)
