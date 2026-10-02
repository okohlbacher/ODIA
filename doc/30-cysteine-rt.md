# The cysteine retention-time problem: whose fault is it?

Answering three questions asked together: is the cysteine issue resolved, is it
systematic in AlphaPeptDeep or is ODIA feeding the model the modification
wrongly, and would fixing it recover identifications when DIA-NN searches our
library.

Short answers: **the alkylation half is resolved and was ours; the retention-time
half is not ours and not fixable inside the model — it is shared, in the same
direction and to within 0.03 min, by DIA-NN's completely independent predictor;
and fixing it would not recover the DIA-NN identification gap, because the gap is
not in cysteine.**

## Two different things were both called "the CYS issue"

1. **Alkylation state of the search space** — whether the library carries
   Carbamidomethyl on C. Getting this wrong makes every cysteine precursor mass
   wrong, so those peptides simply cannot be found. This was ODIA's defect
   (doc/27), it is fixed, and it is now a written-down config key rather than a
   hard-coded constant.
2. **Retention time of cysteine peptides** — whether the RT model predicts them
   correctly once the masses are right. This is the subject of this note. It is
   not ODIA's defect.

## Finding 0: the two benchmarks differ in alkylation state

Same FASTA, same DIA-NN settings, same CAM-free search:

| run | identified precursors | containing C | rate |
|---|---:|---:|---:|
| Astral | 12,308 | 2,268 | **18.43%** |
| IH1 (diaPASEF) | 33,330 | 646 | **1.94%** |

DIA-NN's own Astral output library is 18.49% cysteine-containing, so ~18% is
what an unalkylated sample looks like. IH1 is depleted ~9-fold.

**IH1 is alkylated and the frozen IH1 reference searched it CAM-free**, losing
almost all of its cysteine peptides. Astral is not alkylated. The reference
number 33,330 is therefore a handicapped number, not a ceiling.

Both reviewers rejected the premise as stated, and correctly: a search engine may
omit a *fixed* modification from the modified-sequence string, so "no
modifications in the report" proves nothing. The reported precursor m/z does.
Theoretical bare-cysteine m/z minus reported m/z, over the confident sets:

| | 0 Cys | 1 Cys | 2+ Cys |
|---|---:|---:|---:|
| Astral, assuming bare C | −0.000008 | −0.000006 | −0.000004 Th |
| Astral, assuming CAM | −0.000008 | **+28.5107** | **+57.0214 Th** |
| IH1, assuming bare C | −0.000012 | −0.000010 | +0.000001 Th |
| IH1, assuming CAM | −0.000012 | **+28.5107** | **+38.0143 Th** |

Both reference searches used bare-cysteine masses, to within 10 µTh. So the
cysteine peptides identified in *both* samples carry free thiols, and the IH1
depletion is a property of the sample, not of the mass convention.

What that still does not settle is *what* is on the IH1 cysteines. Alkylation is
the obvious candidate; kimi's alternative — an unreduced sample, where disulfides
remove the same peptides with no alkylation at all — predicts the same 9-fold
depletion. The separating run is DIA-NN on IH1 with `--unimod4`: cysteine returns
to ~15%+ only if carbamidomethyl is what is on the thiol. A rise to 3–4% would
mean something else is.

## Finding A: the modification does reach the model; the model does not use it

`odia_predict_rt` runs the same encoder and the same ONNX session as the library
generator. On one backbone, length held constant (rt_norm units):

```
SAMPLCEKTIDER                     0.341713
SAMPLC(Carbamidomethyl)EKTIDER    0.315007
SAM(Oxidation)PLCEKTIDER          0.219900
S(Acetyl)AMPLCEKTIDER             0.434649
```

At scale, toggling one cysteine at a time (600 single-Cys and 300 two-Cys
peptides from the Astral identifications):

```
single CAM on a single-Cys peptide   median +0.00200  IQR [-0.00934,+0.01535]  44.8% negative
two-Cys, site 1 / site 2 / both      +0.00287 / -0.00148 / -0.00042
additivity  median(both - (s1+s2))   -0.00068     corr(both, s1+s2) 0.968
sign agreement between the two sites  58.7%
```

The response is **additive** — so nothing is being masked by toggling every
cysteine at once, which was kimi's first alternative — but its **direction is a
coin flip**.

The encoding itself is not the problem. Both reviewers asked for the tensor, not
more inference, so here it is, for `PEPC(Carbamidomethyl)TIDEK`:

```
{"rows": 1, "sequence_length": 11,
 "aa_indices": [0, 16, 5, 16, 3, 20, 9, 4, 5, 11, 0],
 "mod_x": {"4": {"0": 2.0, "1": 3.0, "2": 1.0, "3": 1.0}}}
```

Slots 0–3 are C, H, N, O: exactly C2H3NO, nothing in any other slot, on row 4 —
the cysteine at index 3 plus one for the N-terminal token. The element list is
extracted programmatically from AlphaPeptDeep's own `model_const.yaml`
(`scripts/generate_peptdeep_elements.py`, which refuses a list that is not
exactly 109 long) rather than hand-typed, the assignment/accumulation convention
matches upstream, and `test/compare_encoders.py` diffs `mod_x` element-wise
against an independently written Python reference on eight sequences including
CAM, adjacent CAMs, and CAM with terminal modifications.

What this is NOT evidence for: kimi is right that under the mechanism below the
toggle *should* be small and unstable, because bare cysteine is off-manifold for
the model. The toggle result is therefore consistent with the mechanism, not
proof of it, and the "wrong sign and 18× too small" reading of it is withdrawn —
it presumes the toggle is a valid physical probe, which the mechanism denies.
The proof is Finding B.

## Finding B: the bias is real, and DIA-NN has exactly the same one

Astral is unalkylated, so bare cysteine is the correct prediction there.
Calibrating each predictor against observed RT on **non-cysteine peptides only**
and reading the cysteine residual (residual = observed − predicted):

| predictor | 1 Cys | 2 Cys | 3+ Cys |
|---|---:|---:|---:|
| ODIA / AlphaPeptDeep, bare C | +0.896 | +1.558 | +1.607 min |
| ODIA / AlphaPeptDeep, CAM-C | +0.840 | +1.647 | +1.591 min |
| **DIA-NN's own predictor** | **+0.902** | **+1.569** | **+2.214 min** |

Non-cysteine median residual is −0.04 min (ODIA) and −0.01 min (DIA-NN); the
non-cysteine residual sd is 1.10 and 1.14 min. So cysteine peptides elute
**~0.9 min later per cysteine than either model predicts**, and the two agree to
within 0.01 min at one cysteine.

**Replicated on the other instrument and gradient.** IH1's 646 cysteine
identifications are peptides that escaped alkylation, i.e. genuinely free thiols,
on a timsTOF with a 7–30 min gradient:

| | 1 Cys | 2+ Cys |
|---|---:|---:|
| ODIA / AlphaPeptDeep | +0.739 min = **+0.0362 rt_norm** | +1.484 min = +0.0727 |
| DIA-NN's own predictor | +1.034 min | +2.072 min |

Astral in the same units: **+0.0323** and +0.0563 rt_norm. Two instruments, two
gradients, two samples, one number.

**And DIA-NN silently corrects it after the fact.** Its output library is
empirically refined against observed RT. Output-library RT minus predicted
library iRT, calibrated on non-cysteine:

```
 0 Cys n=9964  +0.019 iRT = +0.003 min
 1 Cys n=1780  +5.030 iRT = +0.880 min
 2 Cys n= 413  +8.992 iRT = +1.574 min
 3+Cys n=  68  +11.299 iRT = +1.977 min
```

DIA-NN's own refinement moves cysteine peptides by the same +0.88/+1.57/+1.98
min we measure as its prediction error. Three independent confirmations of one
number.

**This also retracts an earlier finding.** The "our library diverges from
DIA-NN's by −5.8/−10.6/−9.5 iRT on cysteine" measurement compared our
*prediction* against DIA-NN's *empirically refined* output library. At 0.175
min/iRT that is −1.02/−1.86/−1.66 min — the same shared prediction error,
correctly measured and wrongly attributed to a difference between the two
predictors. Recomputed against DIA-NN's *predicted* library, the two agree on
cysteine to within 0.04–0.15 min.

### The mechanism, and the alternative that was tested and refuted

Both models were trained on corpora in which essentially every cysteine carried
carbamidomethyl. The modification is therefore perfectly confounded with the
residue: "C" in the model *means* CAM-C, and free cysteine — which is more
hydrophobic, hence later-eluting — cannot be represented. This is
constant-modification unidentifiability, and "correlated training corpora" is
not a competing explanation but the same one restated.

The competing explanation that *was* live is kimi's: calibration curvature.
Cysteine peptides are more hydrophobic, so they cluster late in the gradient; a
straight-line calibration fitted on non-cysteine peptides would hand any
late-eluting class a one-signed residual, in both predictors, with no chemistry
involved. The saturation pattern (0.90 → 1.56 → 1.61 rather than linear in
count) fits that story too.

It is refuted three ways:

```
                                       Astral 1C   Astral 2C   IH1 1C
linear calibration on non-Cys            +0.896      +1.558     +0.739 min
local-linear (nonparametric) on non-Cys  +0.927      +1.581     +0.860 min
RT-matched to non-Cys at equal pred. RT  +0.956      +1.588     +0.905 min
```

The bias does not shrink; it grows slightly. And the premise fails outright on
Astral, where cysteine peptides elute *earlier* than average, not later
(predicted-RT median 0.381 for 1 Cys against 0.439 for non-Cys) — while on IH1
they elute later (0.659 against 0.586). Opposite positions in the gradient, same
bias, so it cannot be a gradient-position artefact.

There is genuine curvature: nonparametric calibration cuts the non-cysteine
residual sd from 1.101 to 0.877 min on Astral and from 1.029 to 0.685 min on IH1.
That is a separate, real improvement worth taking, and it is not the cysteine
effect.

### Corrections the reviewers forced

**One earlier claim here was simply wrong.** I wrote that DIA-NN's selection
biases DIA-NN's residuals downward "and ours not at all". Codex is right that
this is false: conditioning on DIA-NN identification selects small DIA-NN error,
and because the two predictors' errors are correlated, it selects small ODIA
error too. Both series are biased downward by the same conditioning. Measured, the
two Cys residual distributions have nearly the same shape (ODIA IQR 1.40 sd, p95
+2.89; DIA-NN 1.27 and +2.70), so the selection effect is small — but the
asymmetry argument is withdrawn, not weakened.

**DIA-NN is not an independent replication in the strong sense** either: both
predictors were trained on similar alkylated corpora, which is the mechanism
rather than a competing explanation. What DIA-NN's numbers do establish is that
the bias is not an ODIA implementation defect, and the empirical-refinement table
above is the stronger evidence, because it is DIA-NN correcting DIA-NN.

**The effect is not linear in cysteine count**, and the note above never should
have implied it: 0.032 / 0.056 / 0.058 rt_norm for 1 / 2 / 3+ is well short of
n × 0.032. The correction is bucketed for exactly this reason.

### Registered in advance: what the `--unimod4` runs must show

Written before the runs finished, so the result cannot be reinterpreted after
the fact. Adapted from codex's criteria, in the predictor's own units:

| outcome | conclusion |
|---|---|
| cysteine rescued to ~15%+ at +57.021/z, and CAM-cysteine residual within ±0.008 rt_norm of zero | IH1 is alkylated; free-vs-CAM confirmed as the dominant mechanism |
| cysteine rescued, but a positive residual near the Astral magnitude remains | IH1 is alkylated; the free-vs-CAM explanation is **falsified** — something else about cysteine is mis-modelled |
| cysteine rises only to 3–4% | not carbamidomethyl; kimi's unreduced-disulfide or another reagent |
| CAM residual comes out comparably **negative** | the model is not simply "already predicting CAM"; Finding B fails as stated |

Both reviewers also note what would settle it outright and what we do not have:
the same peptide backbone measured free and carbamidomethylated under one
chromatography. Astral and IH1 differ in tissue, instrument and gradient, so
even the favourable outcome is supportive rather than decisive.

### Outcome: a reciprocal control, and it is unambiguous

Each sample was searched both ways. The two results are mirror images:

| | precursors | cysteine-containing | protein groups |
|---|---:|---:|---:|
| **IH1** CAM-free (the frozen reference) | 33,330 | 646 (1.94%) | 5,500 |
| **IH1** `--unimod4` | **37,334** | **3,454 (9.25%)** | 5,660 |
| **Astral** CAM-free (the frozen reference) | 12,308 | 2,268 (18.43%) | 1,328 |
| **Astral** `--unimod4` | **9,837** | **190 (1.93%)** | 1,219 |

Forcing the wrong alkylation state on either sample collapses its cysteine
identifications to ~1.9%, and the right one restores them. IH1 is alkylated;
Astral is not. Kimi's unreduced-disulfide alternative is excluded — disulfides
would not be rescued by declaring carbamidomethyl.

IH1's rescue reaches 9.25%, not Astral's 18%. That is a 5.3-fold recovery and it
brings 4,004 extra precursors and 160 extra protein groups, so alkylation is
plainly the dominant effect; the shortfall against 18% is unexplained and could
be incomplete alkylation, a different tissue's cysteine content, or a second
cysteine species. Codex's stricter bar ("15–20% compelling, 3–4% not enough")
is met in the middle, so the conclusion is stated at that strength.

**The frozen IH1 reference should be re-cut.** 33,330 was measured with the
wrong alkylation assumption; 37,334 is the honest number for this sample.

### The Astral file's provenance, and a conflict with its paper

`data/astral.mzML` carries no filename — MassIVE's converter renames the source
to `source_file.raw` — but the conversion task id embedded in its `sourceFile`
location resolves the dataset:

* **MassIVE MSV000100943 / ProteomeXchange PXD060573**
* Beimers, Overmyer, Sinitcyn, Lancaster, Quarmby & Coon, *Technical Evaluation
  of Plasma Proteomics Technologies*, J. Proteome Res. (2025),
  doi:10.1021/acs.jproteome.5c00221; preprint bioRxiv 2025.01.08.632035
* Acquired 2024-09-15, Orbitrap Astral serial OA10144, MS1 380–980, 150 windows
  of 4 Th, HCD 25, MS2 150–2000, 2333 s gradient

**The paper says the plasma was alkylated.** Both the neat and the Mag-Net
preparations use 10 mM TCEP with 40 mM 2-chloroacetamide, which gives the same
+57.021464 carbamidomethyl adduct as iodoacetamide.

**The file says otherwise.** A bare-cysteine search finds 2,268 cysteine
precursors (18.4%); a carbamidomethyl search of the same file finds 189 (1.9%),
with the precursor masses confirmed CAM-shifted to −8 µTh. And the free-thiol
identifications are not weak matches that a broken FDR let through — they are
slightly *stronger* than the non-cysteine ones:

```
                median quantity   median evidence   frac q<0.001
Cys      2268       2,052,394           4.17           78.8%
non-Cys 10040       1,815,883           4.11           72.7%
```

What is unresolved is *which* of the study's six preparations and 618 LC-MS runs
this file is; MassIVE's anonymisation removes the only clue, and the quoted
chemistry may not be the arm this run came from. What is not in doubt is how the
file must be searched: **cysteine free, which is what the benchmark does.**

**Consequence for the project.** The Astral benchmark is atypical — nearly all
real DIA data are alkylated, as IH1 is. Tuning against it risks fitting an
unusual case, and the free-cysteine correction above matters far less in ordinary
use than its effect size here suggests: on an alkylated library it is a no-op by
construction. IH1 is the more representative of the two benchmarks.

## Finding C: fixing it would not recover the DIA-NN gap

DIA-NN searching our library against Astral (4,984,739 precursors against
DIA-NN's own 4,954,236, so search-space size is not a confound):

```
DIA-NN own lib 12,308   our lib 11,339   shared 10,653
missed by our lib 1,655, of which cysteine-containing 279 = 16.9%
baseline cysteine rate in the reference set              = 18.4%
found ONLY with our lib: 686
recovery: non-Cys 86.29%   Cys 87.70%
```

Cysteine peptides are recovered slightly *better* and are slightly
*under*-represented among the misses. The difference is not significant
(z ≈ 1.8, p ≈ 0.07) but the direction that would incriminate cysteine is
excluded by the point estimate alone.

The reason is straightforward: DIA-NN refits RT against its own identifications
during the search, so a +0.9 min library error is absorbed. **The cysteine RT
bias costs nothing when the downstream engine recalibrates. It costs whatever
the RT window costs when the engine does not** — which is ODIA's situation, and
therefore where the correction below is worth having.

kimi's caveat stands: this table is conditioned on DIA-NN's own reference set
and so cannot see a blind spot both libraries share. It answers "is cysteine
over-represented among our misses" — no — and not "is cysteine costing us in
absolute terms".

## The correction

Gated on a CAM-free library, in rt_norm units, bucketed by cysteine count —
because the effect saturates and a linear n × 0.034 would overcorrect:

| cysteines | offset (rt_norm) |
|---|---:|
| 1 | +0.0343 |
| 2 | +0.0645 |
| 3+ | +0.0580 |

Cross-validated by fitting on one dataset and applying to the other:

```
Astral  cysteine-peptide |residual| p95   3.193 -> 2.416 min   (offset fitted on IH1)
IH1     cysteine-peptide |residual| p95   3.327 -> 2.678 min   (offset fitted on Astral)
Astral  all-peptide p95                   2.540 -> 2.261 min
IH1     all-peptide p95                   2.237 -> 2.209 min
```

Two gradients, two instruments, held out from each other. That is the transfer
evidence kimi asked for, and it is the minimum that should be required.

It remains a compensation for an off-manifold input, not a model fix. The model
fix is fine-tuning on run data, which is per-run and already the plan (doc/28).
The correction matters for the first pass, before any run data exists.

## What is resolved and what is not

Resolved: the alkylation defect; whether ODIA encodes the modification correctly
(it does); whether the RT bias is ours (it is not); whether it explains the
DIA-NN identification gap (it does not).

Not resolved: the identification gap itself — 1,655 misses whose cause is
elsewhere. kimi's ordering is the right one: check presence in our library
first, then the RT residual distribution of the misses, then fragment intensity,
then proximity to the FDR threshold, then charge/length/m-z structure.
