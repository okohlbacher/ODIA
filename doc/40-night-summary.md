# What the night of 2026-08-18/19 established

One question drove it: **why does ODIA form no candidate at all for 33% of the
precursors DIA-NN identifies from the identical library on the same file?**

The answer turned out to be four layers deep, and three of the first three
answers were wrong. Each was overturned by a measurement closer to the deployed
code than the claim it replaced.

## The chain, as finally measured

```
  extraction        SOUND -- our traces are CLEANER than DIA-NN's
       |
  Gate C            rejects 84% of everything, discriminating nothing
       |            (targets 16.0% admitted, decoys 16.5%)
       |
  removing it       recovers 86.7% of the missing precursors
       |            (candidates 13.6% -> 86.7%, recovered 245 -> 1,236)
       |
  the FDR that      under-calibrated ~12x AT THE OPERATING POINT, because
  judges all of it  constructed decoys do not reach the score tail and
                    real-but-absent peptides do
```

Every identification count quoted tonight -- ours and the sweep's -- was measured
against a q-value that reports 1% and means 11.5%.

## Layer 1: extraction is not the problem

Compared DIA-NN's `--xic` output against ODIA's `-out_chrom` on the SAME 1,000
precursors, the SAME fragments, both aligned on DIA-NN's apex:

| | apex | baseline | prominence |
|---|---:|---:|---:|
| DIA-NN, shared IDs | 16,783 | 1,144 | 0.937 |
| **ODIA, shared IDs** | 6,515 | **1,137** | 0.834 |
| DIA-NN, hard cases | 4,595 | 836 | 0.826 |
| **ODIA, hard cases** | 3,280 | **1,011** | 0.793 |

Target-only, our baseline is 518 vs DIA-NN's 1,144 on shared precursors -- we are
1.9-2.2x CLEANER. Apexes agree to a median of 0.0 s, 96% within 10 s.

Two extraction improvements are real and independent of everything else:
`-mass_width_from_ids apply` (the run measures a 3.62 ppm window and then
discards it; applying it cuts baseline 3x) and the mobility default (confirmed
correct -- widening HURTS).

## Layer 2: Gate C rejects without discriminating

`-gate_log` on a library where tau actually calibrates, 885,045 decisions:

```
  tau = 10.6317      targets admitted 16.0%      decoys admitted 16.5%
```

**Decoys pass more often than targets.** alpha=0.05 nominal delivers 16.5% decoy
admission -- 3.3x off -- because tau is calibrated from the first 20,000 decoys,
which are the earliest-eluting (precursors arrive in RT order), then applied to a
different distribution.

Decoy generation is NOT implicated: decoys share a median 0.0% of fragment m/z
with their targets, and a null built only from zero-overlap decoys moves target
admission by 0.8 points.

## Layer 3: removing the gate works, mechanically

`-gate_alpha 0 -empty_trace_sigma 0` (both flags -- alpha alone falls through to
the excursion gate): candidates for 86.7% of the missing precursors against
13.6%, and 1,236 of them identified against 245. Cost 1.8x wall time, 16.7 GiB,
CPU-bound.

## Layer 4: the FDR, and why the recovery could not be judged

| top N | human | entrapment | decoy | FDP | entrap/decoy vs 0.159 expected |
|---:|---:|---:|---:|---:|---:|
| 1,000 | 992 | 7 | 1 | 4.0% | 44.1x |
| 5,000 | 4,951 | 43 | 6 | 4.9% | 45.2x |
| 14,167 | 13,755 | 282 | 130 | 11.5% | 13.7x |
| 100,000 | 56,180 | 6,423 | 37,397 | 58.7% | 1.1x |

The classifier ranks WELL where it matters -- the top 1,000 are 99.2% human at
4.0% FDP. Its global AUC of 0.529 is dominated by 1.7M mostly-absent precursors
and says nothing about the operating region.

**Decoys and entrapment are exchangeable in the bulk and not in the tail.**
Target-decoy therefore counts 130 false where entrapment implies ~1,614.

## What to do, in order

1. **Calibrate the threshold against entrapment**, or report both. One join over
   a TSV that already exists. Until this is done, no gate or extraction change
   can be evaluated -- tonight's entire sweep was thresholded on a broken
   q-value.
2. **Then** remove Gate C or replace it with the excursion gate. Its removal
   recovers real precursors; that gain was invisible only because the FDP it was
   judged by was itself wrong.
3. `-mass_width_from_ids apply`, pending Astral (doc/13's two-instrument rule).

## Defects found and fixed

| | defect | commit |
|---|---|---|
| 1 | `-out_chrom` applied the RT map twice, emptying the dump | 2a3b3cd |
| 2 | `-out_chrom` merged decoys into their targets under one Precursor.Id | 2438328 |
| 3 | `empty_trace` counter conflated three rejection reasons, so every log ever written reported Gate C rejections as "nothing was extracted" | 17a65fc |
| 4 | `-gate_log` added: the gate's false-negative rate had never been measurable | 9c021aa |

## Claims retracted during the night

Recorded because the pattern matters: every one came from a statistic computed
over a join or a normalisation that had not been validated.

- "Abundance explains the loss" -- a correlate, not a cause.
- "Candidate formation collapses with library size" -- the fixture was below
  `gate_calibration_n`, so the gate was simply off.
- "Our baseline is 1.7-2.8x DIA-NN's" -- backwards; we are 1.9-2.2x cleaner.
- "A constant ~0.13 prominence deficit" -- 0.029-0.033 after removing decoy
  contamination.
- "Gate C rejects 83.6% of true positives" -- the fixture never exercised the
  gate.
- "Ranking is the whole problem" -- corrected by the tail analysis.

## Still open

- `null_calib_` is a file-scope singleton, so tau is calibrated once per process
  and pass 2 inherits pass 1's threshold across different window widths.
  Unmeasured.
- Per-fragment apex spread is 9.0 s against DIA-NN's 1.4 s on shared IDs, with
  66.7% of fragments on-apex against 83.3%. Real, unexplained, and the only
  extraction-side gap left.
- The ODIA + DIA-NN-library arm of the 2x2 was stopped at 7h43, still in pass 1
  at 52% with 170 GiB, against 2h02 for the whole ourlib run. Cause not
  diagnosed.

---

# Adversarial review of the FDR finding, and what it narrows

Codex reviewed the tail analysis. Its verdict, with my audit of the parts that
were cheap to test.

## What survives

**The tail non-exchangeability is established.** At top-1,000, seven of the eight
entrapment-or-decoy hits are entrapments; under the pooled null (p0 = 0.13691)
that is `P(Binom(8, 0.137) >= 7) = 6.4e-6`. The 44.1x point estimate has an exact
95% interval of **5.7x to ~1,990x** -- imprecise upward, but the lower bound
excludes exchangeability decisively. One decoy is enough.

**The entrapment hits are genuinely absent peptides.** Audited on codex's
suggestion:

```
  287 entrapment hits in the top 14,167
  269 distinct stripped sequences        (they do NOT collapse to a few peaks)
    8 also occur as a human sequence     (I/L collapsed) -- physically present
  FDP reported 11.8%  ->  11.4% excluding those 8
```

So codex's alternative #1 (shared or indistinguishable peptides) is worth 0.4
percentage points, and #5 (contamination or carryover concentrating in a few
sequences) is refuted -- 269 distinct sequences from 287 hits, the most repeated
appearing 3 times as charge states.

## What is corrected

**"Exchangeable in the bulk" was sloppy.** AUC 0.510 is not evidence of
exchangeability: its null SE is ~3.6e-4, so 0.510 is many sigma from 0.5. AUC
also measures average pairwise ordering and is largely blind to an extreme-tail
defect. The bulk statement should be dropped; only the tail measurement carries.

**"The classifier is fine" is unsupported.** Good tail behaviour at top-1,000 is
consistent with a classifier that is fine AND with one that is poor but
occasionally lucky. The right claim is narrower: the ranking is USABLE in the
region a threshold occupies.

**"Only threshold placement, not the classifier or candidate formation" does not
follow.** Those are not mutually exclusive, and nothing measured separates them.

## What remains unproven

**Entrapment peptides are demonstrably non-exchangeable with DECOYS. Whether they
correctly model FALSE HUMAN targets is not established**, and that is what an
FDP estimate actually requires. Codex's three surviving alternatives, none
cheaply testable:

- near-human spectral matches (plant-human homologs, near-isobaric precursors
  sharing intense fragments) may occur at a different rate in Arabidopsis than
  among absent human candidates;
- different effective search opportunity -- raw precursor counts do not control
  length, charge, m/z, mobility, window coverage, predicted RT or fragment
  count, and a small covariate shift becomes large after extreme-tail selection;
- false human targets are themselves special (isoforms, homologs, wrong
  charge/modification forms) and may be harder or easier than Arabidopsis. The
  direction is not identifiable from this experiment.

So the honest form of the finding is:

> Constructed decoys and entrapment peptides have severely non-exchangeable
> score tails (>= 5.7x, likely ~44x at top-1,000). Target-decoy therefore
> under-counts false positives at the operating point, by an amount that
> entrapment estimates at ~12x but which is not itself calibrated.

The direction is solid. The magnitude is an estimate with assumptions, and
"11.5% FDP" should be quoted with them.

## Consequence for the plan

Recommendation 1 stands but its wording changes: **report the entrapment FDP
alongside the q-value** rather than replacing the q-value with it. Users get a
number whose assumptions are visible, and the comparison across gate settings
becomes meaningful without asserting that entrapment is ground truth.

The stratified audit codex proposes -- match entrapment and decoys on charge,
length, m/z, RT, IM and fragment count before comparing tails -- is the next
real measurement, and it needs no new run.
