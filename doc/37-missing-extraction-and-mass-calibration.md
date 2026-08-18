# The 33% we never extract, and whether the mass calibration is right

S08, our library, the 2026-08-18 run (14,167 IDs, 2h02, 76.5 GiB). Truth is
DIA-NN's 37,557 confident precursors; 33,749 of them exist in our library.

```
  33,749 truth in library
  22,604 (67.0%)  ODIA scored at least one candidate      -- AVAILABLE
  11,145 (33.0%)  ODIA scored NOTHING                     -- the subject of §1
  11,071 (32.8%)  identified at q <= 0.01
```

## 1. Why 11,145 are never extracted

Tested against the AVAILABLE set as a baseline, so only ENRICHMENT counts.

| reason | missing | available | verdict |
|---|---:|---:|---|
| m/z outside all 24 windows | 0.0% | 0.0% | **not a factor** |
| **charge 1** | **789 (7.1%)** | **10 (0.0%)** | **total loss** |
| library 1/K0 outside the acquired band | 1.1% | 1.1% | not a factor |
| \|1/K0 error\| > 0.05 | 7.5% | 2.7% | **2.8x enriched** |
| \|1/K0 error\| 0.03-0.05 | 16.9% | 16.7% | **not a discriminator** |

The isolation windows span 327.5-1400.6 Th and every truth precursor is inside;
m/z coverage explains nothing. Mobility explains a little, and only in its tail
-- the 0.03-0.05 band is identical between the two groups, so the widely-quoted
"2.8% CCS error" is NOT what loses these precursors.

Those tests leave 67% unexplained. The two that were not tested are:

| | AVAILABLE | MISSING | ratio |
|---|---:|---:|---:|
| **DIA-NN precursor quantity, p10** | 12,762 | 6,403 | **2.0x** |
| **p50** | 43,885 | 18,267 | **2.4x** |
| **p90** | 300,486 | 91,511 | **3.3x** |
| \|RT error\| p50 | 26.7 s | 34.7 s | 1.3x |
| \|RT error\| p90 | 98.6 s | 212.4 s | **2.2x** |
| beyond the pass-1 half-width (+/-70.8 s) | 16.5% | 29.7% | 1.8x |
| **beyond the pass-2 half-width (+/-92.4 s)** | **11.2%** | **24.9%** | **2.2x** |

**Abundance is the dominant reason.** The precursors we never extract are
2.0-3.3x weaker than the ones we do, across the whole distribution rather than
in a tail. These are near-detection-limit peptides that DIA-NN recovers and we
do not: they never reach 3 non-zero points, or never form a peak the picker
accepts.

**Retention time is the clear second, and it is SMALLER than it first looked.**
`rt_window_seconds` is a HALF-width (doc/27 SS1), so pass 2 extracts at
+/-92.4 s, not +/-46.2 s as an earlier draft of this document had it. Against
the correct window, **75.1% of the missing precursors have their true apex
INSIDE the extraction window and still produce no candidate.** The tail is real
-- p90 is 2.2x worse -- but RT can account for at most a quarter of the 11,145.

A check on these numbers: 11.2% of the AVAILABLE set also computes as outside
the window, which is impossible if the map were exact. These residuals use the
SEED map, while pass 2 extracted on the refined one (p50 14.4 s on its anchors),
so both columns overstate the error and only the comparison between groups is
sound. The true RT-attributable share is below 24.9%. (These use the SEED map, so the absolute numbers overstate what
pass 2 saw -- the refined map reaches p50 14.4 s on its own anchors -- but the
comparison between groups is fair because both use the same map.)

Ranked, with what each is worth:

| | reason | scale | fix |
|---|---|---:|---|
| 1 | **abundance / sensitivity** | dominant | lower the bar to forming a candidate; this is a picker and extraction-threshold question, not a calibration one |
| 2 | **RT window tail** | at most 24.9% of missing lie outside +/-92.4 s | the map's tail, not its centre -- a per-precursor window rather than one global width |
| 3 | **charge 1** | 789 | 789 of 11,145 for what is probably a one-line fix |
| 4 | 1/K0 error tail | 838 | only the >0.05 tail matters |

### Charge 1 deserves its own line

**789 missing against 10 available.** Charge-1 precursors are essentially never
extracted, and DIA-NN identifies them. On a diaPASEF run singly-charged ions sit
on a different mobility line and are usually outside the acquisition polygon --
so this may be correct behaviour and the library should simply not carry them,
in which case the fix is to stop generating them rather than to extract them.
Either way it is 7.1% of the loss and is currently silent.

## 2. Is the mass calibration correct? YES -- and there is a proof, not an argument

Pass 1 fitted a CONSTANT **-10.74 ppm at 543.81 Th**, rejecting the m/z shape
because its gain (0.47 ppm) missed the threshold. Centred residual fell
12.885 -> 2.212 ppm (ratio 0.172).

The verification is independent of the fit:

```
fragment window from identifications: 8103 accepted precursors,
  centre 0.08 ppm, per-fragment sigma 1.21 ppm (p95 2.04)
```

**After the correction, the fragments of 8,103 ACCEPTED identifications centre
at 0.08 ppm.** That is a different population, measured a different way, and it
lands on zero. Corroborated over millions of candidates: the per-candidate
deviation sub-score centres at -0.253 ppm (3,070,434 candidates, pass 1) and
-0.226 ppm (5,040,736, pass 2).

The calibration is right to within 0.25 ppm inside a +/-10 ppm window.

### But three things around it are wrong

**a) The seed and the passes calibrate differently, on the same file.**

| | model | at 504 Th | at 1299 Th | residuals |
|---|---|---:|---:|---:|
| CiRT blind search | log_mz, -8.95 ppm + 4.00 ppm per e-fold | -8.95 | **-5.16** | 642 |
| pass 1 and pass 2 | constant | -10.74 | **-10.74** | 5,354 |

**5.6 ppm apart at the top of the range, inside a +/-10 ppm window.** The blind
search therefore hunts the standards at a materially different mass offset than
the passes use, and it is the arm with 8x fewer residuals that chose the more
complex model. Since the seed's anchors decide the RT map, this is not cosmetic.

**b) The gate passed a run in which its own negative control looked better.**
Pass 1: `peakedness 2.15 (control 3.28)`. The m/z-shifted control was MORE peaked
than the data, and the gate passed anyway, with the log conceding it "used to
fail the run and is too noisy to decide anything". The answer is right (see
above), so the gate is not currently doing harm -- but it is not what accepted
the right answer either.

**c) The measured window is still discarded.** Pass 1 concluded a
`suggested half-width 3.62 ppm` from 8,103 accepted precursors; pass 2 extracted
at **+/-10 ppm**. Same defect as doc/34 item 2, now with a better measurement
behind it. A 2.8x tightening of the fragment window is measured and unused, and
narrowing it is the cheapest available purity gain.

## 3. What this changes

The headline is that **our extraction loss is a sensitivity problem, not a
calibration problem.** Mass is right to 0.25 ppm. Mobility's bulk error is not
what loses precursors. RT's centre is fine and only its tail bites. What
separates the 11,145 from the 22,604 is that they are 2-3x dimmer.

That reframes the ranked work from doc/34: narrowing windows and improving
calibration will not recover most of this third. Recovering it means forming
candidates from weaker evidence, which is a picker-threshold question -- and
that must be done together with the FDR work (doc/36, entrapment FDP 11.76%),
because lowering the bar without fixing the null would inflate a number that is
already 12x optimistic.
