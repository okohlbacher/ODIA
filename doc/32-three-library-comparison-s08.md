# Three libraries on S08: DIA-NN predicted, ODIA predicted, DIA-NN empirical

A direct comparison of what two predictors put in a library against what the
instrument actually produced.

## The three

| | what it is | size |
|---|---|---|
| **A** DIA-NN predicted | `--fasta-search --predictor --gen-spec-lib --unimod4` on `bench_reviewed_entrap.fasta` | 4,954,236 precursors |
| **B** ODIA predicted | `DIALibraryGenerator`, default config (CAM), same FASTA | 4,984,739 precursors |
| **C** DIA-NN empirical | `--gen-spec-lib` output of the S08 `--unimod4` run: only precursors identified at 1% FDR, with RT, 1/K0 and fragment intensities refined against the data | **37,193 precursors, 346,531 fragments** |

Both samples and all three libraries are carbamidomethylated — S08 is alkylated
(doc/30), so CAM is correct throughout and no override is involved.

**C is the closest thing to ground truth available**, but it is not neutral: it
was cut FROM A during A's own search. Two consequences run through everything
below and are flagged where they bite.

## Coverage: not a differentiator

Both predicted libraries contain **100.00%** of the 37,193 precursors the
instrument yielded. Neither is losing identifications to a missing entry.

## Precursor m/z

| | median offset | p99 \|diff\| |
|---|---:|---:|
| ODIA | +4.60 µDa | 0.046 mDa |
| DIA-NN | +0.00 µDa | 0.000 mDa |

DIA-NN's zero is **not** agreement — C inherits its m/z from A unchanged, so this
compares a number with itself. ODIA's +4.6 µDa against an independent computation
is 4 ppb at m/z 1000, i.e. nothing.

## Retention time — twice as accurate on a LINE, 1.37x after a monotone fit

Each prediction calibrated against the empirical value by a straight line:

| | Pearson r | residual sd | p95 \|residual\| |
|---|---:|---:|---:|
| ODIA / AlphaPeptDeep | 0.97773 | **8.42** | 18.10 |
| DIA-NN | 0.99494 | **4.03** | 8.03 |

Empirical range −49.8 to 132.9, so DIA-NN's error is 2.2% of the span and ODIA's
4.6%. **DIA-NN's RT model is 2.1× better here**, consistent with the direct
measurement against observed RT on the same run (0.554 min against 1.029 min,
doc/30).

Selection bias favours DIA-NN — these are the precursors its own RT prediction
was good enough to find — but a 2× gap is far larger than that can explain.

**Most of that 2× is calibration shape, not prediction quality.** At the peptide
level (32,595 distinct peptides, charge states collapsed):

| | Spearman rho | Pearson r |
|---|---:|---:|
| ODIA vs DIA-NN | 0.9924 | 0.9761 |
| ODIA vs observed | **0.9936** | 0.9783 |
| DIA-NN vs observed | **0.9956** | 0.9948 |

The two libraries agree on ELUTION ORDER to rho = 0.9924, and ODIA is barely
behind DIA-NN against truth. The gap opens only in Pearson — the signature of a
monotone but non-linear mapping. Replacing the straight line with a monotone
(local-linear) calibration:

| | linear sd | monotone sd | removed |
|---|---:|---:|---:|
| ODIA | 8.224 | **5.382** | **34.6%** |
| DIA-NN | 4.031 | 3.918 | 2.8% |

**The gap narrows from 2.04x to 1.37x.** DIA-NN's iRT scale is already linear
against observed RT and has nothing left to recover; ours is not. Rank
displacement says the same: median 615 places (1.89% of the run) against
DIA-NN's 526 (1.61%) — 17% worse in ordering, not 104%.

So ~35% of ODIA's RT error is available from a monotone recalibration before any
model work, which corroborates doc/26 A7 and puts a number on it.

## Ion mobility — FIXED 2026-08-17; was computed and discarded

S08 is diaPASEF, so C carries **measured** 1/K0 for every precursor.

| | populated | r vs measured | residual sd | as % of mean |
|---|---:|---:|---:|---:|
| DIA-NN predicted IM | 100% | 0.99261 | 0.01574 | **1.56%** |
| ODIA CCS → 1/K0 | 100% *(CCS)* | 0.97560 | 0.02848 | **2.82%** |
| **ODIA IM as shipped** | **0%** | — | — | — |

The middle row is what ODIA's CCS prediction is worth after a Mason–Schamp
conversion at 305 K: 2.8% relative error, about 1.8× DIA-NN's, and unambiguously
useful on an instrument where mobility is a separate separation dimension.
**ODIA shipped `IM = 0` for every precursor** — CCS fully populated, the
conversion never applied. That was deliberate, not an oversight:
`test/ccs_to_mobility.py` argued CCS is a property of the ion and 1/K0 is what
one instrument measures, so the conversion belonged downstream, and that "C is
fitted, not remembered".

Fixed now that there is something to fit against. Fitted on these 37,193
measured values the coefficient is **1039.07** where Mason-Schamp gives
18509/sqrt(305) = **1059.82** — 2.0% apart, with per-charge values of 1045, 1035
and 1031 for z = 2, 3, 4. The form holds, so the textbook constant is used rather
than a number fitted to one file. `DigestParams::derive_ion_mobility` (default
on, in the fingerprint as `;im=`) emits 1/K0 alongside CCS, and
`DIANNLibraryFile::completeMobility` fills whichever of the two a file omits, in
both directions, on load.

## Fragment selection

| | fragments/precursor | share of OBSERVED fragments present | median Jaccard |
|---|---:|---:|---:|
| ODIA | 11.97 | **86.5%** | 0.615 |
| DIA-NN | 11.92 | 100.0% | 0.833 |

DIA-NN's 100% is by construction — C's fragments are a subset of A's — so only
ODIA is genuinely being measured. The real statement is: **ODIA fails to include
13.5% of the fragments the instrument actually produced**, at the same fragment
budget. That is a selection problem, and since selection is by predicted
intensity, it is downstream of the intensity model.

## Fragment intensity — the two models are within noise

Scored against **measured** intensities. Two framings, because the first flatters
ODIA and the second flatters DIA-NN:

**On fragments both libraries contain** (ODIA's natural subset):

| | n | spectral angle | p25 | >0.9 | base peak correct |
|---|---:|---:|---:|---:|---:|
| ODIA | 36,918 | **0.745** | 0.624 | **8.3%** | **54.6%** |
| DIA-NN | 37,193 | 0.731 | 0.617 | 5.6% | 53.7% |

**On every observed fragment**, scoring a missing prediction as zero:

| | n | missing | spectral angle | p25 | >0.9 | base peak correct |
|---|---:|---:|---:|---:|---:|---:|
| ODIA | 37,193 | 13.5% | 0.717 | 0.585 | **6.4%** | 53.5% |
| DIA-NN | 37,193 | 0.0% *(by construction)* | **0.731** | 0.617 | 53.7% |

The gap either way is 0.014–0.015 spectral angle and under a point of base-peak
accuracy. **The two MS2 models are equivalent to within the resolution of this
test.** ODIA's apparent deficit in the second table is entirely the 13.5% of
fragments it does not carry, not the intensities it assigns.

Note what this corrects: on Astral, ODIA scored 0.742 against DIA-NN's
*empirical* library and that read as our model being behind. Here, with both
models measured against the same ground truth, they are level — and 0.73 is
simply where in-silico MS2 prediction sits against real spectra.

## Where this leaves the library

| axis | verdict |
|---|---|
| coverage | tied at 100% |
| precursor m/z | tied |
| fragment intensity | tied |
| **fragment selection** | **ODIA misses 13.5% of observed fragments — all of them ranked out, not un-enumerated** |
| ion mobility | **FIXED** — derived both ways, generator and reader |
| **retention time** | 2.04× on a line, **1.37× after a monotone fit**; rho 0.9936 vs 0.9956 |

Ion mobility is done. Fragment selection is a ranking problem localised to y2+
(see the addendum). RT is two-thirds calibration shape and one-third model.

---

# Addendum: transition completeness per peptide

## Counts

Both libraries cap at 12 and both essentially reach it; the instrument produces
fewer.

| | mean | median | min | at the cap | <8 |
|---|---:|---:|---:|---:|---:|
| ODIA | 11.97 | 12 | 8 | 98.7% | 0.0% |
| DIA-NN | 11.92 | 12 | 6 | 96.9% | 0.0% |
| **observed** | **9.32** | 9 | 4 | 20.4% | 21.7% |

## Agreement between the two libraries

Mean 9.85 transitions shared, 2.12 unique to ODIA, 2.07 unique to DIA-NN;
**median Jaccard 0.714**. Nine or more of the twelve agree for 88.4% of
precursors, and fewer than eight for only 3.3%.

## Completeness against what the instrument produced

**ODIA carries 86.5%** of the observed transitions — per-precursor median 87.5%,
complete for 23.1% of precursors, below half for 0.1%. DIA-NN's 100% is not a
measurement: the empirical library was cut from its own predicted library, so it
cannot fall short by construction, and for the same reason **ODIA's unique picks
can never appear in it** — that comparison is unmeasurable, not unfavourable.

## The missing 13.5% is a RANKING problem, not a coverage problem

Classifying every observed transition ODIA does not carry against ODIA's own
enumeration rules (b/y, fragment m/z 200–1800, charge ≤ min(2, precursor)):

```
  46776  100.0%  ENUMERABLE -- ranked below the cap of 12
     15    0.0%  fragment charge 2 on a charge-1 precursor
```

**Every one of them was a fragment ODIA generated and then discarded.** The
enumeration rules cost 15 transitions out of 346,531. Widening b/y, the m/z
window or the charge limit would gain nothing; the intensity model's ranking is
the whole of it.

## And it is a mild ranking problem, with a sharp tail

| | median observed intensity | share of total observed intensity |
|---|---:|---:|
| kept | 1543 | 94.9% |
| missed | 707 | **5.1%** |

52.6% of the misses sit at observed intensity rank 8 or lower — the weak tail,
where being wrong is cheap. But **2.53% of misses are the observed base peak**:
1,184 precursors, 3.2% of the run, where the single most intense fragment is
absent from our library.

## Where the two disagree: y2+

By ion type and fragment charge:

| | shared | ODIA only | DIA-NN only |
|---|---:|---:|---:|
| b1+ | 19.6% | 37.0% | 38.7% |
| b2+ | 1.4% | **11.2%** | 5.6% |
| y1+ | 66.7% | 48.9% | 17.6% |
| y2+ | 12.3% | **2.9%** | **38.1%** |

The disagreement is not "2+ in general", which is how `LibraryGenerator.h`
records it. It is specific: **ODIA over-picks b2+ (11.2% against 5.6%) and
under-picks y2+ (2.9% against 38.1%)**. Doubly-charged y ions are 38.1% of
everything DIA-NN keeps and we do not — and **60.7% of DIA-NN's unique picks
were actually observed**, against 81.8% for transitions both chose.

So our intensity model systematically demotes y2+ relative to DIA-NN's, and
DIA-NN is right to keep them more often than not. That is a sharper target than
the existing `reserved_doubly_charged` knob, which would promote b2+ as well —
the ion type we already over-pick.
