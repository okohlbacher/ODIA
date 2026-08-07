# A learned precursor pre-filter

Status: design. Measured motivation below; nothing implemented.

## Why a filter and not another sub-score

Our ten sub-scores extract at best a 67% target share in the top 100 against a
48.3% null. Fragment match depth, measured on the engine's reference data,
separates 729.7 against 0.39 true identifications per thousand — a **1,800x
range from one quantity**, before any peak-group scoring exists.

We do not have that signal. Measured on our own S08 run against DIA-NN's 37,247
confident identifications, our `fragment_coverage` gives:

| depth | true per 1000 |
|---:|---:|
| 12 | 16.7 |
| 10 | 23.2 |
| 8 | 20.0 |
| 7 | 29.4 |

Flat, and depth 12 is the worst bin. The reason is that ours is **conditioned on
a candidate we already chose**: it counts fragments above their own median at
that candidate, and in a window where 79% of points are non-zero that is nearly
free. The engine's is **unconditioned** — does there exist *any single spectrum*
in the run where the top library fragments appear together?

That is the difference between a filter and a score, and it is why this belongs
before extraction rather than inside the scorer.

## What it operates on

One pass over the run, before peak-picking and before any RT calibration — it
has to work when no identification yet exists, since supplying the seed is half
its purpose.

For each precursor, over the spectra of its isolation window:

* **`depth`** — the maximum, over spectra, of how many of the top-6 library
  fragments fall within tolerance in that one spectrum. The headline feature.
* **`qualifying_spectra`** — how many spectra reach `depth - 1` or better.
* **`total_matches`** — fragment-spectrum coincidences across the whole run.
  Distinguishes one convincing moment from diffuse noise.
* **`rt_contiguity`** — whether the qualifying spectra are *consecutive cycles*.
  A peptide elutes over a peak; random coincidences scatter. Cheap, and
  orthogonal to everything else here.
* **`intensity_rank_agreement`** — Spearman of the matched fragments' observed
  intensities against their library order, at the best spectrum. Not Pearson:
  at this stage only the ordering is trustworthy.
* **`best_spectrum_rt`** — carried forward, not scored. These are exactly the
  anchors the RT calibration needs, which is why the filter and the calibration
  seed are one job rather than two.

MS1 evidence is deliberately absent: we have no MS1 extraction, and a feature
computed from a placeholder is worse than an absent one.

## The classifier

Reuse the adopted stack (`ODIA::Scoring`, GBT default) rather than inventing a
second one. Same semi-supervised loop, same target/decoy labels, k-fold by
precursor.

The fixed "N of the top-6 in one spectrum" rule is the baseline to beat, and it
is a real baseline: it collapses a 600x dynamic range into one bit, but the bit
is well-placed. **The learned model must be shown to beat it on held-out data
before it replaces it**, on the metric below.

## Three rules that make it safe, and why each is load-bearing

This is selection on data upstream of a target-decoy FDR. Done carelessly it
invalidates the null and every q-value downstream becomes optimistic.

**1. Label symmetry.** Identical code path for targets and decoys, and the same
*number* retained from each class — not the same threshold. A threshold is the
trap: targets clear an evidence bar more often, so the surviving decoys are a
biased, weaker sample, and retained decoys then score systematically lower than
they should. Equal counts make the retained decoy set the *best* decoys, which
biases the FDR conservative — the safe direction.

**2. Cross-fitting against the scorer.** The filter selects on evidence the
scorer later uses. If both see the same fragment-match information, selection
inflates the scorer's apparent separation. Either the filter uses only features
the scorer does not, or it is fitted k-fold so no precursor is filtered by a
model trained on itself. The second is cheaper and stricter.

**3. Report what was dropped.** A filter that silently discards is
indistinguishable from a search that found nothing. The retained fraction and
the depth distribution of the discarded set both go in the run log.

## How to know whether it worked

Not the identification count — that moves for many reasons. The measurement is
the enrichment table above, recomputed on our data:

* **true identifications per thousand retained, by depth**, against DIA-NN's
  q <= 0.01 set. Success is a monotone gradient with a top bin far above the
  base rate. Today ours is flat at 16-29 per thousand; the engine's top bin is
  729.7.
* **retained-set precision at fixed size**: of the top N precursors the model
  keeps, what fraction did DIA-NN identify. Compare directly against the fixed
  depth-6 rule at the same N.
* **decoy retention rate must equal target retention rate** by construction —
  assert it, do not hope for it.

One number worth remembering while reading any of this: of 3,980 target
precursors with a peak group in the current S08 slice, only **84 (2.11%)** are
in DIA-NN's confident set. That is the ceiling. An improvement from 0 to 40 is
half of everything available, not a small number.

---

## MEASURED ON OUR DATA, 2026-08-08: the premise does not survive diaPASEF

`test/tools/odia_prefilter_depth.cpp` computes exactly the `depth` above -- the
maximum, over the spectra of a precursor's isolation window, of how many of its
top-6 library fragments fall within 15 ppm in ONE spectrum. Run over all 32,210
MS2 spectra of S08 against `v6_50k` (50,000 precursors, 738 of them in DIA-NN's
confident set, a 14.8-per-1000 base rate):

    depth   targets    true   true/1000   enrichment
      6      50,000     670        13.4        0.9x

**Every single precursor reaches depth 6.** The feature is saturated and carries
no information at all -- 0.9x enrichment is slightly WORSE than picking at
random.

**Why: on diaPASEF a "spectrum" is not a moment.** S08's frames are
mobility-merged -- `SpectrumSource.h` records frame 1 as 32,570 peaks with 739
m/z descents, i.e. ~600-810 TIMS scans concatenated into one array. Finding six
specific m/z values within 15 ppm somewhere in 32,570 peaks spanning the full
range is near-certain by chance. The unconditioned depth question -- "does there
exist any single spectrum where the top fragments appear together?" -- answers
YES for everything, because the container it searches is ~700 moments, not one.

The 1,800x separation in the section above was measured on the reference
engine's data, which is not ion-mobility-merged. **It does not transfer to
diaPASEF unchanged**, and this was the cheapest possible way to find that out:
one pass over the file, no classifier, no subsystem.

### What the feature has to become

Depth must be computed within a MOBILITY SLICE, not within a frame.
`SpectrumPeaks::ion_mobility` carries the per-peak 1/K0, and the precursor's
expected 1/K0 is in the library, so the natural unit is "the peaks of this frame
within +/-0.025 of this precursor's 1/K0" -- the same window the extractor uses.
That restores "one moment" as the unit and should restore the gradient.

Two consequences for the rest of the design:

* This is no longer free of the mobility calibration. A precursor whose library
  1/K0 is wrong lands in the wrong slice and its depth collapses -- so the
  filter inherits the CCS-conversion scale error, and `best_spectrum_rt` is only
  as good as that.
* On a non-mobility instrument (Astral) the original formulation should work as
  designed. **The filter is therefore instrument-conditional**, and it must be
  re-measured on Astral before either result is generalised -- the same mistake
  the mass calibration made when S08's "narrow window is worse" rule was carried
  to Astral and turned out reversed.

### Also learned

`v6_50k.tsv` contains **no decoy rows** (0 of 50,000). ODIA generates decoys
internally, so the label-symmetry rule -- equal COUNTS retained from each class
-- cannot be checked from the library file and has to be asserted after decoy
generation instead.
