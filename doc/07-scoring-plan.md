# Phase 3 — peak-group scoring and FDR

Status: design. Extraction landed in `04e3c39` (`-stop_after extract`), so the
input to this stage now exists and is reachable from the tool.

## What this stage has to turn into what

In: chromatograms — for each precursor, ~12 transitions sampled on the run's
cycle axis, from `ChromatogramExtractor`.

Out: for each precursor, a **best peak group** with a set of sub-scores, a
discriminant, and a q-value; and the target/decoy pairing that makes the
q-value mean something.

Three separable pieces, and they fail differently:

1. **Peak-group detection** — where in the RT window did this precursor elute,
   if at all. Candidate generation.
2. **Sub-scores** — numbers describing how good each candidate is. ODIA's own
   work, and where the discriminating power actually comes from.
3. **Classifier + FDR** — combine sub-scores into one discriminant, calibrate
   against decoys, report q-values.

Only (3) depends on outside code. (1) and (2) are unblocked.

## 1. Peak-group detection

Candidates from the summed transition signal within the precursor's RT window:
smooth, find local maxima, take the top few as candidates rather than only the
best. Taking only the maximum makes the classifier's job impossible when the
true peak is second — and gives the decoy distribution nothing to be wrong
about, which quietly deflates the FDR.

Boundaries by descending to a fraction of apex height, capped by the window.

## 2. Sub-scores

The OpenSWATH/pyProphet set is the reference, and it is worth being explicit
about which of these ODIA can compute *honestly* today:

| sub-score | what it says | available now |
|---|---|---|
| `xcorr_shape` | mean pairwise cross-correlation of transition shapes at lag 0 | yes |
| `xcorr_coelution` | mean absolute lag of the cross-correlation maxima | yes |
| `library_corr` | Pearson of observed against library intensities | yes |
| `library_dotprod` | normalised dot product, same pair | yes |
| `intensity_score` | peak-group intensity over the window's total | yes |
| `log_sn` | apex over local background | yes |
| `rt_delta` | observed apex minus predicted RT | **needs a real iRT calibration** |
| `mass_error` | mean/variance of observed minus theoretical fragment m/z | needs per-peak m/z, which the extractor currently discards |
| `im_delta` | observed minus predicted ion mobility | needs CCS→1/K0, deliberately downstream |
| `isotope_corr` | precursor isotope pattern agreement | needs MS1 extraction, not built |

The first six are computable from what `Chromatograms` already holds. The rest
each need something that does not exist yet, and **a sub-score computed from a
placeholder is worse than an absent one** — it looks like signal to the
classifier and is noise, so it will be given weight and will degrade the
discriminant. `rt_delta` in particular must wait for the iRT calibration, since
today the tool spreads the library evenly across the run when none is given.

To compute `mass_error` the extractor has to retain the observed m/z of the
matched peak, not just its intensity. That is a small change to the match inner
loop and a second parallel array; noted rather than done, because it widens the
point arrays by a third and the memory budget is already the extraction stage's
binding constraint.

## 3. Classifier and FDR

**Semi-supervised, as pyProphet does it**: train on the confident targets
against decoys, iterate, then compute q-values from the target/decoy
discriminant distributions. XGBoost as the learner.

**Blocked.** The intended implementation is the pyProphet-like XGBoost scoring
already written at `/Users/kohlbach/Claude/mzPeak/OpenDIAlyzer` — a macOS path,
not reachable from the Linux nodes this project builds on. It needs copying to
somewhere shared (`/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/` is the obvious
place) before it can be used or even read.

What *is* on Ceph is upstream pyProphet 3.0.15 with XGBoost 3.2.0 at
`/ceph/ibmi/abi/oliver/envs/pyprophet` — useful as a reference and as a
cross-check, but it is not the code the decision named.

Two things to settle when it arrives: whether it is a library ODIA calls or a
separate step reading a scores table, and whether it expects OpenSWATH's
`.osw` schema (which would make the sub-score table's column names a contract
rather than a choice).

## Ordering

1. Peak-group detection and the six honest sub-scores; write a scores table.
2. Cross-check those sub-scores against OpenSWATH's on the same input, where
   the definitions are supposed to agree.
3. Classifier and q-values, once the code is reachable.
4. **Entrapment before believing any q-value.** Already on the backlog and a
   prerequisite for ProteoBench: an FDR that has not been checked against a
   known-false population is a number, not a rate.
