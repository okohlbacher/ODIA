# 78 — The evidence rebuild: strategy and instrument

2026-08-30. Follows doc/77 (the feature plateau). Full design:
`shared/libv2/analysis77/DESIGN_PLAN.md` (+ the four source lanes beside it);
strategy artifact "The Evidence Rebuild"; verified literature review of
2026-08-30 (105 agents, 11 findings, all 3-0 vote confirmations).

## Why

Doc/77 established that the scorer is finished (62.3% vs a 65.2% feature
ceiling that is a plateau) and that the missing ~35% of DIA-NN's IDs are
low-abundance, interference-dominated precursors whose correct peak IS found
but arrives with collapsed library-agreement and co-elution evidence. The
rebuild changes what each candidate brings to the scorer.

## Literature facts that shaped it (verified, cited in the review)

1. DIA-NN's "interference correction" (clip fragment XICs at 1.5*r*reference)
   is QUANT-only and permanently OFF since v1.9 (QuantUMS). It was never in
   scoring. Its scoring robustness = per-fragment evidence vectors (73 scores
   against a data-driven best-fragment reference, not the library) + a tiny
   FFN ensemble + co-located candidate competition. Tier A targets those.
2. Alpha-XIC is the precedent for the learned module: a BiGRU+attention
   co-elution score on raw traces, appended to hand-crafted features in a
   per-run semi-supervised loop; +16.7-49.1% at 1% FDR, largest where
   chromatographic evidence is weakest — exactly the buried class.
3. DreamDIA's 170-trace matrix uses extraction resolutions r/0.45r/0.2r (the
   same ratios as DIA-NN's tight traces) plus MS1 M+1..M+4 AND M-1 — an
   independent validation of the instrument's trace catalog. AlphaDIA's
   interference handling is hard fragment competition (kmax=1 shared).
4. None of these papers validates a learned scorer with entrapment. Our
   entrapment-FDP protocol (r=0.1469 in dn_pred_cam, free) is the
   circularity arbiter and ahead of the field.

## The five layers

L1 evidence widening (per-fragment vectors; MS1 M+1/M+2, shadow -1 Da, heavy
neighbor-window, unfragmented-precursor, tight-tolerance, b/y ladder) ->
L2 interference-robust recomputation (cleaned columns ALONGSIDE raw, never
replacing) -> L3 learned trace representation (Alpha-XIC-class, promoted from
gated spike; same bar) -> L4 one classifier tournament (HistGB/XGBoost/
logistic/LDA/tiny-FFN ensemble) -> L5 candidate competition (DIA-NN's
deletion rule vs AlphaDIA's kmax=1, judged by entrapment FDP).

Binding evaluation ladder: python prototypes on the frozen 69,279-id cohort
(sha256 816b375f..., rng 77); family integrates only if recall@33,330 lifts
>2.0 points over the 65.2 ceiling; recovered controls must not regress >0.5;
entrapment FDP decides; full-run re-validation before any default flips.

## The instrument (this commit)

- `-out_chrom_ids <tsv>`: restricts `-out_chrom` (and `-out_ms1_iso`) to
  listed Precursor.Ids. A listed id emits BOTH its target and decoy block;
  extraction is untouched — only the output shrinks (writer-side filter in
  `writeChromatogramTsv`).
- `-out_ms1_iso <tsv>`: MS1 M/M+1/M+2 traces for the cohort. `Ms1Traces::build`
  gains a charge-aware isotope offset (`mz + k*1.003355/z`, applied before the
  ppm calibration scaling) and a keep-mask with cohort-restricted allocation
  (the full-library dense matrix is 77.6 GB per isotope at 4.99M precursors;
  the cohort build is hundreds of MB). Writer emits non-zero runs with one
  flanking zero. Requires `-out_chrom_ids`, enforced.
- Derived-library trick (zero C++): shadow/tight/heavy trace channels are
  ordinary `-out_chrom -out_chrom_ids` extract-only runs against transformed
  libraries with the plain run's mass calibration BAKED INTO the m/z columns
  (a -1 Da-shifted library would otherwise feed the mass probe nonsense), run
  with `-mass_calibration off` and a pinned `-fragment_ppm`.

Validation: 76/76 tests green on the instrument build (binary 96485934);
mix10k smoke export checked for id-list equality, both decoy flags, and all
three isotopes.
