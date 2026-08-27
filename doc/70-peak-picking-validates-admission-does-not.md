# Peak picking validates; the deficiency is admission

Measured 2026-08-27 on a 10,000-precursor stratified library (5,000 drawn at
random from DIA-NN's q<=0.01 set on the full 4.99M-target library, 5,000 at
random from the rest, both filtered to fall inside one of the instrument's 24
stated m/z x 1/K0 window boxes), a fixed linear RT map, and one pinned binary.

Three engines, on the 2,553 DIA-NN-confident precursors all three cover.
ODIA's answer is its **top-DScore peak group's RT**; DIA-NN's is its **reported
RT**; OSW's is the nearest of its (unranked) candidates.

    where ODIA ACCEPTS, q<=0.01   n=1,362      median    <=10 s
      ODIA   vs DIA-NN                          0.00 s    94.6%
      ODIA   vs OSW                             0.63 s    93.4%
      DIA-NN vs OSW                             0.75 s    94.1%

    where ODIA DECLINES           n=1,191
      ODIA   vs DIA-NN                         83.11 s    30.1%
      ODIA   vs OSW                            26.37 s    38.1%
      DIA-NN vs OSW                             2.93 s    72.3%

**Where ODIA accepts, its peak choice is indistinguishable from two independent
engines** -- a median of exactly 0.00 s against DIA-NN, and a 10-second
agreement rate (94.6%) that matches what DIA-NN and OSW achieve with each other
(94.1%). Peak selection is not the problem.

**Where ODIA declines, the peak is still there.** DIA-NN and OSW agree with each
other on 72.3% of those same precursors, so a findable peak exists for most of
them; ODIA either never proposes it or does not admit it. That is the whole gap,
and it converges with the project's earlier accounting that Gate C carries ~98.8%
of the emission loss.

## What this replaces

Two claims made earlier the same day, both wrong, both mine:

* *"ODIA's apex sits 131 s from DIA-NN's while DIA-NN and OSW agree to 2.48 s."*
  The 131 s was a measurement artefact. It used the apex of the
  highest-intensity fragment trace (`trace_stats.py`, a raw argmax over a 598 s
  median window) as a proxy for "the peak ODIA chose". ODIA writes its actual
  choice to the `RT` column of `-out` (`OpenDIAlyzer.cpp:4346-4360`, `g.apex_rt`).
  Using it, the same comparison is **1.39 s**.
  DIA-NN's XIC is not a safe proxy for DIA-NN either: it is centred on DIA-NN's
  own answer, so the max-fragment apex misses DIA-NN's reported RT by a median
  13.85 s.

* *"The picker does not degrade with search-window width -- enrichment over a
  same-window chance baseline rises with width."* The test cannot fail. Chance
  scales as 1/W, so enrichment rises unless accuracy falls by more than 3.75x;
  it fell 2.03x. On a fixed population reachable at +/-60 s, top-DScore RT
  against DIA-NN's reported RT goes **71.1% -> 67.2% -> 61.9%** at +/-60, +/-120,
  +/-300 s. The picker does degrade with width.

## The mechanism behind the width effect: generation, not ranking

Of the 203 precursors ODIA gets right at +/-60 s and wrong at +/-300 s, only
**62 (30.5%)** still have a correct peak group anywhere in the +/-300 s candidate
list. **69.5% are never proposed at all.** The winning wrong candidate sits a
median 156.2 s from the window centre while the true answer sits 22.7 s from it,
and these arms ran `-gate_alpha 0`, so nothing penalises distance from the RT
prior. `candidate_min_separation` defaults to 1, so a cap of 3 candidates can be
spent sampling one basin; the option text at `OpenDIAlyzer.cpp:1027` already
records recall-within-cap rising 72.7% -> 75.5% at separation 5.

Ranking is the smaller half: on the 62 where the right group IS present, the
wrong winner's DScore margin is a median 0.657.

## Traps this measurement walked into, for whoever repeats it

* **Select ODIA's row by max DScore, not min QValue.** QValue ties and the first
  row wins, which gives 49.7% instead of 61.1% against OSW on the same file.
* **The +/-300 s arm `acq_*_i020` uses `-precursor_im_window 0.020` while the
  width-sweep arms w060/w120 use 0.025.** Mixing them confounds width with
  aperture by 8.2 points. Use `acq_scores_i025` for any width comparison.
* **OSW candidates are unranked**, so "nearest of 5" is a free minimum rather
  than OSW's answer. It inflates every OSW column equally: DIA-NN against a
  *random* single OSW candidate is 21.2%, against nearest-of-5 it is 84.4%.
* **OSW's coverage is not selection-free** -- it covers 53.9% of the
  DIA-NN-confident stratum against 27.1% of the random one.
* `osw_mix10k.tsv` holds 20,258 features over **4,052** precursors; the 7,744 in
  `osw_extract.log` is ids matched, not precursors with candidates.
