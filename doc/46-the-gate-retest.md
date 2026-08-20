# The Gate C retest: prominence wins at matched FDP, and I read it wrong twice first

2026-08-20. `run_full_v6.sh` against `run_full_v5.sh`. Both carry the fragment
floor, both take the same explicit RT map (RT = 473.77 + 1086.50 x libRT, which
v3 and v5 fitted identically), both pin pass 1 to 75.4069 s, neither runs the
mobility seed. `-gate_mode` is the only difference.

## At each run's own nominal 1%

                          IDs  entrap     FDP   est. true  DIA-NN rec.          cost
    v5  Gate C         13,268     130   5.72%      12,387       11,524  3h00 / 124GB
    v6  prominence     11,007      96   5.08%      10,356        9,521  4h08 / 131GB

Read alone, that says the gate replacement loses badly. It is the wrong read.
The two arms do not sit at the same empirical FDP -- prominence's nominal 1%
lands at 5.08% where Gate C's lands at 5.72% -- so this compares a stricter
threshold against a looser one.

## At MATCHED empirical FDP

    target FDP   Gate C IDs   prominence IDs   Gate C rec   prom rec
         5.08%       12,299           10,936       10,872      9,544
         5.72%       13,140           13,877       11,525     11,930
         7.42%       14,484           16,173       12,503     13,572   <- DIA-NN's point
        10.00%       15,700           18,461       13,205     14,910

Prominence wins everywhere at or above 5.72%: +11.7% identifications and +8.6%
DIA-NN recovery at DIA-NN's own operating point, +17.6% at 10%.

## Two wrong readings, recorded because the pattern matters

1. **From pass 1.** v6's pass 1 identified 8,046 against v5's 8,217 while
   scoring 13.3M peak groups against 2.8M, and I wrote that the gate "still
   fails". Pass-1 counts had ALREADY been shown misleading hours earlier -- the
   mobility seed was +43% in pass 1 and +2% at the end (doc/43) -- and I used
   them anyway.
2. **From the nominal threshold.** The final table above, read at nominal 1%,
   says prominence loses 17% of estimated true identifications. At matched FDP
   it gains 12%. codex pre-registered the matched-FDP rule rounds ago and it is
   the rule this project keeps rediscovering.

## What it does not do

16,173 at DIA-NN's 7.42% is still 2.4x short of DIA-NN's 39,211. The gate is
worth having and is not the answer.

## What to do with it

Not simply flip the default. Prominence's q-values are mis-calibrated relative
to Gate C's -- that is what puts its nominal 1% at 5.08% -- so a user asking for
1% FDR today gets FEWER identifications from it, which is how it was misread
above. The sequence is: fix the q-value calibration under the prominence gate so
nominal tracks empirical, THEN compare at nominal 1%, and only then consider the
default. Adopting it now would ship a gate that is better on a curve nobody
reads and worse at the number everybody reads.

Cost is also real: 4h08 against 3h00 and 131 GB against 124 GB, for 4.7x the
peak groups scored.

## The deeper problem the retest exposed: nominal q does not mean anything

Chasing why prominence's nominal 1% lands at 5.08% while Gate C's lands at
5.72%, the answer turned out not to be about the gates at all. Mapping nominal
q-value to entrapment-measured FDP:

    v5  Gate C
      nominal q       IDs   empirical FDP   inflation
           0.1%     7,811           4.36%       43.6x
           1.0%    13,138           5.72%        5.7x
           5.0%    16,061          10.72%        2.1x
          10.0%    17,818          14.85%        1.5x
      nominal needed for a TRUE 1% FDP: 0.00029  ->    117 IDs

    v6  prominence
           1.0%    10,911           5.08%        5.1x
           5.0%    17,672           8.66%        1.7x
          10.0%    20,575          13.59%        1.4x
      nominal needed for a TRUE 1% FDP: 0.00140  ->    256 IDs

Two things, both worse than the identification gap.

**We cannot deliver a true 1% FDR at any usable depth.** A genuine 1% empirical
FDP costs everything: 117 identifications under Gate C, 256 under prominence.
Every number this project has reported at "1% FDR" is really 5%.

**The failure is concentrated in the tail.** 43.6x inflation at nominal 0.1%,
5.7x at 1%, 1.5x at 10%. The decoy null is least trustworthy exactly where
confident identifications are claimed. This is the 44x entrapment/decoy tail
ratio from doc/43 seen as a calibration curve rather than a single number, and
it is the same mechanism: shuffled decoys do not reproduce the errors that
target-labelled absent precursors actually make.

So prominence's q-values are not "mis-calibrated relative to Gate C's" as
written above -- prominence is the BETTER calibrated of the two (5.1x against
5.7x). Both are broken by roughly the same factor, and the earlier wording
should be read as: at nominal 1% the two arms sit at DIFFERENT empirical FDPs,
so comparing them there compares two different thresholds.

**It also reframes tonight's purity claim.** "5.56% against DIA-NN's 7.42%" is
two tools' nominal 1% each delivering 5-7x the advertised error. We are less
bad; we are not good. Wen et al. 2025 records the same across DIA tools, so this
is a field-wide failure rather than an ODIA defect -- but ours is now measured,
and the measurement should be quoted whenever a 1% number is.

The next correctness work is here, not in the identification count: a null that
tracks the tail. Candidates, none tested: entrapment-calibrated q-values as a
shipped output; a decoy construction whose fragments collide with real
co-eluting peptides the way absent targets' do; or reporting empirical FDP
alongside nominal q so the gap is visible rather than implied.

## The failure is stratified, and all three axes agree

v5's 13,268 accepted at nominal 1%, split by covariates. `r_s` is the
entrapment ratio WITHIN each stratum, so composition is controlled:

           charge   accepted  entrap     r_s      FDP
                2      9,956      94  0.1731    5.51%
                3      3,047      29  0.1732    5.55%
                4        262       7  0.1723   15.93%

        fragments   accepted  entrap     r_s      FDP
             7-11      1,292      30  0.1774   13.40%
               12     11,946      99  0.1713    4.88%

         prec m/z   accepted  entrap     r_s      FDP
          200-400        136       4  0.1731   17.51%
          400-600      3,735      62  0.1728    9.77%
          600-800      6,127      45  0.1742    4.25%
         800-1000      2,743      17  0.1735    3.59%
        1000-1200        484       1  0.1713    1.21%

`r_s` is flat at ~0.173 everywhere, so these are real differences in error rate,
not entrapment composition. And the three axes are correlated -- low-m/z
precursors carry fewer in-range fragments, and charge 4 concentrates at low m/z
-- so this is plausibly ONE axis: the less evidence a precursor carries, the
higher its true error, and a single global threshold cannot see it. The
4.88% / 13.40% split by fragment count is the same statement as the 3.59% /
17.51% split by m/z.

kimi advised against a fragment-count-stratified FDR last round -- one band, one
run, n ~ 1,300, "replicate the band structure before revisiting". That advice
was right for the evidence then. This is three independent stratifications
agreeing, with 3,871 identifications in the low-m/z corner alone, which is a
different situation; the replication kimi asked for should still be done on the
Astral run before anything ships.

The fix is standard and is not a new invention: group-wise or covariate-aware
FDR, which the Percolator family already does. Each stratum gets its own null
rather than borrowing the bulk's. On these numbers that reprices ~5,400
identifications currently admitted at 2-3x the headline error rate.
