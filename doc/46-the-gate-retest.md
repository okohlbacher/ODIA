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
