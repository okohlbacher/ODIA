# The affine 1/K0 correction is a fix for a gate failure, not for the library

> **CORRECTION, 2026-08-28 — the stated mechanism is wrong; the number is a
> different quantity than the document claims it is.** The document says "One
> factor differs: the library's IM column." It does not. Re-reading the logs:
>
>     arm     fragment mass gate          extraction window          candidate groups
>     base    PASSED (ratio 0.186)        +/-10 ppm centred -10.02        511,337
>     imfix   FAILED (ratio 0.252)        +/-50 ppm centred   0.00      1,057,015
>
> The gate requires ratio <= 0.25; the corrected arm reached 0.252 and fell back
> to `fragment_ppm_uncalibrated` = 50 ppm with no centring.
>
> **This is mediation, not confounding, and the distinction is not pedantic.**
> The gate outcome is DOWNSTREAM of the treatment: the corrected library changes
> which cells pass 1 probes, which changes the mass residuals, which changes the
> gate's verdict. So `1,194 -> 780` remains a valid single-draw estimate of the
> TOTAL effect of shipping the corrected library under the current binary — the
> gate flip is part of what shipping it does. Calling the result "confounded"
> would reframe a real, treatment-caused cost as a measurement artefact, which
> happens to rehabilitate the earlier +32.3% claim. That is the direction this
> project's errors have repeatedly leaned, and it is not licensed here.
>
> What IS wrong is the attribution. "-34.7% because the mobility gate fires" is
> dead: the mobility gate was not the differing mechanism, the mass path was.
> The honest form is: the correction's harm on this run operated THROUGH the mass
> gate, and how much of it survives with the mass path held fixed is unmeasured.
>
> Nor does the near-miss license much. The treatment moved the ratio 0.186 ->
> 0.252, an absolute change of 0.066 — the 0.002 distance to the threshold is the
> small number, and quoting only it hides the larger movement. Without a measured
> standard error on the ratio statistic, "0.252 is a knife-edge accident" is
> assertion, not evidence. What the near-miss does license: this run's downstream
> regime is threshold-sensitive, so the total effect is UNSTABLE, and both -34.7%
> and the earlier +32.3% are single draws of a gate that can flip.
>
> What survives measurement: the loss is not decoy-counting arithmetic. Targets
> and decoys expanded almost identically (2.065x and 2.070x), and at a matched
> decoy budget the corrected arm still ranks far fewer targets:
>
>     decoy budget      base     imfix
>               10     1,182       824
>              100     1,592     1,176
>            1,000     2,737     2,418
>
> Read that with care: the classifier is trained per run, so the two arms' scores
> are not on a common scale. It is a statement about each arm's whole pipeline,
> not about extraction alone.
>
> Running now to decompose it: `shared/libv2/run_immass2x2.sh` (library x pass-1
> mobility scale, mass path pinned identically) and `shared/libv2/run_massregime.sh`
> (2 libraries x 3 mass regimes — narrow-centred, wide-CENTRED, wide-uncentred —
> which supplies the missing base-under-wide cell the decomposition needs, and
> separates the fallback's two simultaneous changes, breadth and mis-centring).
>
> The comparability check that would have caught this is now mechanical:
> `shared/libv2/arm_assert.py`, run before any identification count is read.
>
> Reviewed adversarially by codex 0.149.1 and kimi 0.38.0; both independently
> rejected the word "confounded" and the 0.002 argument, and codex identified the
> missing cell.


An earlier measurement on a 10,000-precursor library found that pre-correcting
the library's ion mobility with `IM' = 0.95323*IM + 0.06916` gained **+32.3%**
identifications at the default aperture, and concluded the aperture sweep had
been measuring a calibration defect. The first half is right. The conclusion
does not survive at production scale.

## Measured

500,000 random targets on the full `S08_diaPASEF.mzpeak`, chosen because it is
the smallest scale at which ODIA's own ion-mobility gate clears its 1.25x
peakedness margin -- it cannot at 10k. Same pinned binary
(`8a48387cca4537eb`), same library draw, same RT map, `-rt_window_pass1 400
-rt_window 300`. One factor differs: the library's IM column.

    arm                pass 1   pass 2    gate form
    base (raw IM)         814    1,194    mz_shaped, 56.5% of MSE removed out of fold
    affine-corrected      519      780    constant

**The correction costs 414 identifications, -34.7%.**

## Why

Where the gate fires, ODIA already removes 56.5% of the mean squared 1/K0
error, fitting an m/z-shaped correction of +0.0175 at 407 Th falling to +0.0152
at 1355 Th for charge 2 -- from 1,613 anchors, 814 target and 799 null, out of
fold over 4 folds.

Pre-shifting the library does not replace that fit; it moves the axis underneath
it, and the gate then fits a SECOND correction on top. The gate form changing
from `mz_shaped` to `constant` is that double-shift showing through: on the
pre-corrected axis the residual no longer has the m/z structure the model is
built to find, so it degenerates to an offset.

The +32.3% at 10k was real and is not withdrawn. It was measured where the gate
FAILED -- "peakedness 4.56 against 4.32, a margin of 1.06x where 1.25x is
required" -- leaving the axis raw. There a manual correction is the only
correction there is.

## What the defect actually is

Not the CCS->1/K0 constant: inverting Mason-Schamp per precursor from the
library's own CCS against DIA-NN's observed 1/K0 gives an implied constant of
18,082 against the formula's 18,509, -2.31% and charge-independent (per-charge
spread 0.9%), and `Library.h:64-75` already documents that gap and keeps the
textbook value deliberately.

Not the aperture: the sweep that appeared to favour +/-0.045 was measuring this
same uncorrected axis.

**It is the gate's reachability.** The 1.25x peakedness margin is a threshold on
a measured statistic, and small libraries cannot supply enough anchors to reach
it -- so those runs extract on an axis with a median error of 0.0235 while large
runs self-correct. The margin exists for a good reason (it asks whether the
mobility agreement is precursor-specific rather than interference), so the fix is
not to lower it blindly but to make it attainable: pool anchors across passes, or
widen pass 1's mobility window so the calibration can see the tail it is
correcting (`-im_window_pass1_scale`, added earlier and still untested at this
scale).

## For the record

Applying this correction in production would have cost a third of the
identifications. It was recommended here on the strength of a 10k measurement,
and the only reason it is not in the tree is that the production-scale arm was
run before shipping it.
