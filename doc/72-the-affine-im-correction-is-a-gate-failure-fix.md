# The affine 1/K0 correction is a fix for a gate failure, not for the library

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
