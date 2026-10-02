# Our library against DIA-NN's, axis by axis

Answering "is RT the only difference". It is not. RT is the second-largest one,
fragment intensity is the largest, and there is one outright defect.

## Method

Both libraries are in-silico, from the same FASTA
(`bench_reviewed_entrap.fasta`). DIA-NN's was regenerated as a *predicted*
library with `--fasta-search --predictor --gen-spec-lib` and no run file, then
converted to Parquet — because `diann_lib.parquet`, used in earlier comparisons,
is DIA-NN's **output** library, whose RT and fragment intensities were refined
against the Astral data. Comparing our prediction to that is a truth test, not a
library comparison, and conflating the two is what produced the retracted
cysteine finding in doc/30.

30,000 precursors sampled at random from the 4,954,236 both libraries contain —
not the identified subset, which both engines selected and which flatters both.

Our side is `shared/gpu_libgen/odia_for_diann.tsv`, built CAM-free, before the
free-cysteine RT correction of doc/30.

## Membership: no difference

```
ours 4,984,739   DIA-NN 4,954,236   shared 4,954,236
ours only 30,503 (0.6%)   DIA-NN only 0
```

**Our library is a strict superset.** Every precursor DIA-NN generates, we
generate. The 30,503 extra are spread across charge (1+ 21.7%, 2+ 32.3%,
3+ 26.5%, 4+ 19.5%) with no single filter explaining them.

## Masses: no difference

| | median | p99 \|diff\| |
|---|---:|---:|
| precursor m/z | +4.24 µDa | 0.065 mDa |
| fragment m/z | +3.87 µDa | 0.061 mDa |

A constant ~4 µDa offset, which is a rounding convention, not an error.

## Ion mobility: **a defect, ours**

```
IM populated:  ours 0.0%    DIA-NN 100.0%
CCS populated: ours 100.0%
```

The library predicts CCS and never converts it to 1/K0. `human_gpu.parquet` has
CCS non-null for every row and IM null for every row, so this is in the library
generator, not the TSV exporter.

Consequence, measured on the diaPASEF run: in the IH1 search with our library,
DIA-NN's `iIM` — the library mobility — is zero for **100%** of identifications,
against a median of 0.9964 with its own library. We hand a diaPASEF search engine
no mobility dimension at all. It still reached 94.3% of the reference, so DIA-NN
compensates, but this is free information being discarded, and ODIA's own
mobility scoring cannot compensate the same way.

`test/ccs_to_mobility.py` already characterises the conversion; mass and charge
are both in hand.

## Fragment selection: a real difference

```
fragments per precursor   ours mean 11.93 (7-12)   DIA-NN mean 11.81 (5-12)
median Jaccard on fragment identity   0.600
of DIA-NN's fragments ours also has   77.5%
of ours DIA-NN also has               76.7%
```

Both cap at 12, but **about a quarter of each library's fragments are absent
from the other**. The composition is nearly identical by ion type
(b 32.5% / y 67.5% against b 32.3% / y 67.7%) and neither carries neutral
losses, so the disagreement is not about which *kinds* of fragment to keep.

The one systematic gap is charge: **ours 82.4% singly / 17.6% doubly charged,
DIA-NN 75.5% / 24.5%**. We keep 7 points fewer doubly-charged fragments. This is
already documented as a consequence of the intensity models disagreeing rather
than of the selection rule (`LibraryGenerator.h`, `reserved_doubly_charged`), and
the next section is why.

## Fragment intensity: the largest difference

On the fragments both libraries chose:

```
spectral angle   median 0.816   p25 0.699   >0.7 74.8%   >0.9 19.7%
same most-intense fragment                              63.1%
```

The two MS2 models disagree on which fragment is the *base peak* in more than a
third of precursors. That is the single biggest divergence between the libraries
and it is upstream of the fragment-selection difference above: ranking by
intensity is what chooses the twelve.

Against DIA-NN's **empirical** library, on its 12,225 identified precursors, our
predictions score spectral angle 0.742 and pick the same base peak 57.9% of the
time — so the model is not merely different from DIA-NN's, it is measurably
short of the observed spectra. That is a truth test with a caveat: those
precursors are the ones DIA-NN's own predictions were good enough to find.

## Retention time: real, and the second-largest difference

Both in iRT units, fitted across the random sample:

```
Pearson r 0.9549     DIA-NN = 1.1404 x ours - 1.738
residual sd 15.6 iRT   p95 |residual| 29.5 iRT
```

At the Astral gradient's 0.175 min/iRT that is a 2.7 min spread between the two
predictors — larger than either one's error against observed RT (~1.1 min sd,
doc/30). So the two RT models disagree with each other more than either
disagrees with the truth, which is only possible if their errors are largely
independent, and suggests an ensemble would beat both.

By cysteine count the residual is −1.47 / −1.22 / −0.81 iRT for 0 / 1 / 2+, i.e.
flat. The two predictors agree on cysteine — because, as doc/30 establishes, they
share the same cysteine error.

## Summary

| axis | difference |
|---|---|
| membership | none — ours is a strict superset |
| precursor and fragment m/z | none — 4 µDa |
| ion type composition | none |
| **ion mobility** | **ours is empty; a defect to fix** |
| fragment charge mix | 7 points fewer 2+ fragments, downstream of intensity |
| fragment selection | ~23% of each library's fragments absent from the other |
| retention time | r 0.955, 2.7 min spread; errors largely independent |
| **fragment intensity** | **largest: base peak disagrees 37% of the time** |
