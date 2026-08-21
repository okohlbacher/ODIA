# The decoy model is not the cause, and no shuffled decoy can be

2026-08-21. Both reviewers put "fix the decoy model" ahead of everything else
(review 48): the ruler is 109x wrong at depth 2,000 and nothing can be evaluated
against a null of unknown shape. That branch has now been measured and it is
closed -- not by fixing it, but by establishing that the obvious fix does not
work and why no variant of it can.

## The experiment

`full_v8`: run_full_v5.sh exactly, but the library is targets-only so ODIA
regenerates decoys with `-decoys shuffle` -- composition-preserving (permutes the
interior, termini fixed), DIA-NN 2.x's own default family -- instead of
`mutate`, which maps all 20 residues onto LLLVVLLLLTSSSSLLNDQE: eight residue
masses, no aromatics, no sulfur, no G/A/P.

The hypothesis was that `mutate`'s mass-defect skew puts decoy fragments in less
crowded m/z regions where they absorb less real interference, so they score too
low exactly in the confident tail.

## The result: no difference

                          v5 mutate        v8 shuffle
    nominal 0.1% ->    4.36%  (43.6x)   4.59%  (45.9x)
    nominal 1.0% ->    5.72%  ( 5.7x)   5.66%  ( 5.7x)
    nominal 5.0% ->   10.72%  ( 2.1x)  10.15%  ( 2.0x)
    tail ratio @2,000        108.4x            91.4x
    tail ratio @10,000        24.1x            27.7x

Identical to within the noise on ~70 tail events. And at MATCHED empirical FDP
the two runs are the same run: -1.2% at 5.72%, +1.9% at 7.42%, -0.1% at 10%, all
inside the FDP sigma. **The decoy construction is not the cause.**

Note the instrument mattered here: `fdp_compare.py` uses entrapment for both
arms, so decoys never enter it and it necessarily returned a wash. The decoy
model affects NOMINAL q, and only the nominal-to-empirical mapping tests it.

## The other mechanical explanation is also dead

Winner's curse from candidate multiplicity, checked on v8:

    targets: 679,631 precursors, 2.950 candidates each (max 3)
    decoys : 642,038 precursors, 2.958 candidates each (max 3)

Symmetric. There is a 5.9% OPPORTUNITY asymmetry -- more targets than decoys get
any candidate at all -- but that inflates q by about 6%, not 470%.

## Why no shuffled decoy can work

A decoy is a shuffled or mutated SEQUENCE. Its fragment m/z values are ones that
mostly no real peptide produces, so what it picks up from a crowded window is
essentially random coincidence.

An absent TARGET is a real peptide sequence. Its fragment m/z values coincide
with the fragments of real peptides -- homologues, isotopes, co-eluting species,
modified forms of present peptides. What it picks up is STRUCTURED signal from
things that are actually there.

Those are different null distributions, and the difference is largest exactly in
the tail, where a false positive requires several fragments to agree. `shuffle`
and `mutate` are both sequence permutations and neither changes the property
that matters. The construct that DOES carry real peptide fragment masses while
being absent from the sample is the ENTRAPMENT library -- which is precisely why
the two disagree 100-fold and why entrapment is the correct ruler for this
failure (kimi, review 47).

## What follows

**Stop trying to fix decoys.** Three hours of run time have established that the
decoy method is not the lever, and the argument above says no method in that
family is.

**Entrapment-calibrated reporting becomes the serious option.** kimi called
"report empirical FDP alongside nominal q" partially a trap (review 47): the
curve is a property of sample x library x gate x binary, measured once from ~130
entrapment hits, and a run without an entrapment library cannot recompute it. But
those objections are about SHIPPING A STORED CURVE. They do not apply to
including entrapment IN the searched library and calibrating from it per run,
which is what makes the number a measurement rather than a quotation. The costs
are real and must be measured: entrapment enlarges the search space and changes
the training mixture (vault: *First entrapment measurement* records Astral IDs
754 against 2,078 at r = 3).

**And the sequencing survives the negative.** Both reviewers said the ruler comes
first. It still does -- but the answer is not a better decoy, it is a different
instrument. The detection bucket (42.5%) is next regardless, and the
interference bucket stays behind it.

## Correction to doc/49

Option O8 was "composition-preserving decoys as the shipped default", justified
as the prerequisite for evaluating everything else. **It is measured and it
buys nothing.** It should not ship on FDR grounds. There may be an independent
argument for matching DIA-NN 2.x's default family, but that is a compatibility
argument, not a calibration one, and it was not the argument made.
