#!/usr/bin/env python
"""The modification parser in scripts/finetune_rt.py, case by case.

There was no test for this at all, and the parser had four defects that only a
table like this finds. The worst: anything not matching the UniMod pattern --
a named modification from Spectronaut or MaxQuant -- passed through verbatim
INTO THE SEQUENCE COLUMN and would have been fine-tuned on as though it were a
peptide. Silent corruption of the training data, in the one function whose
whole job is to prevent that.

Note what the coverage claim does and does not mean: the lookup carries 1,524
UniMod ids, but the only dataset ever run through it contains exactly one
modification. Breadth of table is not breadth of testing, which is why the
cases below are hand-written rather than harvested.
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "ft", os.path.join(HERE, "..", "scripts", "finetune_rt.py"))
ft = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ft)
index = ft.build_unimod_index()

failures = []

# (input, expected sequence, expected mods, expected sites)
CASES = [
    ("PEPTIDEK", "PEPTIDEK", "", ""),
    ("AAC(UniMod:4)DEK", "AACDEK", "Carbamidomethyl@C", "3"),
    ("AAC[UniMod:4]DEK", "AACDEK", "Carbamidomethyl@C", "3"),
    ("PEPS(UniMod:21)TIDEK", "PEPSTIDEK", "Phospho@S", "4"),
    ("PEPTIDEM(UniMod:35)K", "PEPTIDEMK", "Oxidation@M", "8"),
    # A peptide N-terminal modification is Any_N-term, not Protein_N-term.
    # Picking whichever came first in alphabase's table meant this flipped on
    # an upstream reordering.
    ("(UniMod:1)SAMPLER", "SAMPLER", "Acetyl@Any_N-term", "0"),
    # alphabase addresses a C-terminal modification as -1. Writing the peptide
    # length put it one residue past the end.
    ("SAMPLER(UniMod:2)", "SAMPLER", "Amidated@Any_C-term", "-1"),
    # Two modifications on ADJACENT residues. The residue lookup used to read
    # the last written chunk, which is empty after a modification, so the
    # second landed on the first one's site.
    ("AM(UniMod:35)C(UniMod:4)K", "AMCK", "Oxidation@M;Carbamidomethyl@C", "2;3"),
]

for raw, seq, mods, sites in CASES:
    try:
        got = ft.parse_modified_sequence(raw, index)
    except SystemExit as e:
        failures.append(f"{raw!r} was refused: {e}")
        continue
    if got != (seq, mods, sites):
        failures.append(f"{raw!r} gave {got}, expected {(seq, mods, sites)}")

# Inputs that MUST be refused. Guessing at any of these corrupts the training
# set in a way nothing downstream could detect.
REFUSE = [
    ("S(Phospho (STY))EQK", "a named modification, not UniMod"),
    ("PEPTIDE(Oxidation)K", "a named modification, not UniMod"),
    ("PEPTIDEC(UniMod:9999)K", "an unknown UniMod id"),
    ("PEPTIDEW(UniMod:21)K", "phospho on tryptophan, which does not exist"),
]
for raw, why in REFUSE:
    try:
        got = ft.parse_modified_sequence(raw, index)
        failures.append(f"{raw!r} was accepted as {got} -- should be refused ({why})")
    except SystemExit:
        pass

# The bare sequence must always be residues only, whatever the input.
for raw, _, _, _ in CASES:
    seq = ft.parse_modified_sequence(raw, index)[0]
    if not seq.isalpha() or not seq.isupper():
        failures.append(f"{raw!r} produced a non-residue sequence {seq!r}")

for f in failures:
    print(f"  FAIL {f}")
print(f"finetune parser: {len(failures)} failures over {len(CASES)} cases "
      f"and {len(REFUSE)} refusals")
sys.exit(1 if failures else 0)
