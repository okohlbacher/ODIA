#!/usr/bin/env python
"""Assert library invariants over a written library, independently of ODIA.

Fragment masses are recomputed from a residue table written here, not from
OpenMS, so this cannot agree with ODIA by sharing its arithmetic. Three rounds
of review found defects that OpenMS-based checks would have missed for exactly
that reason.

Usage: check_invariants.py <library.tsv> [--decoy-table]
Exit 0 if every invariant holds, 1 otherwise.
"""

import csv
import sys
from collections import defaultdict

# Monoisotopic residue masses, independent of OpenMS.
RESIDUE = {
    "G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841,
    "T": 101.04768, "C": 103.00919, "L": 113.08406, "I": 113.08406,
    "N": 114.04293, "D": 115.02694, "Q": 128.05858, "K": 128.09496,
    "E": 129.04259, "M": 131.04049, "H": 137.05891, "F": 147.06841,
    "R": 156.10111, "Y": 163.06333, "W": 186.07931,
}
H2O, PROTON = 18.010565, 1.00727646
LOSS = {
    "noloss": 0.0, "H2O": 18.010565, "NH3": 17.026549,
    "H3PO4": 97.976896, "HPO3": 79.966331, "CO": 27.994915,
}
# Modification deltas, keyed both by UniMod accession (as DIA-NN writes them)
# and by OpenMS name (as ODIA writes them when generating from FASTA).
MOD = {
    "1": 42.010565, "4": 57.021464, "35": 15.994915, "121": 114.042927,
    "Acetyl": 42.010565, "Carbamidomethyl": 57.021464,
    "Oxidation": 15.994915, "GG": 114.042927,
}

# DIA-NN's substitution table; the decoy row stores the target's sequence, so
# reproducing a decoy fragment means applying this.
MUT_FROM = "GAVLIFMPWSCTYHKRQEND"
MUT_TO = "LLLVVLLLLTSSSSLLNDQE"

# Ion-series offsets from the residue sum, before charging.
SERIES = {"b": 0.0, "y": H2O, "a": -27.994915, "c": 17.026549,
          "x": H2O + 25.979265,
          # Classic z = y - NH3, which is what OpenMS's Residue::ZIon gives.
          # The z-dot/z+1 convention (y - NH2) differs by 1.0078 Th.
          "z": H2O - 17.026549}


def tokenise(seq):
    """Split into (residue, mod_delta) pairs.

    The modification NAME must never be scanned for residues -- "UniMod"
    contains an M, so a naive scan invents a methionine and every b ion of an
    N-terminally modified peptide comes out 131 Da heavy. That is the same
    mistake that broke ODIA's own decoy tokeniser, arrived at independently.

    Terminal modifications are folded into the first or last residue's delta,
    which is equivalent for fragment masses: an N-terminal modification appears
    in every b ion and only in the full-length y ion, exactly as residue 0's
    delta does.
    """
    out, nterm, i = [], 0.0, 0
    while i < len(seq):
        c = seq[i]
        if c in "([":
            close = ")" if c == "(" else "]"
            j = seq.index(close, i)
            name = seq[i + 1:j]
            acc = name.split(":")[-1] if name.startswith("UniMod") else name
            delta = MOD.get(acc, MOD.get(name, 0.0))
            if out:
                out[-1] = (out[-1][0], out[-1][1] + delta)
            else:
                nterm += delta
            i = j + 1
            continue
        if c in RESIDUE:
            out.append((c, 0.0))
        i += 1
    if out and nterm:
        out[0] = (out[0][0], out[0][1] + nterm)
    return out


def mutate(tokens):
    """Apply the substitution table the way ODIA does: first and last unmodified."""
    if len(tokens) < 4:
        return None
    n_pos = next((k for k in range(1, len(tokens) - 1) if tokens[k][1] == 0.0), None)
    c_pos = next((k for k in range(len(tokens) - 2, 0, -1)
                  if tokens[k][1] == 0.0 and k != n_pos), None)
    if n_pos is None or c_pos is None:
        return None
    out = list(tokens)
    for pos in (n_pos, c_pos):
        r = out[pos][0]
        idx = MUT_FROM.find(r)
        out[pos] = (MUT_TO[idx] if idx >= 0 else r, out[pos][1])
    return out


def fragment_mz(tokens, ftype, ordinal, charge, loss):
    if ftype not in SERIES or charge <= 0 or ordinal <= 0 or ordinal >= len(tokens):
        return None
    part = tokens[:ordinal] if ftype in "abc" else tokens[len(tokens) - ordinal:]
    # LOSS[loss], not LOSS.get(loss, 0.0). Defaulting to zero made the checker
    # unable to distinguish "loss applied" from "loss ignored" for every label
    # it does not know -- which is exactly the set LossType::Other covers -- and
    # baked a physically wrong mass into the fixtures that import this function.
    neutral = sum(RESIDUE[r] + d for r, d in part) + SERIES[ftype] - LOSS[loss]
    return (neutral + charge * PROTON) / charge


def self_test():
    """Pin the tables against values typed in from a reference.

    The fixtures import fragment_mz() to generate the masses this module later
    verifies, so for target rows the check is circular. These literals are the
    only thing anchoring the tables to reality.
    """
    cases = [
        # (sequence, type, ordinal, charge, loss, expected m/z)
        # Derived by hand from standard monoisotopic residue masses:
        #   b = sum(residues) + proton
        #   y = sum(residues) + H2O + proton
        # PEPTIDEK = P 97.05276, E 129.04259, P, T 101.04768, I 113.08406,
        #            D 115.02694, E, K 128.09496
        ("PEPTIDEK", "b", 2, 1, "noloss", 227.10263),   # P+E
        ("PEPTIDEK", "y", 2, 1, "noloss", 276.15539),   # E+K
        ("PEPTIDEK", "y", 2, 2, "noloss", 138.58133),   # same, doubly charged
        ("PEPTIDEK", "y", 4, 1, "H2O", 486.25583),      # I+D+E+K, water lost
        ("PEPTIDEK", "b", 3, 1, "NH3", 307.12884),      # P+E+P, ammonia lost
    ]
    for seq, ftype, ordinal, charge, loss, want in cases:
        got = fragment_mz(tokenise(seq), ftype, ordinal, charge, loss)
        if got is None or abs(got - want) > 1e-4:
            raise SystemExit(f"self-test failed: {seq} {ftype}{ordinal}^{charge} "
                             f"loss={loss} got {got} want {want}")


def main(path, check_decoys):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    failures = []

    by_precursor = defaultdict(list)
    for r in rows:
        by_precursor[r["Precursor.Id"]].append(r)

    # 1. Precursor.Id must be unique to one precursor, or a reload merges them.
    for pid, rs in by_precursor.items():
        if len({(r["Modified.Sequence"], r["Precursor.Charge"], r["Decoy"]) for r in rs}) > 1:
            failures.append(f"Precursor.Id {pid!r} covers more than one precursor")

    # 2. Rows of a precursor must be contiguous.
    seen, last = set(), None
    for r in rows:
        if r["Precursor.Id"] != last:
            if r["Precursor.Id"] in seen:
                failures.append(f"Precursor.Id {r['Precursor.Id']!r} is not contiguous")
            seen.add(r["Precursor.Id"])
            last = r["Precursor.Id"]

    # 3. Every fragment m/z must be reproducible from the row itself.
    checked = mismatched = 0
    worst = (0.0, "")
    for r in rows:
        tokens = tokenise(r["Modified.Sequence"])
        if r["Decoy"] == "1":
            if not check_decoys:
                continue
            tokens = mutate(tokens)
            if tokens is None:
                continue
        got = float(r["Product.Mz"]) if r["Product.Mz"].strip() else 0.0
        if got == 0.0:
            continue          # unusable m/z, counted separately by the tool
        try:
            expect = fragment_mz(tokens, r["Fragment.Type"],
                                 int(r["Fragment.Series.Number"]),
                                 int(r["Fragment.Charge"]), r["Fragment.Loss.Type"])
        except KeyError as e:
            failures.append(f"unrecognised label {e} on a row carrying m/z {got}")
            continue
        except ValueError:
            continue
        if expect is None:
            # The checker cannot model this row, but the tool wrote a mass for
            # it. Skipping here is what let a mutation that relabels every
            # unknown fragment type as "y" pass the whole suite.
            failures.append(f"row has m/z {got} but type/ordinal/charge "
                            f"({r['Fragment.Type']}/{r['Fragment.Series.Number']}/"
                            f"{r['Fragment.Charge']}) cannot be modelled")
            continue
        checked += 1
        d = abs(expect - got)
        if d > 2e-4:
            mismatched += 1
            if d > worst[0]:
                worst = (d, f"{r['Modified.Sequence']} {r['Fragment.Type']}"
                            f"{r['Fragment.Series.Number']}^{r['Fragment.Charge']} "
                            f"loss={r['Fragment.Loss.Type']} stored={got:.5f} "
                            f"expected={expect:.5f}")
    if mismatched:
        failures.append(f"{mismatched}/{checked} fragment m/z not reproducible; "
                        f"worst {worst[0]:.4f} Th: {worst[1]}")

    # 4. Values that no downstream parser accepts.
    for col in ("RT", "IM", "Precursor.Mz", "Product.Mz", "Relative.Intensity"):
        bad = sum(1 for r in rows if r[col].strip().lower() in ("nan", "-nan", "inf", "-inf"))
        if bad:
            failures.append(f"{bad} rows carry a non-numeric token in {col}")

    print(f"{path}: {len(rows)} rows, {len(by_precursor)} precursors, "
          f"{checked} fragment m/z verified")
    for f in failures:
        print(f"  FAIL {f}")
    return 1 if failures else 0


if __name__ == "__main__":
    self_test()
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    sys.exit(main(args[0], "--decoy-table" in sys.argv))
