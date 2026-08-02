#!/usr/bin/env python
"""Build assay-library fixtures that deliberately exercise the awkward paths.

Every critical defect found in three rounds of review lived in a path the
DIA-NN fixture cannot reach: it has no neutral losses (305,035/305,035 are
`noloss`), no unusual Parquet types, no nulls, and no sequences that collide
under the decoy substitution table. Each round therefore found the next thing
the fixture happened to hide. These fixtures are the fixture's complement.

Run via scripts/make_fixtures.sh; outputs land on node-local scratch.
"""

import os
import sys

import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "test"))
from check_invariants import fragment_mz, mutate, tokenise  # noqa: E402

# Fixtures must carry physically valid masses. A fixture with invented m/z
# cannot detect a mass bug -- it fails against any correct implementation and
# any incorrect one alike, so it discriminates nothing.

COLUMNS = [
    "Precursor.Id", "Modified.Sequence", "Stripped.Sequence", "Precursor.Charge",
    "Decoy", "RT", "IM", "Precursor.Mz", "Product.Mz", "Relative.Intensity",
    "Fragment.Type", "Fragment.Charge", "Fragment.Series.Number",
    "Fragment.Loss.Type", "Protein.Group",
]


def row(seq, z, decoy, mz, prod=None, loss="noloss", ftype="y", ordinal=3,
        fz=1, rt=10.0, im=0.0, pg="P1", intensity=1.0):
    """One transition row, in DIA-NN's dialect.

    @p prod is computed from the sequence unless given explicitly (which the
    malformed-number fixture does on purpose).
    """
    if prod is None:
        tokens = tokenise(seq)
        # A decoy row carries the target's sequence (DIA-NN's convention), so
        # its fragment masses come from the mutated form. Writing target masses
        # on a decoy row would make the fixture internally inconsistent.
        if decoy:
            tokens = mutate(tokens) or tokens
        prod = fragment_mz(tokens, ftype, ordinal, fz, loss)
        prod = "" if prod is None else round(prod, 5)
    stripped = "".join(c for i, c in enumerate(seq) if c.isalpha() and not _in_mod(seq, i))
    return {
        "Precursor.Id": f"{seq}{z}", "Modified.Sequence": seq,
        "Stripped.Sequence": stripped, "Precursor.Charge": z, "Decoy": decoy,
        "RT": rt, "IM": im, "Precursor.Mz": mz, "Product.Mz": prod,
        "Relative.Intensity": intensity, "Fragment.Type": ftype,
        "Fragment.Charge": fz, "Fragment.Series.Number": ordinal,
        "Fragment.Loss.Type": loss, "Protein.Group": pg,
    }


def _in_mod(seq, i):
    depth = 0
    for k, c in enumerate(seq):
        if c in "([":
            depth += 1
        elif c in ")]":
            depth -= 1
            if k == i:
                return True
            continue
        if k == i:
            return depth > 0
    return False


def write_tsv(path, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in COLUMNS) + "\n")


def main(dest):
    os.makedirs(dest, exist_ok=True)

    # --- losses: every label the parser accepts, plus ones it does not -------
    # HPO3 and H3PO4 differ by one water; aliasing them made every metaphosphate
    # decoy fragment 18.0106 Da wrong for a whole round.
    rows = []
    for i, loss in enumerate(
        ["noloss", "H2O", "NH3", "H3PO4", "HPO3", "CO", "H2O+H2O", "CH3SOH", ""]
    ):
        rows.append(row("PEPTIDEK", 2, 0, 466.7, loss=loss, ordinal=3 + i % 4))
    write_tsv(f"{dest}/adv_losses.tsv", rows)

    # --- modifications: N-terminal, C-terminal, several on one peptide -------
    # N-terminal modifications broke the decoy tokeniser, and 1.3% of the DIA-NN
    # fixture carries one.
    rows = []
    for seq in [
        "(UniMod:1)AVVPASLSGQDVGSFAYLTIK",
        "PEPTIDEC(UniMod:4)K",
        "(UniMod:1)PEPC(UniMod:4)TIDEM(UniMod:35)K",
        "PEPTIDEK(UniMod:121)",
        "C(UniMod:4)C(UniMod:4)C(UniMod:4)PEPTIDEK",
    ]:
        for k in range(4):
            rows.append(row(seq, 2, 0, 700.0, ordinal=k + 2))
    write_tsv(f"{dest}/adv_mods.tsv", rows)

    # --- decoy-sequence collisions ------------------------------------------
    # The substitution table is many-to-one, so distinct targets can mutate onto
    # the same sequence; that made the tool unable to reload its own output.
    rows = []
    for seq in ["HLGSLVNSR", "HIGSLVNSR", "HLGSLVNTR", "YVLSLSLR", "YILSLSLR"]:
        for k in range(5):
            rows.append(row(seq, 2, 0, 500.0 + len(seq), ordinal=k + 2))
    write_tsv(f"{dest}/adv_collisions.tsv", rows)

    # --- fragment types beyond b/y, and out-of-range ordinals ----------------
    rows = []
    for k, ftype in enumerate(["a", "b", "c", "x", "y", "z", "p", "?"]):
        rows.append(row("PEPTIDEK", 2, 0, 466.7, ftype=ftype, ordinal=3))
    rows.append(row("PEPTIDEK", 2, 0, 466.7, 450.0, ordinal=99))   # past the peptide
    rows.append(row("PEPTIDEK", 2, 0, 466.7, 451.0, ordinal=0))    # meaningless
    rows.append(row("PEPTIDEK", 2, 0, 466.7, 452.0, fz=0))         # no charge
    rows.append(row("PEPTIDEK", 2, 0, 466.7, 453.0, fz=-1))        # negative
    write_tsv(f"{dest}/adv_fragments.tsv", rows)

    # --- target/decoy pairs sharing sequence, charge and m/z, as DIA-NN emits -
    rows = []
    for decoy in (0, 1):
        for k in range(4):
            rows.append(row("AILVDLEPGTMDSVR", 2, decoy, 808.4216, ordinal=k + 2))
    write_tsv(f"{dest}/adv_target_decoy.tsv", rows)

    # --- numeric edge cases in the m/z and metadata columns ------------------
    rows = [
        row("BADA", 2, 0, "", 500.0),          # empty precursor m/z
        row("BADB", 2, 0, -500.1, 500.0),      # negative
        row("BADC", 2, 0, 1e9, 500.0),         # beyond the fixed-point range
        row("BADD", 2, 0, "n/a", 500.0),       # unparseable
        row("BADE", 2, 0, 466.7, ""),          # empty product m/z
        row("BADF", 2, 0, 466.7, 500.0, loss="H3PO4", ordinal=1, fz=1),  # loss below zero
        row("GOOD", 2, 0, 466.7, 500.0),       # a healthy control
    ]
    write_tsv(f"{dest}/adv_badnumbers.tsv", rows)

    # --- Parquet: unusual but legal column types, and nulls ------------------
    base = [row("PEPTIDEK", 2, 0, 466.7, ordinal=k + 2) for k in range(4)]
    base += [row("PEPTIDER", 2, 1, 470.7, ordinal=k + 2) for k in range(4)]
    cols = {c: [r[c] for r in base] for c in COLUMNS}

    pq.write_table(pa.table({
        "Precursor.Id": pa.array(cols["Precursor.Id"]).dictionary_encode(),
        "Modified.Sequence": pa.array(cols["Modified.Sequence"], pa.large_string()),
        "Stripped.Sequence": pa.array(cols["Stripped.Sequence"]),
        "Precursor.Charge": pa.array(cols["Precursor.Charge"], pa.uint8()),
        "Decoy": pa.array([bool(v) for v in cols["Decoy"]], pa.bool_()),
        "RT": pa.array(cols["RT"], pa.float32()),
        "IM": pa.array(cols["IM"], pa.float16()),
        "Precursor.Mz": pa.array(cols["Precursor.Mz"], pa.float64()),
        "Product.Mz": pa.array(cols["Product.Mz"], pa.float32()),
        "Relative.Intensity": pa.array(cols["Relative.Intensity"], pa.float64()),
        "Fragment.Type": pa.array(cols["Fragment.Type"]).dictionary_encode(),
        "Fragment.Charge": pa.array(cols["Fragment.Charge"], pa.int8()),
        "Fragment.Series.Number": pa.array(cols["Fragment.Series.Number"], pa.int16()),
        "Fragment.Loss.Type": pa.array(cols["Fragment.Loss.Type"], pa.large_string()),
        "Protein.Group": pa.array(cols["Protein.Group"]),
    }), f"{dest}/adv_types.parquet")

    n = len(base)
    pq.write_table(pa.table({
        "Precursor.Id": pa.array(cols["Precursor.Id"]),
        "Modified.Sequence": pa.array([None] * n, pa.string()),
        "Stripped.Sequence": pa.array(cols["Stripped.Sequence"]),
        "Precursor.Charge": pa.array([None] * n, pa.int64()),
        "Decoy": pa.array([None] * n, pa.bool_()),
        "RT": pa.array([None] * n, pa.float64()),
        "IM": pa.array([None] * n, pa.float64()),
        "Precursor.Mz": pa.array(cols["Precursor.Mz"]),
        "Product.Mz": pa.array(cols["Product.Mz"]),
        "Relative.Intensity": pa.array([None] * n, pa.float64()),
        "Fragment.Type": pa.array([None] * n, pa.string()),
        "Fragment.Charge": pa.array([None] * n, pa.int64()),
        "Fragment.Series.Number": pa.array([None] * n, pa.int64()),
        "Fragment.Loss.Type": pa.array([None] * n, pa.string()),
        "Protein.Group": pa.array([None] * n, pa.string()),
    }), f"{dest}/adv_nulls.parquet")

    # --- rows out of precursor order ----------------------------------------
    rows = []
    for k in range(3):
        rows.append(row("AAAAAAAK", 2, 0, 400.0, ordinal=k + 2))
        rows.append(row("BBBBBBBK", 2, 0, 410.0, ordinal=k + 2))
    write_tsv(f"{dest}/adv_interleaved.tsv", rows)

    # --- protein accessions where one is a prefix of another ----------------
    with open(f"{dest}/adv_proteins.fasta", "w") as fh:
        fh.write(">P12345-2 isoform\nPEPTIDEKAAAAAAAAAAR\n")
        fh.write(">P12345 canonical\nPEPTIDEKAAAAAAAAAAR\n")
        fh.write(">P1 short\nPEPTIDEKAAAAAAAAAAR\n")

    print(f"wrote adversarial fixtures to {dest}")
    for f in sorted(os.listdir(dest)):
        if f.startswith("adv_"):
            print(f"  {f}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
