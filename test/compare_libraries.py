#!/usr/bin/env python
"""Compare an ODIA assay library against DIA-NN's, precursor by precursor.

Both are predicted from the same FASTA with the same digest settings, so a
disagreement is about prediction, not about which tool was allowed to enumerate
more peptides.

What is compared, and why each one is here:

* **Precursor overlap.** The digest space should be nearly identical. A large
  asymmetry means the enzymology differs, not the prediction, and would make
  every other number below a comparison of different things.
* **Precursor m/z.** Arithmetic, not prediction: both compute it from the same
  residue masses. A systematic offset is a mass-table or proton-mass bug, and
  it must be at the ppm level or below.
* **Retention time**, by rank correlation. The two are in different scales --
  ODIA emits PeptDeep's normalised output, DIA-NN its own iRT -- so agreement
  can only be monotone, not numeric.
* **Fragment selection**, as the overlap of the chosen (type, ordinal, charge)
  triples. Two libraries can agree on every intensity and still share few
  assays if the cap picks differently.
* **Spectral similarity**, as the normalised dot product over the union of
  fragments. This is the number that says whether the intensities are right.
  It is also the check that catches a b/y or ordinal mapping error: mapping a
  predicted spectrum onto the wrong ions gives a plausible-looking library and
  a similarity near zero.

Usage: compare_libraries.py <odia.tsv> <diann.tsv> [--sample N] [--json out]
"""
import argparse
import collections
import csv
import gzip
import json
import math
import os
import random
import sys
import zlib

# Column synonyms. DIA-NN reads many dialects and writes one of them; ODIA
# writes another. Neither is canonical, so both are resolved by synonym.
FIELDS = {
    "sequence": ["ModifiedPeptide", "ModifiedPeptideSequence", "FullUniModPeptideName",
                 "Modified.Sequence", "ModifiedSequence", "PeptideSequence",
                 "Stripped.Sequence"],
    "charge": ["PrecursorCharge", "Precursor.Charge"],
    "precursor_mz": ["PrecursorMz", "Precursor.Mz"],
    "product_mz": ["ProductMz", "Product.Mz"],
    "intensity": ["LibraryIntensity", "RelativeIntensity", "Relative.Intensity",
                  "Library.Intensity"],
    "rt": ["Tr_recalibrated", "RetentionTime", "iRT", "NormalizedRetentionTime", "RT"],
    "im": ["IonMobility", "IM", "Ion.Mobility"],
    "frag_type": ["FragmentType", "Fragment.Type"],
    "frag_num": ["FragmentSeriesNumber", "Fragment.Series.Number"],
    "frag_charge": ["FragmentCharge", "Fragment.Charge"],
    "frag_loss": ["FragmentLossType", "Fragment.Loss.Type"],
    "decoy": ["decoy", "Decoy"],
    "protein": ["ProteinName", "ProteinGroup", "Protein.Group", "UniprotID"],
}


def resolve(header):
    index = {}
    for key, names in FIELDS.items():
        for name in names:
            if name in header:
                index[key] = name
                break
    missing = [k for k in ("sequence", "charge", "product_mz", "intensity") if k not in index]
    if missing:
        raise SystemExit(f"columns not found: {missing}; header is {sorted(header)}")
    return index


def strip_mods(sequence):
    """Bare residues, so two dialects of modification spelling still match.

    Comparing modified sequences directly would compare notation -- (UniMod:4)
    against [+57.0215] against C[Carbamidomethyl] -- rather than chemistry.
    Both libraries here carry only fixed carbamidomethyl, so stripping loses
    nothing and removes a whole class of false disagreement.
    """
    out = []
    depth = 0
    for c in sequence:
        if c in "([{":
            depth += 1
        elif c in ")]}":
            depth = max(0, depth - 1)
        elif depth == 0 and c.isalpha():
            out.append(c.upper())
    return "".join(out)


def rows_of(path):
    """(header, row iterator) for a TSV or a Parquet library.

    DIA-NN 2.x writes Parquet whatever extension is asked for, so both have to
    be read. Parquet is streamed a row group at a time: the reference library
    is 25.6 M rows and materialising it whole costs several gigabytes for no
    reason.
    """
    if path.endswith(".parquet"):
        import pyarrow.parquet as pq
        pf = pq.ParquetFile(path)
        header = [f.name for f in pf.schema_arrow]

        def stream():
            for batch in pf.iter_batches(batch_size=200000):
                columns = [batch.column(i).to_pylist() for i in range(batch.num_columns)]
                for i in range(batch.num_rows):
                    yield [c[i] for c in columns]
        return header, stream()

    opener = gzip.open if path.endswith(".gz") else open
    handle = opener(path, "rt", newline="")
    reader = csv.reader(handle, delimiter="\t")
    header = next(reader)
    return header, reader


def in_sample(key, fraction):
    """Deterministic, content-based subsampling applied identically to both.

    Taking every Nth row would sample differently on two files with different
    row orders. Hashing the key means the same precursors are drawn from both,
    which is the only way a subsample can be compared at all.
    """
    if fraction >= 1.0:
        return True
    h = zlib.crc32(f"{key[0]}/{key[1]}".encode()) & 0xFFFFFFFF
    return h < fraction * 0x100000000


def load(path, want=None, fraction=1.0):
    """{(sequence, charge): {rt, im, mz, protein, frags {(t,n,z): intensity}}}."""
    precursors = {}
    header, reader = rows_of(path)
    if True:
        col = resolve(set(header))
        at = {k: header.index(v) for k, v in col.items()}
        for row in reader:
            if len(row) <= max(at.values()):
                continue
            if "decoy" in at and str(row[at["decoy"]]) not in ("0", "False", "false", ""):
                continue
            key = (strip_mods(row[at["sequence"]]), int(float(row[at["charge"]])))
            if want is not None and key not in want:
                continue
            if not in_sample(key, fraction):
                continue
            entry = precursors.get(key)
            if entry is None:
                entry = precursors[key] = {
                    "rt": float(row[at["rt"]]) if "rt" in at and row[at["rt"]] else None,
                    "im": float(row[at["im"]]) if "im" in at and row[at["im"]] else None,
                    "mz": float(row[at["precursor_mz"]]) if "precursor_mz" in at else None,
                    "protein": row[at["protein"]] if "protein" in at else "",
                    "frags": {},
                }
            loss = row[at["frag_loss"]] if "frag_loss" in at else "noloss"
            ion = (row[at["frag_type"]].lower() if "frag_type" in at else "?",
                   int(float(row[at["frag_num"]])) if "frag_num" in at else 0,
                   int(float(row[at["frag_charge"]])) if "frag_charge" in at else 1,
                   "noloss" if loss in ("", "noloss", "None", "none", "0") else loss)
            entry["frags"][ion] = float(row[at["intensity"]])
    return precursors


def dot_similarity(a, b):
    """Normalised dot product over the union of the two fragment sets.

    Over the union, not the intersection: scoring only shared fragments would
    reward a library that selected two fragments confidently over one that
    selected twelve, which is backwards.
    """
    keys = set(a) | set(b)
    va = [a.get(k, 0.0) for k in keys]
    vb = [b.get(k, 0.0) for k in keys]
    na = math.sqrt(sum(x * x for x in va))
    nb = math.sqrt(sum(x * x for x in vb))
    if na == 0 or nb == 0:
        return 0.0
    return sum(x * y for x, y in zip(va, vb)) / (na * nb)


def spearman(pairs):
    """Rank correlation. The two RT scales differ, so only order can agree."""
    if len(pairs) < 3:
        return float("nan")

    def ranks(values):
        order = sorted(range(len(values)), key=lambda i: values[i])
        r = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
                j += 1
            mean = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = mean
            i = j + 1
        return r

    xs = ranks([p[0] for p in pairs])
    ys = ranks([p[1] for p in pairs])
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return num / (dx * dy) if dx and dy else float("nan")


def quantiles(values, qs=(0.05, 0.25, 0.5, 0.75, 0.95)):
    if not values:
        return {q: float("nan") for q in qs}
    s = sorted(values)
    return {q: s[min(len(s) - 1, int(q * len(s)))] for q in qs}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("odia")
    ap.add_argument("diann")
    ap.add_argument("--sample", type=int, default=0,
                    help="compare only this many shared precursors (0 = all)")
    ap.add_argument("--fraction", type=float, default=1.0,
                    help="load only this fraction of precursors, chosen by a hash "
                         "of the key so both files yield the same ones")
    ap.add_argument("--json", default="")
    args = ap.parse_args()

    sys.stderr.write("loading ODIA library...\n")
    odia = load(args.odia, fraction=args.fraction)
    sys.stderr.write(f"  {len(odia)} target precursors\n")
    sys.stderr.write("loading DIA-NN library...\n")
    diann = load(args.diann, fraction=args.fraction)
    sys.stderr.write(f"  {len(diann)} target precursors\n")

    shared = sorted(set(odia) & set(diann))
    only_odia = len(odia) - len(shared)
    only_diann = len(diann) - len(shared)

    keys = shared
    if args.sample and args.sample < len(shared):
        random.seed(20260803)
        keys = random.sample(shared, args.sample)

    similarities, frag_overlaps, mz_ppm, rt_pairs, im_pairs = [], [], [], [], []
    top1_agree = 0
    for key in keys:
        a, b = odia[key], diann[key]
        similarities.append(dot_similarity(a["frags"], b["frags"]))
        shared_frags = set(a["frags"]) & set(b["frags"])
        union = set(a["frags"]) | set(b["frags"])
        frag_overlaps.append(len(shared_frags) / len(union) if union else 0.0)
        if a["frags"] and b["frags"]:
            if max(a["frags"], key=a["frags"].get) == max(b["frags"], key=b["frags"].get):
                top1_agree += 1
        if a["mz"] and b["mz"]:
            mz_ppm.append(1e6 * (a["mz"] - b["mz"]) / b["mz"])
        if a["rt"] is not None and b["rt"] is not None:
            rt_pairs.append((a["rt"], b["rt"]))
        if a["im"] and b["im"]:
            im_pairs.append((a["im"], b["im"]))

    n = len(keys)
    report = {
        "odia_precursors": len(odia),
        "diann_precursors": len(diann),
        "shared": len(shared),
        "only_odia": only_odia,
        "only_diann": only_diann,
        "jaccard": len(shared) / len(set(odia) | set(diann)),
        "compared": n,
        "mean_fragments_odia": sum(len(odia[k]["frags"]) for k in keys) / n if n else 0,
        "mean_fragments_diann": sum(len(diann[k]["frags"]) for k in keys) / n if n else 0,
        "spectral_similarity": quantiles(similarities),
        "mean_spectral_similarity": sum(similarities) / n if n else 0,
        "fragment_jaccard": quantiles(frag_overlaps),
        "base_peak_agreement": top1_agree / n if n else 0,
        "precursor_mz_ppm": quantiles(mz_ppm),
        "rt_spearman": spearman(rt_pairs),
        "rt_pairs": len(rt_pairs),
        "im_spearman": spearman(im_pairs),
        "im_pairs": len(im_pairs),
    }

    def pct(x):
        return f"{100 * x:.1f}%"

    print(f"precursors      ODIA {report['odia_precursors']:,}  "
          f"DIA-NN {report['diann_precursors']:,}  "
          f"shared {report['shared']:,} (Jaccard {pct(report['jaccard'])})")
    print(f"                only ODIA {only_odia:,}   only DIA-NN {only_diann:,}")
    print(f"compared        {n:,} precursors")
    print(f"fragments/prec  ODIA {report['mean_fragments_odia']:.2f}  "
          f"DIA-NN {report['mean_fragments_diann']:.2f}")
    q = report["spectral_similarity"]
    print(f"spectral dot    median {q[0.5]:.3f}  "
          f"(5% {q[0.05]:.3f}, 25% {q[0.25]:.3f}, 75% {q[0.75]:.3f}, 95% {q[0.95]:.3f})")
    q = report["fragment_jaccard"]
    print(f"fragment overlap median {q[0.5]:.3f}  (5% {q[0.05]:.3f}, 95% {q[0.95]:.3f})")
    print(f"base peak agree {pct(report['base_peak_agreement'])}")
    q = report["precursor_mz_ppm"]
    print(f"precursor m/z   median {q[0.5]:+.4f} ppm  "
          f"(5% {q[0.05]:+.4f}, 95% {q[0.95]:+.4f})")
    print(f"RT rank corr    {report['rt_spearman']:.4f} over {report['rt_pairs']:,} pairs")
    if report["im_pairs"]:
        print(f"IM rank corr    {report['im_spearman']:.4f} over {report['im_pairs']:,} pairs")
    else:
        print("IM              absent from at least one library")

    if args.json:
        with open(args.json, "w") as f:
            json.dump({k: (v if not isinstance(v, dict) else {str(a): b for a, b in v.items()})
                       for k, v in report.items()}, f, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
