#!/usr/bin/env python
"""Why do ODIA and DIA-NN keep different fragments for the same precursor?

Two very different causes produce the same symptom, and they need different
fixes:

* **The intensities disagree.** Our model ranks the ions differently, so the
  top twelve differ. The fix would be in prediction.
* **The rules disagree.** Both rank the same way, but one tool never considered
  an ion the other kept -- because of a fragment m/z window, a charge cap, an
  ordinal floor, or an ion series the other does not enumerate. The fix would be
  in the digest parameters, and no amount of better prediction would close it.

This separates them. For every fragment DIA-NN kept that ODIA did not, it asks
whether that ion was even admissible under ODIA's rules. An ion outside our
window was never in the running; an ion inside it was considered and rejected on
predicted intensity, which is the only case that implicates the model.

Usage: analyse_fragment_choice.py <odia.tsv> <diann.parquet> [--sample 20000]
"""
import argparse
import collections
import csv
import math
import sys
import zlib

# ODIA's generation rules, as DigestParams defaults them. Stated here rather
# than inferred, so a mismatch between this and the generator shows up as a
# nonsense result instead of a plausible one.
FRAGMENT_MZ_MIN = 200.0
FRAGMENT_MZ_MAX = 1800.0
MAX_FRAGMENT_CHARGE = 2
SERIES = {"b", "y"}
MAX_FRAGMENTS = 12


def strip_mods(s):
    out, depth = [], 0
    for c in s:
        if c in "([{":
            depth += 1
        elif c in ")]}":
            depth = max(0, depth - 1)
        elif depth == 0 and c.isalpha():
            out.append(c.upper())
    return "".join(out)


def sampled(key, fraction):
    if fraction >= 1.0:
        return True
    return (zlib.crc32(f"{key[0]}/{key[1]}".encode()) & 0xFFFFFFFF) < fraction * 0x100000000


def load_tsv(path, fraction):
    out = {}
    with open(path, newline="") as f:
        r = csv.reader(f, delimiter="\t")
        h = next(r)
        i = {n: h.index(n) for n in
             ("Modified.Sequence", "Precursor.Charge", "Decoy", "Product.Mz",
              "Relative.Intensity", "Fragment.Type", "Fragment.Charge",
              "Fragment.Series.Number")}
        for row in r:
            if row[i["Decoy"]] != "0":
                continue
            key = (strip_mods(row[i["Modified.Sequence"]]),
                   int(float(row[i["Precursor.Charge"]])))
            if not sampled(key, fraction):
                continue
            out.setdefault(key, {})[
                (row[i["Fragment.Type"]].lower(),
                 int(float(row[i["Fragment.Series.Number"]])),
                 int(float(row[i["Fragment.Charge"]])))] = (
                     float(row[i["Relative.Intensity"]]), float(row[i["Product.Mz"]]))
    return out


def load_parquet(path, fraction):
    import pyarrow.parquet as pq
    out = {}
    pf = pq.ParquetFile(path)
    cols = ["Modified.Sequence", "Precursor.Charge", "Decoy", "Product.Mz",
            "Relative.Intensity", "Fragment.Type", "Fragment.Charge",
            "Fragment.Series.Number"]
    for batch in pf.iter_batches(batch_size=300000, columns=cols):
        d = batch.to_pydict()
        for k in range(batch.num_rows):
            if d["Decoy"][k]:
                continue
            key = (strip_mods(d["Modified.Sequence"][k]), int(d["Precursor.Charge"][k]))
            if not sampled(key, fraction):
                continue
            out.setdefault(key, {})[
                (d["Fragment.Type"][k].lower(), int(d["Fragment.Series.Number"][k]),
                 int(d["Fragment.Charge"][k]))] = (
                     float(d["Relative.Intensity"][k]), float(d["Product.Mz"][k]))
    return out


def admissible(ion, mz, precursor_charge):
    """Could ODIA's generator have produced this ion at all?"""
    series, ordinal, charge = ion
    if series not in SERIES:
        return "ion series not enumerated"
    if charge > MAX_FRAGMENT_CHARGE:
        return "fragment charge above the cap"
    if charge > max(1, precursor_charge - 1):
        return "fragment charge >= precursor charge"
    if not (FRAGMENT_MZ_MIN <= mz <= FRAGMENT_MZ_MAX):
        return "outside the fragment m/z window"
    if ordinal < 1:
        return "ordinal below 1"
    return None


def pearson(pairs):
    if len(pairs) < 3:
        return float("nan")
    n = len(pairs)
    mx = sum(p[0] for p in pairs) / n
    my = sum(p[1] for p in pairs) / n
    num = sum((a - mx) * (b - my) for a, b in pairs)
    dx = math.sqrt(sum((a - mx) ** 2 for a, _ in pairs))
    dy = math.sqrt(sum((b - my) ** 2 for _, b in pairs))
    return num / (dx * dy) if dx and dy else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("odia")
    ap.add_argument("diann")
    ap.add_argument("--fraction", type=float, default=0.01)
    args = ap.parse_args()

    sys.stderr.write("loading...\n")
    ours = load_tsv(args.odia, args.fraction)
    theirs = (load_parquet if args.diann.endswith(".parquet") else load_tsv)(
        args.diann, args.fraction)
    shared_keys = sorted(set(ours) & set(theirs))
    sys.stderr.write(f"  {len(ours)} ODIA, {len(theirs)} DIA-NN, {len(shared_keys)} shared\n")

    reasons = collections.Counter()
    theirs_only_series = collections.Counter()
    ours_only_series = collections.Counter()
    both_series = collections.Counter()
    intensity_pairs = []
    rank_of_their_picks = collections.Counter()
    n_shared = n_ours_only = n_theirs_only = 0
    their_pick_intensity_in_ours = []

    for key in shared_keys:
        a, b = ours[key], theirs[key]
        z = key[1]
        # Our own ranking of the ions we kept, so a fragment DIA-NN also kept
        # can be placed within it.
        ours_ranked = sorted(a, key=lambda k: -a[k][0])
        rank = {ion: i for i, ion in enumerate(ours_ranked)}

        for ion in set(a) & set(b):
            n_shared += 1
            both_series[ion[0]] += 1
            intensity_pairs.append((a[ion][0], b[ion][0]))
            rank_of_their_picks[rank[ion]] += 1
        for ion in set(a) - set(b):
            n_ours_only += 1
            ours_only_series[(ion[0], ion[2])] += 1
        for ion in set(b) - set(a):
            n_theirs_only += 1
            theirs_only_series[(ion[0], ion[2])] += 1
            why = admissible(ion, b[ion][1], z)
            reasons[why or "admissible: we ranked it outside our top %d" % MAX_FRAGMENTS] += 1

    total_theirs_only = max(1, n_theirs_only)
    print(f"shared precursors      {len(shared_keys):,}")
    print(f"fragments in both      {n_shared:,}")
    print(f"only ODIA              {n_ours_only:,}")
    print(f"only DIA-NN            {n_theirs_only:,}")
    print()
    print("Why DIA-NN kept a fragment ODIA did not:")
    for why, n in reasons.most_common():
        print(f"  {100 * n / total_theirs_only:5.1f}%  {n:>9,}  {why}")
    print()
    print("Ion series and charge of the fragments each tool kept alone:")
    print("  only ODIA:  ", dict(ours_only_series.most_common(8)))
    print("  only DIA-NN:", dict(theirs_only_series.most_common(8)))
    print()
    r = pearson(intensity_pairs)
    print(f"On the {n_shared:,} fragments BOTH kept, intensity correlation r = {r:.4f}")
    print("  If this is high, the models agree about the ions they both chose,")
    print("  and the difference is in which ions were considered or ranked.")
    print()
    print("Where DIA-NN's shared picks sit in ODIA's own ranking:")
    total_rank = max(1, sum(rank_of_their_picks.values()))
    cum = 0
    for i in range(MAX_FRAGMENTS):
        cum += rank_of_their_picks.get(i, 0)
        print(f"  rank {i + 1:>2}: {100 * rank_of_their_picks.get(i, 0) / total_rank:5.1f}%"
              f"   cumulative {100 * cum / total_rank:5.1f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
