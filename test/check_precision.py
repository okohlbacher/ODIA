#!/usr/bin/env python
"""Assert that m/z survives load-store to the precision D6 claims.

D6 stores m/z as uint32 fixed-point at 1e-5 Th, so a value must round-trip to
within half a quantum. Nothing else in the suite asserts this: a mutation making
toFixed truncate, or the writer emit two fewer significant digits, changed every
mass in the library and passed every test.

Usage: check_precision.py <input.tsv> <output.tsv>
"""
import csv
import sys

TOLERANCE = 6e-6      # half a 1e-5 Th quantum, with room for the decimal text


def load(path):
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        key = (r["Modified.Sequence"], r["Precursor.Charge"], r["Decoy"],
               r["Fragment.Type"], r["Fragment.Series.Number"], r["Fragment.Charge"],
               r["Fragment.Loss.Type"])
        out.setdefault(key, []).append(r)
    return out


def main(src, dst):
    a, b = load(src), load(dst)
    worst = (0.0, "")
    compared = 0
    for key, rows in a.items():
        if key not in b or len(b[key]) != len(rows):
            continue
        for x, y in zip(rows, b[key]):
            for col in ("Precursor.Mz", "Product.Mz"):
                if not x[col].strip() or not y[col].strip():
                    continue
                d = abs(float(x[col]) - float(y[col]))
                compared += 1
                if d > worst[0]:
                    worst = (d, f"{col} {x[col]} -> {y[col]} ({key[0]})")
    if compared == 0:
        print("nothing comparable — the two files share no rows")
        return 1
    print(f"{compared} m/z values compared, worst deviation {worst[0]:.3e} Th")
    if worst[0] > TOLERANCE:
        print(f"  FAIL exceeds {TOLERANCE:.1e} Th: {worst[1]}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
