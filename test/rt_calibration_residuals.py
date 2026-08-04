#!/usr/bin/env python
"""How much of ODIA's retention-time error a calibration can remove.

Compares three ways of mapping a predicted iRT onto a run's observed retention
time. Two are post-hoc calibrations measurable from data already on disk; the
third, retraining, cannot be done here and is measured separately.

Method matters more than the numbers here, so it is stated up front.

**Evaluation is on held-out peptides.** A calibration fitted and evaluated on
the same precursors reports its own flexibility, not its accuracy -- and a
flexible calibration can drive in-sample residual to nearly zero while being no
better out of sample. The split is by STRIPPED SEQUENCE, so a peptide's charge
states never straddle it: 2+ and 3+ of one peptide elute together, and splitting
between them leaks the answer.

**The evaluation SET is chosen against our own interest.** Residuals measured on
precursors that OUR search identified are selected on our RT being right, which
flatters us. Where both are available the DIA-NN-identified set is used, and
both are reported so the size of that bias is visible rather than argued about.

**The predicted value comes from the LIBRARY, not from the search report.**
DIA-NN's `Predicted.iRT` column is its own run-refined estimate, not the value
it read from the library. Calibrating that against observed RT is circular and
reports 0.04 min for every library alike -- which is how this script was wrong
the first time. The library file is the only honest source for "what did we
predict before seeing the run".

Usage: rt_calibration_residuals.py <label> <library> <report.parquet> [...]
"""
import bisect
import math
import statistics
import sys
import zlib


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


def load_library(path):
    """{(stripped sequence, charge): predicted iRT} for the targets."""
    out = {}
    if path.endswith(".parquet"):
        import pyarrow.parquet as pq
        pf = pq.ParquetFile(path)
        for b in pf.iter_batches(batch_size=300000,
                                 columns=["Modified.Sequence", "Precursor.Charge",
                                          "Decoy", "RT"]):
            d = b.to_pydict()
            for k in range(b.num_rows):
                if d["Decoy"][k]:
                    continue
                out.setdefault((strip_mods(d["Modified.Sequence"][k]),
                                int(d["Precursor.Charge"][k])), float(d["RT"][k]))
        return out
    import csv
    with open(path, newline="") as f:
        r = csv.reader(f, delimiter="\t")
        h = next(r)
        i = {n: h.index(n) for n in ("Modified.Sequence", "Precursor.Charge", "Decoy", "RT")}
        for row in r:
            if row[i["Decoy"]] != "0" or not row[i["RT"]]:
                continue
            out.setdefault((strip_mods(row[i["Modified.Sequence"]]),
                            int(float(row[i["Precursor.Charge"]]))), float(row[i["RT"]]))
    return out


def load_observed(path):
    """{(stripped sequence, charge): observed RT in minutes} at 1% FDR."""
    import pyarrow.parquet as pq
    out = {}
    pf = pq.ParquetFile(path)
    for b in pf.iter_batches(batch_size=200000,
                             columns=["Stripped.Sequence", "Precursor.Charge",
                                      "RT", "Q.Value"]):
        d = b.to_pydict()
        for k in range(b.num_rows):
            if d["Q.Value"][k] > 0.01:
                continue
            out[(d["Stripped.Sequence"][k], int(d["Precursor.Charge"][k]))] = float(d["RT"][k])
    return out


def split(keys):
    """Held-out by stripped sequence, so charge states never straddle."""
    train, test = [], []
    for key in keys:
        h = zlib.crc32(key[0].encode()) & 0xFFFFFFFF
        (test if h % 5 == 0 else train).append(key)
    return train, test


def fit_linear(xs, ys):
    n = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    d = n * sxx - sx * sx
    slope = (n * sxy - sx * sy) / d
    return lambda x: slope * x + (sy - slope * sx) / n


def fit_monotone(xs, ys, knots=200):
    """Piecewise-linear through quantile knots -- a LOESS stand-in.

    Local medians rather than local regressions: the same shape-freedom, no
    bandwidth to tune, and robust to the mis-assigned identifications that a
    1% FDR guarantees are present. Interpolation between knots, and the end
    segments extended, so a test point outside the training range still maps.
    """
    pairs = sorted(zip(xs, ys))
    n = len(pairs)
    step = max(1, n // knots)
    kx, ky = [], []
    for i in range(0, n, step):
        chunk = pairs[i:i + step]
        if len(chunk) < 3:
            continue
        kx.append(statistics.median(p[0] for p in chunk))
        ky.append(statistics.median(p[1] for p in chunk))
    if len(kx) < 2:
        return fit_linear(xs, ys)

    def f(x):
        i = bisect.bisect_left(kx, x)
        if i == 0:
            i = 1
        if i >= len(kx):
            i = len(kx) - 1
        x0, x1, y0, y1 = kx[i - 1], kx[i], ky[i - 1], ky[i]
        if x1 == x0:
            return y0
        return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    return f


def stats(residuals):
    a = sorted(abs(r) for r in residuals)
    return {
        "n": len(a),
        "sd": statistics.pstdev(residuals),
        "mad": statistics.median(a),
        "p95": a[int(0.95 * len(a))],
        "p99": a[int(0.99 * len(a))],
    }


def main(argv):
    for label, lib, report in zip(argv[0::3], argv[1::3], argv[2::3]):
        predicted = load_library(lib)
        observed = load_observed(report)
        data = {k: (predicted[k], observed[k]) for k in set(predicted) & set(observed)}
        path = label
        if len(data) < 1000:
            print(f"{path}: only {len(data)} usable rows, skipping")
            continue
        keys = sorted(data)
        train, test = split(keys)
        tx = [data[k][0] for k in train]
        ty = [data[k][1] for k in train]

        print(f"\n=== {path} ===")
        print(f"{len(data):,} precursors at 1% FDR; "
              f"{len(train):,} train / {len(test):,} held out, split by sequence")

        rows = []
        for name, fit in (("linear", fit_linear), ("monotone (LOESS-like)", fit_monotone)):
            f = fit(tx, ty)
            for label, subset in (("in-sample", train), ("held out", test)):
                res = [data[k][1] - f(data[k][0]) for k in subset]
                s = stats(res)
                rows.append((name, label, s))

        print(f"{'calibration':<24} {'set':<10} {'n':>7} {'sd':>7} {'median|r|':>10} "
              f"{'p95':>7} {'p99':>7}   (minutes)")
        for name, label, s in rows:
            print(f"{name:<24} {label:<10} {s['n']:>7,} {s['sd']:>7.3f} "
                  f"{s['mad']:>10.3f} {s['p95']:>7.3f} {s['p99']:>7.3f}")

        # Saturation lives at the end of the gradient, so a global number can
        # hide it. Deciles of TRUE retention time, held-out only.
        f = fit_monotone(tx, ty)
        obs = sorted((data[k][1], data[k][1] - f(data[k][0])) for k in test)
        print("\nheld-out residual sd by decile of observed RT "
              "(monotone calibration):")
        per = max(1, len(obs) // 10)
        for i in range(10):
            chunk = obs[i * per:(i + 1) * per]
            if len(chunk) < 20:
                continue
            print(f"  decile {i + 1:>2}  RT {chunk[0][0]:6.1f}-{chunk[-1][0]:6.1f} min  "
                  f"n={len(chunk):>6,}  sd={statistics.pstdev([c[1] for c in chunk]):.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
