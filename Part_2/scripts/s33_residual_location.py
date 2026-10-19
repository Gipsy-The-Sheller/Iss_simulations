"""Are the residuals located at zero?

The question is about location, not about a distributional family.  Whether
the residual is normal, t, or skewed is a separate matter, and for this
pipeline nothing generates a t in any case: the Monte Carlo error is a
binomial count averaged over thousands of replicates, so it is normal; the
theory error is a displacement, not a distribution.  What can be asked, and
what matters, is where the residuals sit.

STANDARDISE FIRST.  The per-point uncertainty varies by more than an order
of magnitude along this sweep, because it is sigma_P / |dP/ds| and dP/ds is
steep on some configurations and nearly flat on others -- in the Farris zone
the transition is shallow enough that a raw-residual average is dominated by
whichever points happened to be worst conditioned.  Dividing each residual
by its own sigma is what makes them comparable, and it is why every point
carries an error bar.

THEN report the mean AND the median, because their disagreement is itself
the diagnostic:

    both near zero      calibrated; no bias
    median 0, mean off  a few points pulling the tail -- look for a
                        specific L or shape where the theory breaks
    both displaced      a systematic offset, i.e. the shape is wrong

Both are tested without a distributional assumption: the sign test (exact
binomial) for the median, a bootstrap for the mean.

POWER.  For the sign test at 80% power and alpha = 0.05, detecting a
departure from an even split needs roughly 47 points at a 70/30 split, 85 at
65/35 and 200 at 60/40.  Below about 20 points the sign test cannot reject
anything and should not be quoted.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "results" / "issc_definition"

NEEDED = {0.70: 47, 0.65: 85, 0.60: 200}


def report(d: pd.DataFrame, label: str) -> dict:
    d = d[np.isfinite(d["iss_err"]) & (d["iss_err"] > 0)]
    if len(d) == 0:
        print(f"  {label}: no usable points")
        return {}
    z = (d["diff"] / d["iss_err"]).to_numpy(float)
    n = len(z)
    mean = float(z.mean())
    med = float(np.median(z))
    boot = np.array([np.random.default_rng(i).choice(z, n, replace=True).mean()
                     for i in range(4000)])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    signs = int((z > 0).sum())
    bt = binomtest(signs, n, 0.5)
    print(f"  {label}   n = {n}")
    print(f"    standardised residuals: mean {mean:+.3f} "
          f"[{lo:+.3f}, {hi:+.3f}] (bootstrap 95%)")
    print(f"                            median {med:+.3f}")
    print(f"                            sd {z.std(ddof=1):.3f}")
    print(f"    sign test: {signs} of {n} positive, p = {bt.pvalue:.3f}")
    if n >= 10:
        w = wilcoxon(z)
        print(f"    Wilcoxon signed-rank: p = {w.pvalue:.3f}")
    else:
        print("    Wilcoxon: too few points")
    if n < 20:
        print("    NOTE: below ~20 points the sign test cannot reject "
              "anything; do not quote it")
    print(f"    worst standardised residual: "
          f"{z[np.argmax(np.abs(z))]:+.2f} at "
          f"{d.iloc[int(np.argmax(np.abs(z)))]['shape']} "
          f"L={int(d.iloc[int(np.argmax(np.abs(z)))]['L'])}")
    return dict(label=label, n=n, mean=mean, median=med, lo=float(lo),
                hi=float(hi), sd=float(z.std(ddof=1)), signs=signs,
                sign_p=float(bt.pvalue))


def main() -> None:
    src = OUT / "issc_definition.csv"
    if not src.exists():
        raise SystemExit(f"no {src}; run s32 first")
    d = pd.read_csv(src)
    d = d[np.isfinite(d["diff"])]
    print("RESIDUAL LOCATION.  diff = Iss.c^MC - Iss.c^pred, on the")
    print("resampling definition of Iss.c.\n")
    rows = []
    print("  every point")
    r = report(d, "all")
    if r:
        rows.append(r)
    for sh in d["shape"].unique():
        sub = d[d["shape"] == sh]
        if len(sub) >= 3:
            print()
            r = report(sub, sh[:34])
            if r:
                rows.append(r)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "residual_location.csv", index=False)

    print("\n  What to look for: a mean far from zero with a median near")
    print("  zero means a few points carry the discrepancy, and the worst")
    print("  standardised residual names them.  A median far from zero")
    print("  means the whole curve is displaced.  Quoting either before")
    print("  the point count reaches the table below is not informative.")
    for frac, need in NEEDED.items():
        print(f"    sign-test power 80% at a {int(frac*100)}/"
              f"{int((1-frac)*100)} split needs about {need} points")


if __name__ == "__main__":
    main()
