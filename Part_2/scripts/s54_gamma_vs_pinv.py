"""Does rate heterogeneity move Iss.c, and does it matter whether it is +G or +I?

Part I states both of the following:

  * Proposition 3: Iss.c is a strictly monotone function of the Gamma shape
    alpha, so the tabulated values correspond to alpha = infinity alone.
  * Figure 9 caption: a fifteen-fold change in the critical tree length has
    "no visible effect" on Iss.c, because "the index is nearly at its
    ceiling, not because it is well calibrated".

Those are not contradictory in the strict sense -- a monotone function can be
flat -- but the second offers a mechanism that Part II rejects: Iss.c is set
by the KL distance to the star, and nothing in that account makes it
insensitive to the substitution process.

There is a sharper reason to distrust the asymmetry. Invariable sites lower
Iss, because an invariant site has zero entropy and a fraction P of them
scales Hbar by (1-P). Gamma heterogeneity raises Iss, because h is concave
and spreading the rate over sites raises mean entropy by Jensen. Opposite
signs, but both must act on Iss.c.

And the two are the same axis: a two-category discrete Gamma at small alpha
*is* a +I model, with one near-zero rate carrying most of the mass and the
other rate compensating to keep the mean at one. A claim that one moves
Iss.c and the other does not is therefore a claim that the same axis behaves
differently at its two ends.

Both are computed here with the machinery of Part II, exactly. For each
process the pattern distribution is built as a mixture over rate categories,
the two wrong topologies are projected onto the star as in s36, and Iss.c is
located as the falling crossing of 95% recovery, the definition s47b shows is
the correct one.

Outputs
    results/gamma_pinv/gamma_pinv.csv
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import gammainc
from scipy.stats import gamma as G

HERE = Path(__file__).resolve().parent
PART2 = HERE.parent
ROOT = PART2.parent
sys.path.insert(0, str(ROOT / "src"))

OUT = PART2 / "results" / "gamma_pinv"

_spec = importlib.util.spec_from_file_location("s36", str(HERE / "s36_factors.py"))
s36 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s36)
m = s36.m

_s34 = importlib.util.spec_from_file_location("s34", str(HERE / "s34_iss_null_path.py"))
s34 = importlib.util.module_from_spec(_s34)
_s34.loader.exec_module(s34)

FARRIS = (1.0, 1.0, 0.05, 0.05, 0.02)
L = 30_000
LABS = ("12|34", "13|24", "14|23")

# the four all-same patterns, one base each
CONST = np.zeros(256)
for _d in range(4):
    CONST[85 * _d] = 0.25


def discrete_gamma(alpha, k):
    """Category mean rates, mean 1, Yang (1994) mean discretisation.

    For X ~ Gamma(alpha, scale=1/alpha) the first moment within a bin is
    E[X 1_{lo<X<hi}] = P(alpha+1, alpha*hi) - P(alpha+1, alpha*lo), with P the
    regularised lower incomplete gamma; the rate is k times that, because each
    bin carries mass 1/k.
    """
    if not np.isfinite(alpha):
        return np.ones(1)
    inner = G.ppf(np.arange(1, k) / k, a=alpha, scale=1.0 / alpha)
    edges = np.r_[0.0, inner, np.inf]
    cum = gammainc(alpha + 1.0, alpha * edges)
    return k * np.diff(cum)


def pattern(shape, alpha=np.inf, pinv=0.0, k=8):
    """The quartet pattern distribution under a rates-across-sites model."""
    if pinv > 0:
        p = (1.0 - pinv) * m.pattern_batch(
            tuple(v / (1.0 - pinv) for v in shape), m.SPLITS["12|34"])[0]
        return p + pinv * CONST
    rates = discrete_gamma(alpha, k)
    P = np.zeros(256)
    for r in rates:
        P += m.pattern_batch(tuple(v * r for v in shape),
                             m.SPLITS["12|34"])[0]
    return P / len(rates)


def fits_from_p(p, pend4):
    out = [p]
    for lab in LABS[1:]:
        th = s36.star_fit(p, m.SPLITS[lab], pend4)
        out.append(m.pattern_batch(th, m.SPLITS[lab])[0])
    return out


def recovery(shape, L, **kw):
    p = pattern(shape, **kw)
    return m.analytic(p, fits_from_p(p, np.array(shape[:4])), L)


def iss_c(**kw):
    """Falling crossing of 95% recovery, and Iss there."""
    grid = np.geomspace(0.02, 60.0, 15)
    vals = np.array([recovery(tuple(v * s for v in FARRIS), L, **kw)
                     for s in grid])
    j = int(np.nanargmax(vals))
    if vals[j] < 0.95:
        return float("nan"), float("nan"), float(vals[j])
    for i in range(j, len(grid) - 1):
        if not (vals[i] >= 0.95 > vals[i + 1]):
            continue
        lo, hi = grid[i], grid[i + 1]
        for _ in range(22):
            mid = 0.5 * (lo + hi)
            if recovery(tuple(v * mid for v in FARRIS), L, **kw) >= 0.95:
                lo = mid
            else:
                hi = mid
        s = 0.5 * (lo + hi)
        return s, s34.e_iss_vector(pattern(tuple(v * s for v in FARRIS), **kw)), \
            float(vals.max())
    return float("nan"), float("nan"), float(vals.max())


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"Farris quartet {FARRIS}, L = {L}\n")
    print(f"  {'model':>16} {'critical TL':>12} {'Iss.c':>8} {'max rec':>8}")
    rows = []

    def go(label, **kw):
        s, ic, mx = iss_c(**kw)
        tl = sum(FARRIS) * s if np.isfinite(s) else float("nan")
        rows.append(dict(model=label, scale=s, tree_length=tl, iss_c=ic,
                         max_rec=mx, **kw))
        print(f"  {label:>16} {tl:12.3f} {ic:8.4f} {mx:8.3f}", flush=True)

    for alpha in (0.5, 1.0, 2.0, 4.0, 16.0, np.inf):
        go(f"G a={alpha:g} k=8", alpha=alpha, k=8)
    for pinv in (0.0, 0.2, 0.4, 0.6):
        go(f"+I P={pinv:g}", pinv=pinv)
    # the two-category Gamma, which at small alpha is a +I model
    for alpha in (0.05, 0.1, 0.2, 0.5, 1.0, 2.0, np.inf):
        go(f"G a={alpha:g} k=2", alpha=alpha, k=2)

    d = pd.DataFrame(rows)
    d.to_csv(OUT / "gamma_pinv.csv", index=False)

    print("\n  the two-category case, with the implied size of its slow class:")
    for _, x in d[d.model.str.contains("k=2")].iterrows():
        a = x.get("alpha", np.nan)
        if "alpha" not in x or not np.isfinite(a):
            continue
        r = discrete_gamma(float(a), 2)
        print(f"    a={float(a):<5g}  rates = {r[0]:.4f}, {r[1]:.4f}"
              f"   (mass 0.5 on each)")

    g8 = d[d.model.str.contains("k=8")].dropna(subset=["iss_c"])
    g2 = d[d.model.str.contains("k=2")].dropna(subset=["iss_c"])
    ii = d[d.model.str.startswith("+I")].dropna(subset=["iss_c"])
    if len(g8) > 1:
        print(f"\n  Iss.c over the Gamma shapes (k=8): "
              f"{g8.iss_c.min():.4f} to {g8.iss_c.max():.4f}"
              f"   spread {g8.iss_c.max() - g8.iss_c.min():.4f}")
    if len(g2) > 1:
        print(f"  Iss.c over the Gamma shapes (k=2): "
              f"{g2.iss_c.min():.4f} to {g2.iss_c.max():.4f}"
              f"   spread {g2.iss_c.max() - g2.iss_c.min():.4f}")
    if len(ii) > 1:
        print(f"  Iss.c over the +I fractions:      "
              f"{ii.iss_c.min():.4f} to {ii.iss_c.max():.4f}"
              f"   spread {ii.iss_c.max() - ii.iss_c.min():.4f}")


if __name__ == "__main__":
    main()
