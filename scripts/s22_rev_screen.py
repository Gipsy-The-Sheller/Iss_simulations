"""How far the JC pair screen generalises, model by model.

The screen has four moving parts and they do not generalise to the same
extent.  Separating them is the point of this script.

  (a) p-hat is binomial.  For ANY model: given a pair at path length d the
      indicator that a site differs is Bernoulli with probability
      wp(d) = 1 - sum_j pi_j P_jj(d), and sites evolve independently, so the
      count of differing sites over L sites is exactly Binomial(L, wp(d)).
      No reversibility, no exchangeabilities and no state count enters.  The
      binomial half of the threshold therefore generalises verbatim.

  (b) The boundary.  wp(d) -> 1 - sum_j pi_j^2 as d -> infinity, which is 3/4
      only at uniform composition.  So the threshold is 1 - sum pi^2 at the
      alignment's own composition, for every reversible model including GTR.

  (c) The corrected distance.  d = -(3/4) ln(1 - 4p/3) is the JC case of a
      log-of-linear inverse, available exactly when the pairwise process has
      ONE distinct non-zero eigenvalue.

  (d) What p alone can recover.  For a fixed model d -> wp(d) is monotone and
      inverts; the question is whether p carries enough to do it, which
      depends on how many summaries the pairwise divergence has.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench.models import build_Q                                   # noqa: E402
from issbench.spectral import factorize                               # noqa: E402

OUT = ROOT / "results" / "rev_screen"
UNIF = np.ones(4) / 4
SKEW = np.array([0.40, 0.30, 0.20, 0.10])
R_GTR = np.array([2.1, 4.7, 0.9, 1.3, 3.6, 1.0])

#            name              build_Q model   kwargs                          pi
MODELS = [
    ("JC69",            "JC69",  dict(),                                  UNIF),
    ("K80 k=5",         "K80",   dict(kappa=5.0),                         UNIF),
    ("F81 skewed",      "F81",   dict(),                                  SKEW),
    ("F84 k=0.5",       "F84",   dict(kappa=0.5),                         SKEW),
    ("HKY85 k=5",       "HKY85", dict(kappa=5.0),                         SKEW),
    ("TN93 k=3,8",      "TN93",  dict(kappa=3.0, kappa2=8.0),             SKEW),
    ("GTR skewed",      "GTR",   dict(R=R_GTR),                           SKEW),
    ("GTR uniform",     "GTR",   dict(R=R_GTR),                           UNIF),
]
# free pairwise quantities beyond the branch length
SUMMARIES = {"JC69": (0, "p"),
             "K80 k=5": (1, "p and ti/tv"),
             "F81 skewed": (0, "p"),
             "F84 k=0.5": (1, "p and ti/tv"),
             "HKY85 k=5": (1, "p and ti/tv"),
             "TN93 k=3,8": (2, "p and the two transitions"),
             "GTR skewed": (5, "the whole 4x4 divergence matrix"),
             "GTR uniform": (5, "the whole 4x4 divergence matrix")}


def threshold(boundary: float, L: int = 1000, z: float = 1.6449) -> float:
    """Largest true p whose observed value stays below the boundary."""
    lo, hi = 0.0, boundary
    for _ in range(200):
        mid = (lo + hi) / 2
        if boundary - mid > z * np.sqrt(mid * (1 - mid) / L):
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------ (a) the binomial half
    print("(a) is the differing-site count binomial for every model?")
    print("    the per-site probability of differing is 1 - sum_j pi_j P_jj(d);")
    print("    this checks it against the entrywise sum over P(d)")
    worst = 0.0
    for name, m, kw, pi in MODELS:
        Q = build_Q(m, pi=pi, **kw)
        sp = factorize(Q, pi)
        for d in (0.05, 0.5, 2.0, 8.0):
            direct = 1.0 - float(np.sum(pi * np.diag(expm(Q * d))))
            worst = max(worst, abs(direct - sp.w_p(d)))
    print(f"    max discrepancy over the whole family and four d values: "
          f"{worst:.2e}")
    print("    it is an identity, not a coincidence: sum_{x != y} pi_x P_xy(d)")
    print("    = 1 - sum_j pi_j P_jj(d).  Sites being independent, the count")
    print("    is therefore Binomial(L, wp) for every model, so the binomial")
    print("    half of the JC threshold carries over with nothing changed.")

    # ------------------------------------------- (b,c) boundary and inverse
    print("\n(b,c) boundary and corrected distance, by model")
    print(f"    {'model':14s} {'#distinct':>9} {'boundary':>9} "
          f"{'wp(d=1)':>8} {'wp(d=10)':>9}  inverse")
    rows = []
    for name, m, kw, pi in MODELS:
        Q = build_Q(m, pi=pi, **kw)
        sp = factorize(Q, pi)
        n_distinct = len(np.unique(np.round(sp.lam[1:], 10)))
        boundary = 1.0 - float(np.sum(pi ** 2))
        inv = ("closed form, one eigenvalue" if n_distinct == 1
               else "1-D monotone, no closed form")
        rows.append(dict(model=name, n_distinct=n_distinct, boundary=boundary,
                         wp_d1=float(sp.w_p(1.0)), wp_d10=float(sp.w_p(10.0)),
                         inverse=inv, threshold=threshold(boundary)))
        print(f"    {name:14s} {n_distinct:9d} {boundary:9.4f} "
              f"{sp.w_p(1.0):8.4f} {sp.w_p(10.0):9.4f}  {inv}")
    pd.DataFrame(rows).to_csv(OUT / "model_family.csv", index=False)
    print("    the boundary is 1 - sum pi^2 at each model's own composition;")
    print("    only the two uniform rows have 3/4.")

    # ------------------------------------------- (d) what the correction needs
    print("\n(d) the corrected distance is closed form where the pairwise")
    print("    process splits into independent two-state processes, and")
    print("    numerical otherwise.  Which summaries it must consume:")
    print(f"    {'model':14s} {'#distinct':>9}  {'correction consumes':<38}")
    consumes = {
        "JC69": "p",
        "K80 k=5": "p and the transversion proportion (Kimura)",
        "F81 skewed": "p",
        "F84 k=0.5": "numerical: three rates",
        "HKY85 k=5": "numerical, or the HKY closed form given kappa",
        "TN93 k=3,8": "the two transition proportions (Tamura-Nei)",
        "GTR skewed": "numerical: the whole divergence matrix",
        "GTR uniform": "numerical: the whole divergence matrix"}
    for name, m, kw, pi in MODELS:
        Q = build_Q(m, pi=pi, **kw)
        sp = factorize(Q, pi)
        n = len(np.unique(np.round(sp.lam[1:], 10)))
        print(f"    {name:14s} {n:9d}  {consumes[name]:<38}")
    print("    the split is not JC against GTR: it is one eigenvalue")
    print("    (JC69, F81) against the rest, with K80 and TN93 rescued by")
    print("    their decomposition and F84, HKY85 and GTR not.")

    # ---------------------------------------- the threshold, restated once
    print("\n    the threshold in one line per model: replace 3/4 by the")
    print("    model's own boundary, and nothing else changes (L = 1000, 95%)")
    for name, m, kw, pi in MODELS:
        b = 1.0 - float(np.sum(pi ** 2))
        print(f"      {name:14s} {b:.4f} -> {threshold(b):.4f}")
    print("    so the SCREEN generalises to every REV model, GTR included:")
    print("    (a) is an identity and (b) is a limit, and neither asks for")
    print("    anything a submodel has and GTR lacks.")


if __name__ == "__main__":
    main()
