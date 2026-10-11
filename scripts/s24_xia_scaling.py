"""Xia's own topologies, scaled: is Iss affine in the four-point margin?

The refinement: do not vary the model, the topology or the tree's profile --
vary only the overall tree length, which is how the Iss.c table was generated.
Under that restriction both quantities are functions of one variable, so the
question is whether Iss is affine in the additivity margin M, and if so what
Iss is at M = 0.

Two trees, taken from Xia's symmetrical and asymmetrical families.  Both are
written in the paper's rooted notation and unrooted here, which matters: the
root is suppressed and its two branches merge into the internal edge, so the
balanced quartet ((A:.1,B:.1):.1,(C:.1,D:.1):.1) has internal length 0.2, not
0.1.

    balanced      ((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1)   -> a=b=c=d=0.1, e=0.2
    asymmetrical  (((A:0.1,B:0.1):0.1,C:0.2):0.1,D:0.3)   -> a=b=0.1, c=0.2,
                                                             d=0.4, e=0.1

Iss is the exact E[Iss] by pattern enumeration, not the second-order form,
because the second-order form is invalid at the small-divergence end that
this sweep has to cover.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench import exact, trees                                    # noqa: E402
from issbench.models import build_Q                                   # noqa: E402

OUT = ROOT / "results" / "xia_scaling"
Q = build_Q("JC69")
PI = np.ones(4) / 4


def margin(a, b, c, d, e):
    A, B, C, D, E = (np.exp(-4 * x / 3) for x in (a, b, c, d, e))
    return 0.75 * (A * B + C * D - E * (A * C + B * D))


def iss(a, b, c, d, e):
    return exact.E_Iss(4, trees.quartet(a, b, c, d, e), Q, PI, convention="A")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    TREES = {"balanced (a=b=c=d=0.1, e=0.2)": (0.1, 0.1, 0.1, 0.1, 0.2),
             "asymmetrical (a=b=0.1, c=0.2, d=0.4, e=0.1)":
                 (0.1, 0.1, 0.2, 0.4, 0.1)}

    s = np.exp(np.linspace(np.log(0.02), np.log(8.0), 220))
    print(f"{'s':>7} {'M':>11} {'Iss':>9}   (both trees)")
    rows = []
    for name, base in TREES.items():
        for si in s:
            v = [x * si for x in base]
            rows.append(dict(tree=name, s=si, M=margin(*v), iss=iss(*v)))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "scaling.csv", index=False)

    for name in TREES:
        g = d[d.tree == name].sort_values("s")
        print(f"\n{name}")
        print(f"  M over the sweep: {g.M.min():+.3e} to {g.M.max():+.3e}"
              f"   sign changes: {int((np.diff(np.sign(g.M)) != 0).sum())}")
        # is Iss affine in M?
        x, y = g.M.to_numpy(), g.iss.to_numpy()
        b, a0 = np.polyfit(x, y, 1)
        r = y - (a0 + b * x)
        print(f"  Iss against M: slope {b:+.4f}  intercept {a0:+.4f}  "
              f"max |resid| {np.abs(r).max():.2e} "
              f"({100 * np.abs(r).max() / (y.max() - y.min()):.2f}% of range)")
        # the same against the mean p-distance, for comparison
        pp = []
        for si in g.s:
            v = [q * si for q in TREES[name]]
            a_, b_, c_, d_, e_ = v
            D = [a_ + b_, c_ + d_, a_ + c_ + e_, b_ + d_ + e_,
                 a_ + d_ + e_, b_ + c_ + e_]
            pp.append(np.mean([0.75 * (1 - np.exp(-4 * z / 3)) for z in D]))
        pp = np.array(pp)
        b2, a2 = np.polyfit(pp, y, 1)
        r2 = y - (a2 + b2 * pp)
        print(f"  Iss against mean p: slope {b2:+.4f}  max |resid| "
              f"{np.abs(r2).max():.2e} "
              f"({100 * np.abs(r2).max() / (y.max() - y.min()):.2f}% of range)")

    # ------------------------------------------------------- the s -> 0 end
    print("\nwhat the two quantities do at small s (asymmetrical tree)")
    base = TREES["asymmetrical (a=b=0.1, c=0.2, d=0.4, e=0.1)"]
    print(f"    {'s':>7} {'M':>12} {'M/s^2':>10} {'Iss':>10}")
    for si in (0.02, 0.05, 0.1, 0.2, 0.5):
        v = [x * si for x in base]
        m = margin(*v)
        print(f"    {si:7.3f} {m:12.4e} {m / si ** 2:10.4f} {iss(*v):10.5f}")
    print("    M goes as s^2 at the origin: the linear terms cancel because")
    print("    both trees sit in the balanced combination.  So M has a")
    print("    double zero at s = 0 and no crossing at positive s.")

    # ----------------------------------------- so what IS Xia's threshold?
    # If NJ never fails systematiically on these topologies, the tabulated
    # critical value has to be a detectability threshold: the tree length at
    # which the margin stops standing out from sampling noise.  On the
    # p-distance scale the observed margin is p13 + p24 - p12 - p34, whose
    # variance is (1/L) times a sum of p(1-p) over the four terms plus
    # covariances; the diagonal part is computed here and the covariances
    # are dropped, so the number is a guide rather than an identity.
    print("\nif not M = 0, then M = noise: the length at which the margin")
    print("stops being detectable, and the Iss there, for both trees")
    print(f"    {'tree':<12} {'L':>6} {'end':<6} {'s*':>9} {'M*':>9} "
          f"{'Iss at s*':>10}")
    rows2 = []
    for name, base in TREES.items():
        for L in (350, 1000, 4000, 20000):
            def moments(si):
                """(E[M], SD[M]) of the OBSERVED margin, exactly.

                M-hat is a mean over sites of a per-site variable, so its
                variance is Var(m)/L with m taken under the pattern law --
                no diagonal approximation.  The four pairwise indicators
                share taxa and are correlated; dropping those covariances
                shifts the threshold by a lot.
                """
                from issbench.exact import tree_patterns
                v = [x * si for x in base]
                a_, b_, c_, d_, e_ = v
                D = {"13": a_ + c_ + e_, "24": b_ + d_ + e_,
                     "12": a_ + b_, "34": c_ + d_}
                m = (0.75 * (1 - np.exp(-4 * D["13"] / 3))
                     + 0.75 * (1 - np.exp(-4 * D["24"] / 3))
                     - 0.75 * (1 - np.exp(-4 * D["12"] / 3))
                     - 0.75 * (1 - np.exp(-4 * D["34"] / 3)))
                pats = tree_patterns(trees.quartet(a_, b_, c_, d_, e_),
                                     Q, PI, (1.0,), (1.0,))
                vals, probs = [], []
                for pat, pr in pats.items():
                    x1, x2, x3, x4 = pat
                    vals.append(float(x1 != x3) + float(x2 != x4)
                                - float(x1 != x2) - float(x3 != x4))
                    probs.append(pr)
                vals, probs = np.array(vals), np.array(probs)
                var = float((vals ** 2) @ probs) - m ** 2
                return m, np.sqrt(max(var, 0.0) / L)

            def gap(si):
                m, sd = moments(si)
                return m - sd
            # M rises to a maximum and then decays to zero as everything
            # saturates, while the sampling noise of the margin tends to a
            # constant.  So the margin is detectable on a WINDOW, and the
            # first crossing is the one that matters.
            grid = np.exp(np.linspace(np.log(1e-4), np.log(60.0), 4000))
            vals = np.array([gap(x) for x in grid])
            cross = np.where(np.diff(np.sign(vals)) != 0)[0]
            if not len(cross):
                print(f"    {name.split()[0]:<12} {L:6d}   never detectable")
                continue
            got = []
            for ci in cross:
                lo, hi = grid[ci], grid[ci + 1]
                for _ in range(200):
                    mid = (lo + hi) / 2
                    if gap(mid) * gap(lo) <= 0:
                        hi = mid
                    else:
                        lo = mid
                got.append((lo + hi) / 2)
            for tag, sstar in zip(("lower", "upper"), got):
                v = [x * sstar for x in base]
                m = margin(*v)
                I = iss(*v)
                rows2.append(dict(tree=name, L=L, end=tag, s_star=sstar,
                                  M=m, iss=I))
                print(f"    {name.split()[0]:<12} {L:6d} {tag:<6} "
                      f"{sstar:9.4f} {m:9.5f} {I:10.5f}")
    pd.DataFrame(rows2).to_csv(OUT / "detectability.csv", index=False)
    g = pd.DataFrame(rows2)
    if len(g):
        for name in TREES:
            sub = g[g.tree == name]
            if len(sub) > 1:
                for end in ("lower", "upper"):
                    e2 = sub[sub.end == end]
                    if len(e2) > 1:
                        print(f"    {name.split()[0]:<12} {end:<6} Iss runs "
                              f"{e2.iss.min():.4f} to {e2.iss.max():.4f} "
                              f"over L = {e2.L.min()} to {e2.L.max()}")


if __name__ == "__main__":
    main()
