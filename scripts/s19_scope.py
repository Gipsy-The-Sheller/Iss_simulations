"""Two scope questions raised by the preprint review, both answerable.

(1) THE AFFINITY TEST OUTSIDE A STAR.  Theorem 13 was verified in s15 on an
equidistant star, where every pair has the same path length and the pair
average is trivial.  Nothing in the theorem asks for that: Psi and wp are pair
averages of tr P(d) and of the single-site identity, so any topology just
replaces the single d by the pair-distance multiset.  This part runs the test
on an unbalanced quartet (internal fractions 0.10 and 0.50), on a balanced one
at N = 8 and N = 16, and under GTR, to see whether the criterion's verdict ---
and the size of the failure --- depend on the tree.

(2) DOES THE CRITERION SURVIVE THE PARSIMONY-INFORMATIVE RESTRICTION?
Duchene et al. (2022) recommend restricting to parsimony-informative sites,
and the manuscript repeats the recommendation.  The review asks whether
Theorem 13 still holds there.  The tempting argument --- the weights h(c)
change but xi_mu does not --- does not settle it: conditioning on a column
being informative changes the column law, so E[Hbar] is no longer an i.i.d.
mean of the unconditioned W and Lemma (variance aggregation) does not apply as
written.  What can be checked directly is whether the result is still affine
in the mean p-distance of the SAME conditioned ensemble, which is what this
part measures by enumeration.

    python scripts/s19_scope.py
"""
from __future__ import annotations

import sys
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench import trees                                          # noqa: E402
from issbench.entropy import h                                      # noqa: E402
from issbench.models import build_Q                                 # noqa: E402
from issbench.spectral import factorize                             # noqa: E402

OUT = ROOT / "results" / "scope"
UNIF = np.ones(4) / 4
SKEW = np.array([0.40, 0.30, 0.20, 0.10])
R_GTR = np.array([2.1, 4.7, 0.9, 1.3, 3.6, 1.0])


def scalars(sp, dmat, rates=(1.0,), weights=(1.0,)):
    """Psi and the mean p-distance for any tree, from its pair distances."""
    iu = np.triu_indices(dmat.shape[0], 1)
    d = dmat[iu]
    xi = np.zeros_like(sp.lam)
    for r, w in zip(rates, weights):
        xi += w * np.mean(np.exp(np.outer(r * d, sp.lam)), axis=0)
    Psi = float(np.dot(sp.tau[1:], xi[1:]))
    wp = float((1.0 - np.sum(sp.pi ** 2)) - np.dot(sp.w[1:], xi[1:]))
    return Psi, wp


def balanced(n: int, a: float, r: float, labels=None):
    """A balanced binary tree of n tips with pendant a and internal r."""
    labels = labels or [f"T{i+1}" for i in range(n)]
    if n == 1:
        return trees.Node(labels[0], a)
    half = n // 2
    return trees.Node("", r / 2.0,
                      [balanced(half, a, r, labels[:half]),
                       balanced(n - half, a, r, labels[half:])])


def affinity(sp, node, rates=(1.0,), weights=(1.0,), n_pt=120):
    """Residual of Psi against wp along a scaling family of the whole tree."""
    def scale(v, s):
        if v.is_leaf:
            return trees.Node(v.name, v.length * s)
        return trees.Node(v.name, v.length * s,
                          [scale(c, s) for c in v.children])

    s = np.exp(np.linspace(np.log(0.02), np.log(1.5), n_pt))
    psi, wp = [], []
    for si in s:
        P, Q = scalars(sp, trees.pairwise_distances(scale(node, si)),
                       rates, weights)
        psi.append(P)
        wp.append(Q)
    psi, wp = np.array(psi), np.array(wp)
    b, a = np.polyfit(wp, psi, 1)
    r = psi - (a + b * wp)
    return dict(slope=float(b), max_abs=float(np.abs(r).max()),
                span=float(psi.max() - psi.min()),
                pct=float(100 * np.abs(r).max() / (psi.max() - psi.min())))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------- (1) topology scope
    print("(1) Psi vs p-bar affinity, outside the equidistant star")
    print(f"    {'model':16s} {'topology':34s} {'tau/w':26s} "
          f"{'slope':>8}  {'max|resid|':>11}  verdict")
    rows = []
    cases = [
        ("HKY k=5, skewed", build_Q("HKY85", kappa=5.0, pi=SKEW), SKEW),
        ("GTR, skewed", build_Q("GTR", pi=SKEW, R=R_GTR), SKEW),
        ("GTR, uniform", build_Q("GTR", pi=UNIF, R=R_GTR), UNIF),
    ]
    tops = [("star, N=4", trees.star(4, 1.0)),
            ("quartet fi=0.35", trees.quartet(1.0, 1.0, 1.0, 1.0, 1.08)),
            ("quartet fi=0.10", trees.quartet(1.0, 1.0, 1.0, 1.0, 0.22)),
            ("quartet fi=0.50", trees.quartet(1.0, 1.0, 1.0, 1.0, 2.0)),
            ("balanced, N=8", balanced(8, 1.0, 0.6)),
            ("balanced, N=16", balanced(16, 1.0, 0.6))]
    for mname, Q, pi in cases:
        sp = factorize(Q, pi)
        ratios = f"{[round(float(x), 3) for x in sp.ratios]}"
        for tname, node in tops:
            a = affinity(sp, node)
            holds = sp.affine_ok()
            rows.append(dict(model=mname, topology=tname, affine=holds, **a))
            print(f"    {mname:16s} {tname:34s} {ratios:26s} "
                  f"{a['slope']:8.3f}  {a['max_abs']:11.2e}  "
                  f"{'holds' if holds else 'FAILS'}")
    pd.DataFrame(rows).to_csv(OUT / "affinity_scope.csv", index=False)

    # the failure is a property of the model, so one topology suffices to
    # report its size; print it once for the load-bearing case
    bad = [r for r in rows if not r["affine"]]
    if bad:
        print(f"\n    failure size across the six topologies, "
              f"max|resid| in Psi: "
              f"{min(r['max_abs'] for r in bad):.2e} to "
              f"{max(r['max_abs'] for r in bad):.2e}")

    # ------------------------------------- (2) informative-site restriction
    # The comparison has to be two-sided.  Under conditioning the columns stop
    # being i.i.d. draws from the model, so Lemma (variance aggregation) does
    # not apply and there is no second-order reason for E[Hbar] to be a
    # function of p-bar alone; the test is therefore whether the restriction
    # ADDS a departure, measured against the unrestricted residual for the
    # same model and the same tree.  A model whose criterion holds and one
    # whose criterion fails are both run, so the two effects are separated.
    print("\n(2) E[Hbar] against the mean p-distance, unrestricted and "
          "restricted to\n    parsimony-informative columns (N=6 star, "
          "L -> infinity)")
    print(f"    {'model':22s} {'tau/w':22s} {'condition':12s} "
          f"{'slope':>8}  {'max|resid|':>11}  {'%of range':>9}")
    N = 6
    comps = _comps(N, 4)
    w_all = np.array([h(np.array(c, dtype=float) / N) for c in comps])
    keep = np.array([sum(1 for x in c if x >= 2) >= 2 for c in comps])
    rows2 = []
    for mname, Q, pi in (("GTR, uniform (holds)", build_Q("GTR", pi=UNIF, R=R_GTR), UNIF),
                         ("HKY k=5, skewed (fails)", build_Q("HKY85", kappa=5.0, pi=SKEW), SKEW)):
        sp = factorize(Q, pi)

        grid = np.exp(np.linspace(np.log(0.05), np.log(4.0), 300))
        for what in ("all columns", "informative"):
            pts = [t for x in grid
                   for t in curve_for(sp, pi, float(x), comps, w_all, keep, N)
                   if t[2] == what]
            pts = [(a, b) for a, b, _ in pts if np.isfinite(a)]
            y = np.array([a for a, _ in pts])
            x = np.array([b for _, b in pts])
            b, a0 = np.polyfit(x, y, 1)
            res = y - (a0 + b * x)
            pct = 100 * np.abs(res).max() / (y.max() - y.min())
            rows2.append(dict(model=mname, condition=what, slope=float(b),
                              max_abs=float(np.abs(res).max()), pct=float(pct)))
            print(f"    {mname:22s} {str([round(float(r), 3) for r in sp.ratios]):22s} "
                  f"{what:12s} {b:8.3f}  {np.abs(res).max():11.2e}  {pct:8.3f}%")
    pd.DataFrame(rows2).to_csv(OUT / "restriction.csv", index=False)
    print("      the residuals there are dominated by the O(3) remainder that\n"
          "      Figure 1B already shows masks the criterion at the E level,\n"
          "      so that table cannot separate the two effects.")

    # ------------------------------- (2b) two models matched on p-distance
    # The criterion's content is that E[Iss] is fixed by the mean p-distance,
    # so two models that BOTH satisfy it must agree at matched p-bar.  This is
    # the test the restriction can fail cleanly: JC69 and GTR with uniform
    # composition both have tau/w = K for every non-zero eigenvalue, so any
    # disagreement between them at equal p-bar is either O(3) --- the
    # unrestricted case, already characterised --- or the restriction.
    print("\n    matched p-bar: JC69 against GTR with uniform composition,\n"
          "    both of which satisfy the criterion exactly")
    grid = np.exp(np.linspace(np.log(0.05), np.log(4.0), 400))
    traces = {}
    for mname, Q in (("JC69", build_Q("JC69")),
                     ("GTR uni", build_Q("GTR", pi=UNIF, R=R_GTR))):
        sp = factorize(Q, UNIF)
        for what in ("all columns", "informative"):
            x, y = [], []
            for d in grid:
                for a, b, tag in curve_for(sp, UNIF, float(d), comps,
                                           w_all, keep, N):
                    if tag == what and np.isfinite(a):
                        x.append(b)
                        y.append(a)
            traces[(mname, what)] = (np.array(x), np.array(y))
    for what in ("all columns", "informative"):
        xa, ya = traces[("JC69", what)]
        xb, yb = traces[("GTR uni", what)]
        common = np.linspace(max(xa.min(), xb.min()),
                             min(xa.max(), xb.max()), 200)
        ja = np.interp(common, xa, ya)
        gb = np.interp(common, xb, yb)
        gap = np.abs(ja - gb)
        rows2.append(dict(model="JC69 vs GTR-uniform", condition=what,
                          slope=np.nan, max_abs=float(gap.max()),
                          pct=float(100 * gap.max() / (ja.max() - ja.min()))))
        print(f"    {'JC69 vs GTR-uniform':22s} {'matched':22s} {what:12s} "
              f"{'':8s}  {gap.max():11.2e}  "
              f"{100 * gap.max() / (ja.max() - ja.min()):8.3f}%")
    pd.DataFrame(rows2).to_csv(OUT / "restriction.csv", index=False)


def curve_for(sp, pi, d, comps, w_all, keep, N):
    """(E[Hbar], p-bar, tag) for all columns and for informative ones."""
    P = sp.P(d)
    g = np.zeros(len(comps))
    for rho in range(4):
        q = P[rho]
        for i, c in enumerate(comps):
            cls = (comb(N, c[0]) * comb(N - c[0], c[1])
                   * comb(N - c[0] - c[1], c[2]))
            g[i] += pi[rho] * cls * np.prod(q ** np.array(c))
    pbar = lambda z: sum(z[i] * (1.0 - sum(x * (x - 1) for x in c)
                                 / (N * (N - 1)))
                         for i, c in enumerate(comps))
    out = [(float(w_all @ g), float(pbar(g)), "all columns")]
    z = g[keep].sum()
    if z > 1e-12:
        gc = g * keep / z
        out.append((float(w_all @ gc), float(pbar(gc)), "informative"))
    else:
        out.append((np.nan, np.nan, "informative"))
    return out


def _comps(N, K, prefix=()):
    if K == 1:
        return [prefix + (N,)]
    out = []
    for i in range(N + 1):
        out += _comps(N - i, K - 1, prefix + (i,))
    return out


if __name__ == "__main__":
    main()
