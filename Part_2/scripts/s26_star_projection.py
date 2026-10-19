"""Is the KL projection of a quartet onto a wrong topology always a star?

Claim, found numerically and now to be tested over a wide domain: if P lies on
the quartet topology (A,B)|(C,D) and Q is its KL projection onto the
alternative (A,C)|(B,D), then Q has zero internal branch.  The consequence
would be that the two wrong topologies fit identically -- at zero internal
branch there is no split and the star belongs to all three manifolds -- so
the three-way ML choice on a quartet would really be two-way, and the error
probability a one-dimensional normal tail rather than a bivariate orthant
probability.

That is a strong claim and it needs a search for counterexamples before any
proof.  A grid of pendant lengths crossed with internal branch lengths is
scanned and every projection's internal branch recorded.

NOTE ON THE FIRST RUN OF THIS SCAN.  It flagged 11 of 230 projections as
interior and the result was committed as "the star-projection claim is false
in general".  That was wrong, and the reason is worth keeping: the flag was
a threshold on the fitted branch length, `internal_t > 0.02`, which is a
statement about the optimiser's stopping point and not about the geometry.
Near saturation the objective is flat along the internal branch -- once the
internal edge has washed out, extending it further changes the pattern
distribution by nothing -- so an optimiser that starts away from the star
stops at whatever branch length it happens to reach, with a KL exactly equal
to the star's.  Checking two of the eleven by hand gave KL differences of
1e-14.  The verdict printed below therefore compares objectives, and s30
repeats the whole scan with a properly constrained star fit over a domain
that extends to e_true = 12, where a failure would be most likely.

The earlier version of this scan was killed by its own cost: it called the
pattern enumeration inside the optimiser.  P(x) is multilinear in the five
x_e = exp(-4 t_e / 3), so it is written out directly here, in O(K) per
pattern, and checked against tree_patterns before it is used.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "src"))

from issbench import trees                                          # noqa: E402
from issbench.exact import tree_patterns                           # noqa: E402
from issbench.models import build_Q                                 # noqa: E402

OUT = ROOT / "results" / "star_projection"
Q = build_Q("JC69")
PI = np.ones(4) / 4
SPLITS = {"12|34": ("A", "B", "C", "D"),
          "13|24": ("A", "C", "B", "D"),
          "14|23": ("A", "D", "B", "C")}


def to_x(t):
    return np.exp(-4.0 * np.asarray(t, float) / 3.0)


def to_t(x):
    return -0.75 * np.log(np.clip(x, 1e-15, 1.0))


def pattern_fast(lengths, labels):
    """P over the 256 patterns, packed in the fixed leaf order A,B,C,D.

    Rooting at the node carrying the first two leaves of `labels`, the
    pattern probability is

        P(x) = sum_{s,s'} pi_s L_{ab}(s) P_{ss'}(e) R_{cd}(s')

    with P_ij(t) = (1 - x_t)/4 + x_t * 1[i = j].  Writing P_{ss'}(e) =
    alpha + beta * 1[s = s'] turns the double sum into two single ones, so
    each pattern costs O(K) and no tree is built.  `labels` only decides
    which pendant length and which leaf state sit in which of the four
    slots, which is exactly how the topology is meant to change.
    """
    xa, xb, xc, xd, xe = to_x(lengths)
    alpha = (1.0 - xe) / 4.0
    beta = xe
    order = [labels.index(n) for n in ("A", "B", "C", "D")]
    states = np.arange(4)

    def pair(xi, xj):
        """(16, 4): entry (i*4+j, s) = P_{s,i}(xi) P_{s,j}(xj)."""
        out = np.zeros((16, 4))
        for i in range(4):
            for j in range(4):
                out[i * 4 + j] = (((1 - xi) / 4 + xi * (states == i))
                                  * ((1 - xj) / 4 + xj * (states == j)))
        return out

    pab, pcd = pair(xa, xb), pair(xc, xd)
    ssab = pab @ PI                              # sum_s pi_s L(s)   (pi at u)
    sscd = pcd.sum(axis=1)                       # sum_s' R(s')       (no pi)
    cross = (pab * PI) @ pcd.T                   # (16,16): sum_s pi_s L R
    out = np.zeros(256)
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for l in range(4):
                    idx = 0
                    for p in order:
                        idx = (idx << 2) | (i, j, k, l)[p]
                    out[idx] = (alpha * ssab[i * 4 + j] * sscd[k * 4 + l]
                                + beta * cross[i * 4 + j, k * 4 + l])
    return out


def pattern_tree(lengths, labels):
    """Reference implementation through tree_patterns, for the check."""
    node = trees.quartet(*lengths, labels=labels)
    pats = tree_patterns(node, Q, PI, (1.0,), (1.0,))
    order = [labels.index(n) for n in ("A", "B", "C", "D")]
    v = np.zeros(256)
    for pat, pr in pats.items():
        idx = 0
        for p in order:
            idx = (idx << 2) | pat[p]
        v[idx] = pr
    return v


def kl(p, q, floor=1e-300) -> float:
    m = p > 0
    return float(np.sum(p[m] * np.log(p[m] / np.maximum(q[m], floor))))


def project(p_true, labels, x0):
    """(KL, internal branch) of the best fit, in two parametrisations.

    The second reparametrises t = -(3/4) ln x with x in (0,1), so the
    zero-length boundary is the interior point x = 1 rather than a bound the
    optimiser can stick to.
    """
    def obj_t(th):
        if np.any(th <= 1e-10) or np.any(th > 25):
            return 1e3
        return kl(p_true, pattern_fast(th, labels))

    best = (np.inf, np.asarray(x0, float))
    for trial in (x0, list(np.asarray(x0) * 0.3), [0.4, 0.4, 0.4, 0.4, 0.05]):
        r = minimize(obj_t, trial, method="L-BFGS-B", bounds=[(1e-10, 25)] * 5)
        if r.fun < best[0]:
            best = (float(r.fun), r.x)

    def obj_x(x):
        return kl(p_true, pattern_fast(to_t(x), labels))

    xs = to_x(best[1])
    bx = (np.inf, xs)
    for trial in (xs, np.full(5, 0.6)):
        r = minimize(obj_x, trial, method="L-BFGS-B",
                     bounds=[(1e-9, 1.0)] * 5)
        if r.fun < bx[0]:
            bx = (float(r.fun), r.x)
    return best[0], float(best[1][4]), bx[0], float(to_t(bx[1])[4])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    print("checking the fast pattern formula against tree_patterns")
    worst = 0.0
    for b1 in ((0.7, 0.2, 0.5, 0.15, 0.3), (1.0, 0.05, 1.0, 0.05, 0.8),
               (0.3, 0.3, 0.3, 0.3, 0.02)):
        for lab in SPLITS.values():
            worst = max(worst, np.abs(pattern_fast(b1, lab)
                                      - pattern_tree(b1, lab)).max())
    print(f"    max |fast - tree_patterns| over 9 cases: {worst:.2e}")
    if worst > 1e-12:
        sys.exit("the fast formula disagrees; fix it before scanning")

    rng = np.random.default_rng(20261014)
    cases = [(*rng.uniform(0.02, 1.5, 4), float(rng.uniform(1e-4, 3.0)))
             for _ in range(80)]
    for e in (1e-4, 0.01, 0.1, 0.5, 1.0, 2.0, 4.0):
        for pend in ((1.0, 0.05, 1.0, 0.05), (1.0, 0.05, 0.7, 0.2),
                     (0.8, 0.3, 0.5, 0.1), (1.5, 1.2, 0.08, 0.05),
                     (0.5, 0.5, 0.5, 0.5)):
            cases.append((*pend, e))

    print(f"\nscanning {len(cases)} configurations x 2 alternatives")
    rows = []
    for n, b1 in enumerate(cases):
        p = pattern_fast(b1, SPLITS["12|34"])
        for lab in ("13|24", "14|23"):
            kt, it, kx, ix = project(p, SPLITS[lab], b1)
            rows.append(dict(case=n, topology=lab, e_true=b1[4],
                             pend=",".join(f"{v:.3f}" for v in b1[:4]),
                             kl=kt, internal_t=it, kl_x=kx, internal_x=ix))
        if (n + 1) % 20 == 0:
            print(f"    {n + 1}/{len(cases)}", flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "projections.csv", index=False)

    interior = d[(d.internal_x > 0.02) | (d.internal_t > 0.02)]
    print(f"\n    max internal (bounded t): {d.internal_t.max():.5f}")
    print(f"    max internal (x-param)  : {d.internal_x.max():.5f}")
    print(f"    off-star branch length  : {len(interior)} of {len(d)}")

    # A fitted branch length is not evidence.  The objective near saturation
    # is flat along the internal branch, so an optimiser that starts away
    # from the star stops wherever it runs out of gradient without the KL
    # being any lower.  A projection counts as interior only if it BEATS the
    # star fit in KL, so the star is fitted here properly -- four pendants
    # with the internal branch held at zero -- rather than by taking an
    # unconstrained fit and zeroing its internal branch, which leaves the
    # pendants chosen for a tree that had one.
    print("\n  testing whether any off-star point is actually better:")
    rng0 = np.random.default_rng(4)
    sample = d[d.internal_x > 0.02].index if len(interior) else []
    if len(sample) == 0:
        print("    nothing to test; every projection already sits on the star")
    worst = (-np.inf, None)
    for i in list(sample)[:40]:
        row = d.loc[i]
        pend = np.array([float(v) for v in row.pend.split(",")])
        lab = SPLITS[row.topology]
        q = pattern_fast((*pend, row.e_true), SPLITS["12|34"])
        star = np.array([pend[0], pend[2], pend[1], pend[3], 0.0])

        def neg(th):
            if np.any(th <= 1e-12) or np.any(th > 40):
                return 1e3
            return kl(q, pattern_fast(th, lab))

        best = (np.inf, star)
        for t in (star, np.array([*star[:4] * 0.3, 0.1])):
            r = minimize(neg, t, method="L-BFGS-B", bounds=[(1e-12, 40)] * 5)
            if r.fun < best[0]:
                best = (float(r.fun), r.x)
        # the star itself, pendants refitted with the internal branch at zero

        def neg_star(pp):
            return kl(q, pattern_fast((pp[0], pp[2], pp[1], pp[3], 0.0), lab))

        rs = minimize(neg_star, star[:4], method="L-BFGS-B",
                      bounds=[(1e-12, 40)] * 4)
        kl_star = float(rs.fun)
        gain = kl_star - best[0]
        if gain > worst[0]:
            worst = (gain, row.topology)
        print(f"    {row.topology} e_true={row.e_true:7.3f}  KL(star)="
              f"{kl_star:.6e}  KL(best)={best[0]:.6e}  gain={gain:+.2e}")
    if worst[1] is not None:
        print(f"\n    largest KL improvement over the star, over the sample: "
              f"{worst[0]:+.2e}")
        print("    a gain at the level of arithmetic is the optimiser sliding")
        print("    along a flat direction, not a better fit; s30 repeats this")
        print("    properly and over a wider domain.")


if __name__ == "__main__":
    main()
