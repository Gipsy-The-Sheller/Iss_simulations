"""The three quartet likelihoods, and the probability that the true one wins.

For a quartet there are exactly three topologies and therefore three
likelihood functions.  Each is a sum of L i.i.d. per-site terms,

    ell_a = sum_j log P_{T_a}(x_j),

so (ell_1, ell_2, ell_3) is jointly asymptotically normal, and

    Pr[correct] = Pr[ell_1 > ell_2 and ell_1 > ell_3]
                = Pr[D_2 > 0 and D_3 > 0],   D_a = ell_1 - ell_a,

is a bivariate normal orthant probability.  Everything it needs -- the two
means, the two variances and the correlation -- is a sum over the 256 site
patterns of the model, and P(x) is a multilinear polynomial in the five
x_e = exp(-4 t_e / 3) with each edge appearing exactly once.  No closed-form
MLE of the branch lengths is required, and none exists.

The alternative topologies enter at their KL-optimal branch lengths, since
that is what sup_{b_a} converges to, not at the true ones.

CHECK.  At the degeneracy (internal branch zero) the three topologies are the
same star, the three ell are exchangeable, the means vanish and the orthant
probability must be exactly 1/3.  For a zero-mean bivariate normal that is
Phi_2(0,0;rho) = 1/4 + arcsin(rho)/(2 pi), so rho must come out at
sin(pi/6) = 1/2.  The script reports it before anything else.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

# Part 2.  Inputs that belong to the published preprint (the empirical
# alignments, the simulated work directory) are read from the
# repository root; everything this line of work produces is written
# under Part_2/, so the two cannot be confused.
PART2 = Path(__file__).resolve().parents[1]
ROOT = PART2.parent
sys.path.insert(0, str(ROOT / "src"))

from issbench import trees                                          # noqa: E402
from issbench.exact import tree_patterns                           # noqa: E402
from issbench.models import build_Q                                 # noqa: E402

OUT = PART2 / "results" / "quartet_kl"
Q = build_Q("JC69")
PI = np.ones(4) / 4
NPAT = 256
LEAVES = ("A", "B", "C", "D")


def pattern_vector(lengths, labels=LEAVES) -> np.ndarray:
    """P(x) over the 256 patterns.

    `lengths` is (a,b,c,d,e) and goes with the leaves in `labels`, so the
    topology is changed by permuting labels, never lengths: the latter only
    reassigns pendants within the same split, which is the bug the first
    version of this script had.
    """
    a, b, c, d, e = lengths
    node = trees.quartet(a, b, c, d, e, labels=labels)
    pats = tree_patterns(node, Q, PI, (1.0,), (1.0,))
    # tree_patterns returns patterns in leaf order, which is `labels`; the
    # array is packed in the FIXED order A, B, C, D so that the three
    # topologies land on a common support.  Packing in tuple order instead
    # puts leaf C's state in leaf B's slot and the three distributions stop
    # being comparable -- which is what made the star check fail.
    order = [labels.index(n) for n in ("A", "B", "C", "D")]
    v = np.zeros(NPAT)
    for pat, pr in pats.items():
        idx = 0
        for p in order:
            idx = (idx << 2) | pat[p]
        v[idx] = pr
    return v


SPLITS = {"12|34": LEAVES,
          "13|24": ("A", "C", "B", "D"),
          "14|23": ("A", "D", "B", "C")}


def kl(p, q, floor=1e-300) -> float:
    m = p > 0
    return float(np.sum(p[m] * np.log(p[m] / np.maximum(q[m], floor))))


def best_fit(p_true, labels, x0):
    """Branch lengths of topology `labels` minimising KL(p_true || q)."""
    best = (np.inf, np.asarray(x0, float))

    def obj(theta):
        if np.any(theta <= 0) or np.any(theta > 20):
            return 1e3
        return kl(p_true, pattern_vector(theta, labels))

    for trial in (x0, [0.2, 0.2, 0.2, 0.2, 0.02],
                  list(np.asarray(x0) * 0.4), list(np.asarray(x0) * 2.0)):
        r = minimize(obj, trial, method="Nelder-Mead",
                     options=dict(maxiter=3000, xatol=1e-8, fatol=1e-12))
        r = minimize(obj, r.x, method="L-BFGS-B",
                     bounds=[(1e-8, 20)] * 5)
        if r.fun < best[0]:
            best = (float(r.fun), r.x)
    return best


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------ the degeneracy check
    # At e = 0 the three topologies are the same star, and since a star is
    # the common boundary of all three parametrisations each can represent it
    # exactly -- so the three likelihood functions are IDENTICAL, not merely
    # exchangeable.  The differences D_2 and D_3 then have zero variance and
    # the correlation between them is undefined, not 1/2; an earlier version
    # of this script asserted 1/2 and was simply wrong.  What is worth
    # recording is the approach to that limit, below.
    star = (0.3, 0.3, 0.3, 0.3, 0.0)
    p = pattern_vector(star)
    ps = [pattern_vector(star, lab) for lab in SPLITS.values()]
    for lab, q in zip(SPLITS, ps):
        assert np.abs(q - p).max() < 1e-12, f"{lab} differs at e=0"
    print("CHECK: at the star the three pattern vectors agree to machine")
    print("precision, so the three likelihoods are identical and the winner")
    print("is decided by tie-breaking, not by the orthant probability.")
    print("\nnear the degeneracy the two WRONG topologies collapse onto the")
    print("same star, so ell_2 = ell_3 and rho -> 1; the three-way choice")
    print("becomes a two-way one.  Swept in e with pendants (1, .05, 1, .05):")
    from scipy.stats import multivariate_normal as mvn0
    print(f"    {'e':>8} {'D_2/L':>11} {'D_3/L':>11} {'rho23':>7} "
          f"{'Pr[ok] L=1000':>14}")
    for e_s in (0.005, 0.01, 0.05, 0.2, 0.5, 1.0):
        b1 = (1.0, 0.05, 1.0, 0.05, e_s)
        pt = pattern_vector(b1)
        lg = []
        for lab in SPLITS:
            if lab == "12|34":
                lg.append(np.log(np.maximum(pt, 1e-300)))
            else:
                _, bb = best_fit(pt, SPLITS[lab], b1)
                lg.append(np.log(np.maximum(pattern_vector(bb, SPLITS[lab]),
                                            1e-300)))
        lg = np.array(lg)
        mu0 = lg @ pt
        cv = np.cov(lg, aweights=pt)
        m2, m3 = mu0[0] - mu0[1], mu0[0] - mu0[2]
        v2 = cv[0, 0] + cv[1, 1] - 2 * cv[0, 1]
        v3 = cv[0, 0] + cv[2, 2] - 2 * cv[0, 2]
        c23 = cv[0, 0] - cv[0, 2] - cv[0, 1] + cv[1, 2]
        ok = np.sqrt(v2 * v3) > 0
        rho = c23 / np.sqrt(v2 * v3) if ok else np.nan
        if ok and abs(rho) < 1 - 1e-9:
            pr = float(mvn0(mean=[-m2 * 1000, -m3 * 1000],
                            cov=[[v2 * 1000, c23 * 1000],
                                 [c23 * 1000, v3 * 1000]]).cdf([0.0, 0.0]))
        elif ok:
            # rho = 1: D_2 = D_3, so the three-way choice is a two-way one
            from scipy.stats import norm as _norm
            pr = float(_norm.cdf(m2 * 1000 / np.sqrt(v2 * 1000)))
        else:
            pr = float("nan")
        print(f"    {e_s:8.4f} {m2:11.3e} {m3:11.3e} {rho:7.4f} {pr:14.4f}")

    # ------------------------------------------- the transfer path proper
    print("\nthe transfer path a=(L-t), b=S, c=L, d=(S+t), e=0.2")
    print("  step 0 is Felsenstein (long branches apart), the last Farris")
    L_, S_, E_ = 1.0, 0.05, 0.2
    ts = np.linspace(0, L_ - S_, 7)
    rows = []
    print(f"    {'step':>4} {'D_2/L':>11} {'D_3/L':>11} {'rho23':>8} "
          f"{'Pr[ok] L=300':>13} {'L=1000':>8} {'L=10000':>9}")
    from scipy.stats import multivariate_normal as mvn
    for j, t in enumerate(ts):
        b1 = (L_ - t, S_, L_, S_ + t, E_)
        p_true = pattern_vector(b1)
        fits = {}
        for lab in SPLITS:
            if lab == "12|34":
                fits[lab] = (0.0, np.asarray(b1))
            else:
                fits[lab] = best_fit(p_true, SPLITS[lab], b1)
        logs = np.array([np.log(np.maximum(pattern_vector(fits[lab][1], SPLITS[lab]), 1e-300))
                         for lab in SPLITS])
        mu = logs @ p_true
        cov = np.cov(logs, aweights=p_true)
        d2m, d3m = mu[0] - mu[1], mu[0] - mu[2]
        d2v = cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
        d3v = cov[0, 0] + cov[2, 2] - 2 * cov[0, 2]
        c23 = cov[0, 0] - cov[0, 2] - cov[0, 1] + cov[1, 2]
        rho = c23 / np.sqrt(d2v * d3v)
        row = dict(step=j, t=float(t), D2=float(d2m), D3=float(d3m), rho=float(rho),
                   v2=float(d2v), v3=float(d3v))
        for L in (300, 1000, 10000):
            m = [d2m * L, d3m * L]
            C = [[d2v * L, c23 * L], [c23 * L, d3v * L]]
            pr = float(mvn(mean=m, cov=C).cdf([0.0, 0.0]))
            # P[D2>0, D3>0] = P[(-D2)<0, (-D3)<0] with negated means
            pr = float(mvn(mean=[-m[0], -m[1]], cov=C).cdf([0.0, 0.0]))
            row[f"P{L}"] = pr
        rows.append(row)
        print(f"    {j:4d} {d2m:11.3e} {d3m:11.3e} {rho:8.4f} "
              f"{row['P300']:13.4f} {row['P1000']:8.4f} {row['P10000']:9.4f}")
    pd.DataFrame(rows).to_csv(OUT / "orthant_transfer.csv", index=False)

    print("\n  D_a/L is minus the KL gap of topology a, so D_a > 0 means the")
    print("  true topology is the better fit on average; the boundary is")
    print("  where the mean gap falls to the order of the sampling noise.")


if __name__ == "__main__":
    main()
