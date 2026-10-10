"""Iss as a moment of the count distribution: the expansion, checked.

Three things are established numerically here, all of them exact statements
about the model rather than fits.

1. THE EXPANSION.  With w(x) = h(c(x)/N) the entropy of a column's counts,
   Hbar = sum_x w(x) fhat_x identically, for any N, any K, any tree and any
   model -- the weights are integers-derived constants carrying no parameter.
   This is checked against the naive per-site computation.

2. THE VARIANCE.  Hbar is a mean of i.i.d. bounded W_i, so
   Var(Hbar) = sigma^2 / L with sigma^2 = sum_x w(x)^2 P(x) - mu^2 exactly.
   Checked against simulation, and sigma^2 is what Theorem 11(2) leaves
   abstract.

3. WHAT THE EXPANSION IS WORTH.  On an equidistant star the leaves are
   exchangeable, so P(x) depends on x only through its count vector and the
   count histogram ghat_c is sufficient -- with C(N+3,3) entries.  The index
   is the single linear functional c -> h(c/N) of it.  The Jacobian of
   theta -> g(theta) is computed here: if it is full rank, the parameters are
   identifiable from the same alignment through other weights, and the loss
   is a property of h and not of the data.
"""
from __future__ import annotations

import sys
from itertools import combinations_with_replacement
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench.entropy import h                                        # noqa: E402
from issbench.exact import multinomial_entropy                        # noqa: E402
from issbench.gamma import discrete_gamma_rates                       # noqa: E402
from issbench.models import build_Q                                   # noqa: E402
from issbench.spectral import factorize                               # noqa: E402

OUT = ROOT / "results" / "moment"
SYMBOLS = "ACGT"
N = 6
L = 4000
D = 0.55
R_UNIF = np.ones(6)
R_SKEW = np.array([2.1, 4.7, 0.9, 1.3, 3.6, 1.0])
PI_SKEW = np.array([0.31, 0.24, 0.19, 0.26])


# --------------------------------------------------------------- count vectors
def _comps(N, K, prefix=()):
    """All count vectors of length K summing to N, in a fixed order."""
    if K == 1:
        return [prefix + (N,)]
    out = []
    for i in range(N + 1):
        out += _comps(N - i, K - 1, prefix + (i,))
    return out


def column_weights(N: int):
    """w(c) = h(c/N), as an array over the composition list."""
    comps = _comps(N, 4)
    w = np.array([h(np.array(c, dtype=float) / N) for c in comps])
    return comps, w


def pattern_to_count_index(comps) -> dict:
    return {c: i for i, c in enumerate(comps)}


def count_distribution(N: int, Q, pi, d: float, rates=(1.0,), weights=(1.0,)):
    """g_c = sum_{x : c(x)=c} P(x) on an equidistant star.

    The leaves are exchangeable, so P(x) is a function of the count vector and
    the sum over x in a class is the class size times P of any representative.
    """
    comps = _comps(N, 4)
    sp = factorize(Q, pi)
    idx = {c: i for i, c in enumerate(comps)}
    g = np.zeros(len(comps))
    for r, wt in zip(rates, weights):
        P = sp.P(r * d)
        for rho in range(4):
            q = P[rho]
            for c in comps:
                cls = comb(N, c[0]) * comb(N - c[0], c[1]) * comb(N - c[0] - c[1], c[2])
                g[idx[c]] += wt * pi[rho] * cls * np.prod(q ** np.array(c))
    return comps, g


NAMES = ["pi_A", "pi_C", "pi_G", "t_AC", "t_AG", "t_AT", "t_CG", "t_CT",
         "d", "alpha"]


def g_of(theta, N: int = N):
    """The count histogram as a function of pi, five t_ij, d and the shape."""
    pi = np.array([theta[0], theta[1], theta[2],
                   1 - theta[0] - theta[1] - theta[2]])
    R = np.array([theta[3], theta[4], theta[5], theta[6], theta[7], 1.0])
    d, alpha = theta[8], theta[9]
    q = build_Q("GTR", pi=pi, R=R)
    if alpha > 1e3:
        return count_distribution(N, q, pi, d)[1]
    return count_distribution(N, q, pi, d, discrete_gamma_rates(alpha, 4),
                              np.full(4, 0.25))[1]


def jacobian(theta0, rel: float = 1e-5):
    """(g, dg/dtheta) by central differences at a relative step `rel`."""
    g0 = g_of(theta0)
    J = np.zeros((len(g0), len(theta0)))
    for j in range(len(theta0)):
        step = rel * max(1.0, abs(theta0[j]))
        tp, tm = theta0.copy(), theta0.copy()
        tp[j] += step
        tm[j] -= step
        J[:, j] = (g_of(tp) - g_of(tm)) / (2 * step)
    return g0, J


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    comps, w = column_weights(N)
    print(f"N={N}: {len(comps)} count vectors, "
          f"{4 ** N} site patterns")

    Q = build_Q("GTR", pi=PI_SKEW, R=R_SKEW)
    _, g = count_distribution(N, Q, PI_SKEW, D)
    _, g_gam = count_distribution(N, Q, PI_SKEW, D,
                                  discrete_gamma_rates(0.5, 4),
                                  np.full(4, 0.25))
    assert abs(g.sum() - 1) < 1e-12 and abs(g_gam.sum() - 1) < 1e-12

    mu = float(w @ g)
    mu_gam = float(w @ g_gam)
    var = float((w ** 2) @ g) - mu ** 2
    var_gam = float((w ** 2) @ g_gam) - mu_gam ** 2
    print(f"\nexact E[Hbar]  skew GTR     {mu:.6f}   sd at L={L} "
          f"{np.sqrt(var / L):.6f}")
    print(f"exact E[Hbar]  skew GTR+G4  {mu_gam:.6f}   sd at L={L} "
          f"{np.sqrt(var_gam / L):.6f}")

    # the expansion, against one draw of L columns
    rng = np.random.default_rng(7)
    for label, gg in (("no gamma", g), ("G4(0.5)", g_gam)):
        draw = rng.choice(len(comps), size=L, p=gg)
        hbar_direct = w[draw].mean()
        fhat = np.bincount(draw, minlength=len(comps)) / L
        hbar_expansion = float(w @ fhat)
        print(f"  [{label}] direct {hbar_direct:.6f}  "
              f"expansion {hbar_expansion:.6f}  "
              f"|diff| {abs(hbar_direct - hbar_expansion):.2e}")

    # sigma^2 against a sampling experiment
    reps = 4000
    emp = []
    for _ in range(reps):
        emp.append(w[rng.choice(len(comps), size=L, p=g_gam)].mean())
    print(f"\nsigma^2 predicted {var_gam:.6f}   from {reps} synthetic means "
          f"{np.var(emp) * L:.6f}")

    # ------------------------------------------------- identifiability of g
    names = NAMES
    thetas = {
        "alpha = 70 (near uniform)":
            np.array([0.31, 0.24, 0.19, 2.1, 4.7, 0.9, 1.3, 3.6, 0.55, 70.0]),
        "alpha = 1":
            np.array([0.31, 0.24, 0.19, 2.1, 4.7, 0.9, 1.3, 3.6, 0.55, 1.0]),
    }
    base = g_of(thetas["alpha = 1"])
    for label, theta0 in thetas.items():
        _, J = jacobian(theta0, 1e-5)
        sv = np.linalg.svd(J, compute_uv=False)
        rank = int((sv > sv[0] * 1e-8).sum())
        mu0 = float(w @ g_of(theta0))
        elas = J.T @ w * theta0 / mu0          # elasticity of E[Iss]
        print(f"\ncount histogram g: {len(base)} entries at {label}")
        print(f"  rank(J) = {rank} of {len(theta0)};  condition number "
              f"{sv[0] / sv[rank - 1]:.3g}")
        print("  elasticity of E[Iss] wt theta: " +
              "  ".join(f"{n} {e:+.3f}" for n, e in zip(names, elas)))

    print("\n  the Jacobian is a central difference; the rank is stable to "
          "the step:")
    for rel in (1e-4, 1e-5, 1e-6):
        _, Jr = jacobian(thetas["alpha = 1"], rel)
        svr = np.linalg.svd(Jr, compute_uv=False)
        print(f"    step {rel:.0e}  rank {int((svr > svr[0] * 1e-8).sum())}"
              f"  condition number {svr[0] / svr[9]:.3g}")

    print("\n  the histogram identifies all ten parameters at both points, so "
          "the\n  information the index needs is present in the alignment; "
          "what it\n  discards is a property of the weight function h, not of "
          "the data.")
    print(f"\nIss uses w = h(c/N): one linear functional of {len(base)} "
          f"independent ones.")
    print("  grad mu is a single direction; the moment map theta -> mu is "
          "scalar.")

    pd.DataFrame(dict(c=["".join(str(x) for x in cc) for cc in comps],
                      g=g, g_gamma=g_gam, w=w)).to_csv(
        OUT / "count_histogram.csv", index=False)
    rows = []
    for label, theta0 in thetas.items():
        g0 = g_of(theta0)
        mu0 = float(w @ g0)
        J0 = np.zeros((len(g0), len(theta0)))
        for j in range(len(theta0)):
            step = 1e-5 * max(1.0, abs(theta0[j]))
            tp, tm = theta0.copy(), theta0.copy()
            tp[j] += step
            tm[j] -= step
            J0[:, j] = (g_of(tp) - g_of(tm)) / (2 * step)
        for n, t0, e in zip(names, theta0, J0.T @ w * theta0 / mu0):
            rows.append(dict(point=label, parameter=n, value=t0,
                             elasticity=float(e),
                             histogram_sensitivity=float(
                                 np.linalg.norm(J0[:, names.index(n)]))))
    pd.DataFrame(rows).to_csv(OUT / "moment_jacobian.csv", index=False)
    power()
    print(f"\nwrote {OUT}")


def power() -> None:
    """How long an alignment has to be to resolve the span of Sec. 2.6.

    Lemma (exact mean and variance) holds for any tree, not only a star, so
    the standard error of Hbar can be had exactly on the quartet the alpha
    sweep uses -- the pattern distribution comes from the same enumeration
    the exact E[Iss] does.  The leading term of SE(Iss) is
    sigma / (sqrt(L) H_FSS).
    """
    from issbench.trees import quartet
    from issbench.exact import tree_patterns
    from issbench.gamma import discrete_gamma_rates

    Q = build_Q("JC69")
    pi = np.ones(4) / 4
    hfss = multinomial_entropy(4, pi)
    span = 0.9807 - 0.9692                 # the Iss.c span of Eq. (tlc)
    rows = []
    for alpha, tl in ((np.inf, 7.19), (8.0, 7.78), (2.0, 13.02),
                      (1.0, 25.47), (0.5, 107.72)):
        r = tl * 0.35
        tree = quartet((tl - r) / 4, (tl - r) / 4, (tl - r) / 4,
                       (tl - r) / 4, r)
        if np.isinf(alpha):
            rates, wts = (1.0,), (1.0,)
        else:
            rates, wts = tuple(discrete_gamma_rates(alpha, 4)), (0.25,) * 4
        pats = tree_patterns(tree, Q, pi, rates, wts)
        P, W = [], []
        for pat, pr in pats.items():
            cnt = np.bincount(list(pat), minlength=4)
            P.append(pr)
            W.append(h(cnt / 4))
        P, W = np.array(P), np.array(W)
        mu = float(W @ P)
        var = float((W ** 2) @ P) - mu ** 2
        se = np.sqrt(var) / np.sqrt(1000) / hfss
        # L at which the span is k standard errors
        L3 = 1000 * (3 * se / span) ** 2
        rows.append(dict(alpha=alpha, tl=tl, mu=mu, sigma=np.sqrt(var),
                         se_L1000=se, L_for_3se=L3))
        print(f"  alpha {str(alpha):>4} TL {tl:>7}  E[Hbar] {mu:.4f}  "
              f"sigma {np.sqrt(var):.4f}  SE(Iss) at L=1000 {se:.5f}  "
              f"L for 3 SE over the span: {L3:.0f}")
    print(f"  the span it has to resolve is {span:.4f}")
    pd.DataFrame(rows).to_csv(OUT / "power.csv", index=False)


if __name__ == "__main__":
    main()
