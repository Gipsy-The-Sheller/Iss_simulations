"""The closed form for a general REV model, written out and checked.

The question this answers: with the six exchangeabilities t_AC, t_AG, t_AT,
t_CG, t_CT, t_GT and the four frequencies pi_A, pi_C, pi_G, pi_T, can E[Iss] be
written down?  On an equal-length star, yes, and this script is that formula,
evaluated two independent ways.

The construction, in the manuscript's notation (Sec. 2.2).  With Q_ij = t_ij
pi_j and c = sum_{i<j} 2 pi_i pi_j t_ij (so one substitution per unit time),

    S = Pi^{1/2} Q Pi^{-1/2},   S_ij = sqrt(pi_i pi_j) t_ij / c   (i != j)
                                S_ii = -sum_{j != i} pi_j t_ij / c

is symmetric, has S sqrt(pi) = 0, and therefore S = U Lambda U^T with
orthonormal U and u_1 = sqrt(pi).  Note the diagonal is NOT the one that makes
the row sums vanish -- only S sqrt(pi) = 0 is required, and the row sums are
generally non-zero.  Then

    P_ij(t) = pi_j + sum_{k>=2} pi_i^{-1/2} u_ki u_kj pi_j^{1/2} e^{lambda_k t}.

On a star of height d, given the root state rho and the rate class c the N
leaf states are conditionally i.i.d. with law q = P_{rho.}(r_c d), so the
column counts are Mult(N, q) and the exact multinomial entropy applies:

    E[Iss] = 1/H_FSS * sum_c w_c sum_rho pi_rho
             [ ln N - (1/N) sum_j sum_{m=1}^{N} m ln m C(N,m)
               q_{rho j}^m (1 - q_{rho j})^{N-m} ],

with H_FSS = E_{Mult(N,pi)}[h] under convention (A) and h(pi) under (B).  The
whole parameter dependence is carried by three numbers per rate class -- the
non-zero eigenvalues -- and their eigenvectors.  Nothing is iterative, no
matrix exponential is formed, and this holds for any K.

What it is NOT is a formula in the p-distance.  The last section prints
tau_mu/w_mu, which Theorem 13 requires to be constant across the distinct
non-zero eigenvalues; for a generic REV model it is not, so E[Iss] depends on
the whole 9-dimensional parameter vector and no reduction to p-bar exists.
"""
from __future__ import annotations

import sys
from math import comb
from pathlib import Path

import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench import exact, trees                                    # noqa: E402
from issbench.entropy import h as entropy                            # noqa: E402
from issbench.models import build_Q                                  # noqa: E402
from issbench.spectral import factorize                              # noqa: E402

PAIRS = [(i, j) for i in range(4) for j in range(i + 1, 4)]

# states are A, C, G, T and build_Q takes R ordered AC, AG, AT, CG, CT, GT, so
# this is exactly the parameter list of the question
R = np.array([2.1, 4.7, 0.9, 1.3, 3.6, 1.0])
PI = np.array([0.31, 0.24, 0.19, 0.26])
N = 8
D = 0.6
RATES = (0.25, 0.9, 2.4)
WEIGHTS = (1 / 3, 1 / 3, 1 / 3)


def sym_spectrum(R: np.ndarray, pi: np.ndarray):
    """S = Pi^{1/2} Q Pi^{-1/2} from the exchangeabilities alone.

    Returns (S, lam, U, c): eigenvalues descending, U's columns orthonormal.
    """
    K = len(pi)
    c = sum(2 * pi[i] * pi[j] * R[k] for k, (i, j) in enumerate(PAIRS))
    S = np.zeros((K, K))
    for k, (i, j) in enumerate(PAIRS):
        S[i, j] = S[j, i] = np.sqrt(pi[i] * pi[j]) * R[k] / c
    for i in range(K):
        # q_ii = -sum_{j != i} t_ij pi_j, scaled by the same c
        S[i, i] = -sum(pi[j] * R[k] for k, (a, b) in enumerate(PAIRS)
                       if i in (a, b) for j in ((b,) if a == i else (a,))) / c
    lam, U = np.linalg.eigh(S)                       # ascending
    order = np.argsort(-lam)
    return S, lam[order], U[:, order], c


def closed_form(R, pi, N, d, rates, weights, convention="A", base="nats"):
    """E[Iss] on an equal-length star, from the formula above."""
    pi = np.asarray(pi, float)
    _, lam, U, _ = sym_spectrum(R, pi)
    K = len(pi)
    m = np.arange(1, N + 1)
    ln_coef = m * np.log(m) * np.array([comb(N, int(mm)) for mm in m])

    num = 0.0
    for r, w in zip(rates, weights):
        for rho in range(K):
            q = np.array([
                pi[j] + sum(pi[rho] ** -0.5 * U[rho, k] * U[j, k] * pi[j] ** 0.5
                            * np.exp(lam[k] * r * d) for k in range(1, K))
                for j in range(K)])
            num += w * pi[rho] * (np.log(N) - sum(
                float(ln_coef @ (qj ** m * (1 - qj) ** (N - m)))
                for qj in np.clip(q, 0.0, 1.0)) / N)
    if base == "bits":
        num /= np.log(2)
    hfss = (exact.multinomial_entropy(N, pi, base) if convention == "A"
            else entropy(pi, base))
    return num / hfss


def main() -> None:
    Q = build_Q("GTR", pi=PI, R=R)
    S, lam, U, c = sym_spectrum(R, PI)
    print("exchangeabilities  AC AG AT CG CT GT =", R)
    print("frequencies        A C G T =", PI, " (sum %.3f)" % PI.sum())
    print(f"normalisation c = sum 2 pi_i pi_j t_ij = {c:.6f}")

    print("\nS, off-diagonals sqrt(pi_i pi_j) t_ij / c:")
    with np.printoptions(precision=5, suppress=True):
        print(S)
    print("  eigenvalues      ", np.round(lam, 8))
    print("  S sqrt(pi)       ", np.round(S @ np.sqrt(PI), 12),
          "  (must be 0)")
    print("  eigenvalues of Q ", np.round(np.sort(np.linalg.eigvals(Q).real)[::-1], 8))
    print("  max|S - S^T| =", np.abs(S - S.T).max())

    print("\nThe spectral resolution must reproduce exp(Qt) before it is used "
          "to build q:")
    sp = factorize(Q, PI)
    for t in (0.05, 0.6, 3.0):
        print(f"  t = {t:<5} max |P_spectral(t) - exp(Qt)| = "
              f"{np.abs(sp.P(t) - expm(Q * t)).max():.2e}")

    print(f"\nE[Iss] on an equal-length star, N={N}, d={D}, rates={RATES}:")
    for conv in ("A", "B"):
        mine = closed_form(R, PI, N, D, RATES, WEIGHTS, convention=conv)
        ref = exact.E_Iss(N, trees.star(N, D), Q, PI, RATES, WEIGHTS,
                          convention=conv)
        print(f"  convention ({conv}):  formula {mine:.12f}   "
              f"enumeration {ref:.12f}   |diff| {abs(mine - ref):.2e}")

    print("\nTheorem 13, for this parameter set:")
    print("  distinct tau/w ratios:", [round(float(r), 4) for r in sp.ratios],
          "->", "CONSTANT" if sp.affine_ok() else "NOT constant")
    print("  so Psi is not affine in p-bar here, and the formula above is as "
          "short as\n  E[Iss] gets for this model: it is a function of the "
          "whole parameter vector,\n  not of the mean p-distance.")

    print("\nThe same formula under the two cases where the dependence "
          "does collapse:")
    cases = [
        ("uniform pi, any R",
         np.ones(4) / 4, R),
        ("F81 (all t equal), skewed pi",
         PI, np.ones(6)),
    ]
    for label, pi, r in cases:
        q = build_Q("GTR", pi=pi, R=r)
        s = factorize(q, pi)
        mine = closed_form(r, pi, N, D, RATES, WEIGHTS)
        ref = exact.E_Iss(N, trees.star(N, D), q, pi, RATES, WEIGHTS)
        print(f"  {label:32s} tau/w {[round(float(x), 3) for x in s.ratios]}"
              f"  {'affine' if s.affine_ok() else 'NOT affine'}"
              f"  |formula - enum| {abs(mine - ref):.1e}")

    estimation()


def estimation() -> None:
    """Can the index be inverted to a parameter?

    On an equal-length star, E[Iss] is monotone increasing in d for a fixed
    model -- it runs from 0 at d = 0, where every column is a single state, to
    its ceiling as the columns randomise -- so the index IS invertible, to a
    tree length, under the assumption that pi, the exchangeabilities and the
    shape are already known.  This section inverts it at a fixed reading and
    shows what that assumption costs: the same value implies tree lengths that
    differ by a factor of several across models the alignment alone cannot
    tell apart.
    """
    from issbench.gamma import discrete_gamma_rates

    # a grid that reaches the asymptote even for the slowest rate class: at
    # alpha = 0.5 the slow category has rate 0.03, so d = 6 leaves it far from
    # randomised and the curve has not levelled off there
    grid = np.concatenate([np.linspace(1e-3, 6.0, 2000),
                           np.geomspace(6.0, 400.0, 1500)])
    grid = np.unique(grid)
    q_uni = np.ones(4) / 4
    g_rates = np.array(discrete_gamma_rates(0.5, 4))
    g_w = np.full(4, 0.25)
    models = [
        ("JC69, uniform", build_Q("JC69"), q_uni, (1.0,), (1.0,)),
        ("GTR, uniform pi", build_Q("GTR", pi=q_uni, R=R), q_uni,
         (1.0,), (1.0,)),
        ("GTR, skewed pi", build_Q("GTR", pi=PI, R=R), PI, (1.0,), (1.0,)),
        ("GTR, skewed pi, G4(0.5)", build_Q("GTR", pi=PI, R=R), PI,
         tuple(g_rates), tuple(g_w)),
    ]
    curve = {label: np.array([exact.E_Iss(N, d, q, pi, r, w)
                              for d in grid])
             for label, q, pi, r, w in models}

    print("\nInverting the index: what tree length does a fixed reading "
          f"imply?  (N={N})")
    print("  " + "Iss".ljust(8) + "".join(k.ljust(24) for k in curve))
    for target in (0.85, 0.75, 0.65):
        row = f"  {target:<8.2f}"
        for label, y in curve.items():
            # E[Iss] increases in d, so y is already the increasing abscissa;
            # a reading past the ceiling has no tree length to invert to
            row += (f"{float(np.interp(target, y, grid)):<24.3f}"
                    if target <= y.max() else f"{'unattainable':<24s}")
        print(row)
    print("  ceilings: " + ", ".join(f"{k} {v.max():.4f}"
                                     for k, v in curve.items()))
    print("  a spread of this size between models the alignment alone cannot "
          "\n  distinguish is the mis-calibration stated as an estimator: the "
          "index\n  identifies a tree length only under a model assumed in "
          "advance.")


if __name__ == "__main__":
    main()
