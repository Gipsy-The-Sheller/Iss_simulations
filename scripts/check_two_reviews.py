"""Independent checks of the two review rounds' load-bearing claims.

R1-A  Review 1's internal consistency: for N=4, uniform pi, the affine identity
      forces  d E[Iss] / d(wp) = K(N-1)/(2N h) = 1.08200  exactly.
R1-B  Review 1's mechanism: how large CAN the alpha-sensitivity of Iss be, as a
      function of the internal-branch fraction?  (max over tree length.)
R2-A  Review 2's headline: the Lagrange spectral projector degenerates for K>4.
"""
from __future__ import annotations

import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench import exact, gamma, models, trees                       # noqa: E402
from issbench.entropy import h                                        # noqa: E402

K4 = 4
PI4 = np.ones(4) / 4
Q = models.build_Q("JC69")


def hr(t):
    print("\n" + "=" * 74 + f"\n{t}\n" + "=" * 74)


# ---------------------------------------------------------------- R1-A
hr("R1-A  affine identity forces Delta(Iss) = 1.082 * Delta(wp)  [N=4, uniform pi]")
slope = K4 * (4 - 1) / (2 * 4 * np.log(4))
print(f"  dE[Iss]/d(wp) = K(N-1)/(2N ln K) = {slope:.5f}")
print("  Review 1's table:")
for lab, dIss, dwp in [("f_i=0.10", 0.193, 0.217),
                       ("f_i=0.20", 0.027, 0.028),
                       ("f_i=0.35", 0.0047, 0.0047)]:
    print(f"    {lab}:  dIss={dIss:.4f}  dwp={dwp:.4f}  ->  ratio {dIss/dwp:.3f}"
          f"   (should be {slope:.3f}, off by {100*abs(dIss/dwp-slope)/slope:.0f}%)")


# ---------------------------------------------------------------- R1-B
hr("R1-B  how alpha-sensitive can E[Iss] be?  max over TL, by internal fraction")
print("  quartet, pendants l, internal f_i*TL  (TL = l*(4+...));  alpha in {inf, 0.5}")
print(f"  {'f_i':>6} {'max|Iss(inf)-Iss(0.5)|':>24} {'at TL':>9} {'Iss(inf) there':>15}")
for fi in (0.05, 0.10, 0.20, 0.35, 0.60):
    best, bestTL, bestv = 0.0, 0.0, 0.0
    for TL in np.exp(np.linspace(np.log(0.3), np.log(400), 90)):
        r = fi * TL
        l = (TL - r) / 4.0
        if l <= 0:
            continue
        tr = trees.quartet(l, l, l, l, r)
        a = exact.E_Iss(4, tr, Q, PI4, convention="B")
        b = exact.E_Iss(4, tr, Q, PI4,
                        rates=tuple(gamma.discrete_gamma_rates(0.5, 4)),
                        weights=tuple(gamma.gamma_weights(4)), convention="B")
        if abs(a - b) > best:
            best, bestTL, bestv = abs(a - b), TL, a
    print(f"  {fi:>6.2f} {best:>24.4f} {bestTL:>9.3f} {bestv:>15.4f}")


# ---------------------------------------------------------------- R2-A
hr("R2-A  Lagrange projector vs orthogonal projector, growing K")
rng = np.random.default_rng(7)


def random_reversible(K):
    pi = rng.dirichlet(np.ones(K))
    R = rng.uniform(0.2, 3.0, size=K * (K - 1) // 2)
    Qm = np.zeros((K, K))
    k = 0
    for i in range(K):
        for j in range(i + 1, K):
            Qm[i, j] = R[k] * pi[j]
            Qm[j, i] = R[k] * pi[i]
            k += 1
    for i in range(K):
        Qm[i, i] = -Qm[i].sum()
    Qm /= -np.sum(pi * np.diag(Qm))
    return Qm, pi


print(f"  {'K':>4} {'sum(tau) Lagrange':>18} {'sum(tau) orthogonal':>20} "
      f"{'max|w err| Lagrange':>21}")
for K in (4, 5, 6, 8, 12, 20):
    Qm, pi = random_reversible(K)
    S = (np.sqrt(pi)[:, None] * Qm) / np.sqrt(pi)[None, :]
    S = (S + S.T) / 2
    lam, U = np.linalg.eigh(S)
    mus = np.unique(np.round(lam, 10))
    tL = wL = 0.0
    for mu in mus:
        if abs(mu) < 1e-12:
            continue
        Pm = np.eye(K)
        for nu in mus:
            if nu != mu:
                Pm = Pm @ (S - nu * np.eye(K)) / (mu - nu)
        tL += np.trace(Pm).real
        wL += np.trace(np.diag(pi) @ Pm).real
    Porth = (U * (np.isclose(lam, mu))) @ U.T if False else None
    tO = 0.0
    wO = 0.0
    for mu in mus:
        if abs(mu) < 1e-12:
            continue
        m = np.isclose(lam, mu, atol=1e-9)
        P = U[:, m] @ U[:, m].T
        tO += np.trace(P)
        wO += float(pi @ np.diag(P))
    print(f"  {K:>4} {tL:>18.4g} {tO:>20.4g} {abs(wL-wO):>21.4g}")
