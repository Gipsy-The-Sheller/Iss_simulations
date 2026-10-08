"""Spectral resolution of the transition matrix -- Part II, Lemma 9 and Lemma 10.

    S = Pi^{1/2} Q Pi^{-1/2} = U Lambda U^T ,   0 = l1 > l2 >= ... >= lK

    P_ij(t) = pi_j + sum_{k>=2} A^(k)_ij e^{lambda_k t}
    A^(k)_ij = pi_i^{-1/2} u_ki u_kj pi_j^{1/2} ,   A^(k)_jj = u_kj^2

and the three scalars of Lemma 10, which every statistic is a linear functional of:

    tau_k = tr A^(k)            (multiplicity-like weight)
    w_k   = sum_j pi_j A^(k)_jj
    xi_k  = < exp(lambda_k d) >  (pair-averaged, and rate-averaged under Gamma)

    Psi = sum_k tau_k xi_k
    wp  = (1 - sum_j pi_j^2) - sum_k w_k xi_k
    Delta = sum_k (w_k - tau_k / K) xi_k

Using this form instead of scipy.linalg.expm makes the code read like the proofs,
gives an analytic dP/dt, and is ~10x faster in the tight loops (pattern
enumeration, design-density inversion).
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass

from .entropy import h


@dataclass
class Spectral:
    pi: np.ndarray           # (K,)   stationary frequencies
    lam: np.ndarray          # (K,)   eigenvalues, descending, lam[0] = 0
    U: np.ndarray            # (K,K)  orthonormal eigenvectors of S (columns)
    A: np.ndarray            # (K,K,K) A[k, i, j]
    tau: np.ndarray          # (K,)
    w: np.ndarray            # (K,)

    # ---------------------------------------------------------------- P(t)
    def P(self, t: float) -> np.ndarray:
        """P(t) at a scalar t (Lemma 9)."""
        e = np.exp(self.lam * t)                       # (K,)
        return np.einsum("k,kij->ij", e, self.A)

    def P_stack(self, ts) -> np.ndarray:
        """P(t) at an array of t -> (len(ts), K, K)."""
        ts = np.asarray(ts, dtype=float)
        e = np.exp(np.outer(ts, self.lam))             # (T, K)
        return np.einsum("tk,kij->tij", e, self.A)

    def dP(self, t: float) -> np.ndarray:
        """dP/dt at scalar t (analytic)."""
        e = self.lam * np.exp(self.lam * t)
        return np.einsum("k,kij->ij", e, self.A)

    # ------------------------------------------------------------ scalars
    def xi(self, d: float | np.ndarray, rates=(1.0,), weights=(1.0,)) -> np.ndarray:
        """xi_k = << exp(lambda_k d) >>, averaged over pairs AND rate classes."""
        d = np.asarray(d, dtype=float)
        out = np.zeros_like(self.lam)
        for rc, wc in zip(rates, weights):
            out += wc * np.exp(np.outer(rc * d, self.lam)).mean(axis=0) if d.ndim else \
                   wc * np.exp(self.lam * rc * float(d))
        return out

    def Psi(self, d, rates=(1.0,), weights=(1.0,)) -> float:
        return float(np.dot(self.tau, self.xi(d, rates, weights)))

    def w_p(self, d, rates=(1.0,), weights=(1.0,)) -> float:
        x = self.xi(d, rates, weights)
        return float((1.0 - np.sum(self.pi ** 2)) - np.dot(self.w, x))

    def Delta(self, d, rates=(1.0,), weights=(1.0,)) -> float:
        x = self.xi(d, rates, weights)
        return float(np.dot(self.w - self.tau / len(self.pi), x))

    @property
    def ratios(self) -> np.ndarray:
        """tau_k / w_k, grouped by distinct eigenvalue -- the Theorem 6 criterion."""
        out = []
        for mu in np.unique(np.round(self.lam[1:], 12)):
            m = np.isclose(self.lam, mu)
            out.append(float(self.tau[m].sum() / self.w[m].sum()))
        return np.array(out)

    def affine_ok(self, tol: float = 1e-9) -> bool:
        r = self.ratios
        return bool(np.allclose(r, r[0], rtol=tol, atol=tol))


def factorize(Q: np.ndarray, pi) -> Spectral:
    """Build the spectral resolution of a reversible Q (Lemma 9)."""
    pi = np.asarray(pi, dtype=float)
    pi = pi / pi.sum()
    K = len(pi)
    S = (np.sqrt(pi)[:, None] * Q) / np.sqrt(pi)[None, :]
    S = (S + S.T) / 2.0                                   # symmetrise away round-off
    lam, U = np.linalg.eigh(S)
    order = np.argsort(lam)[::-1]
    lam, U = lam[order], U[:, order]
    A = np.zeros((K, K, K))
    for k in range(K):
        A[k] = (U[:, k][:, None] * U[:, k][None, :]) * (np.sqrt(pi)[None, :] / np.sqrt(pi)[:, None])
    tau = np.array([np.trace(A[k]) for k in range(K)])
    w = np.array([float(np.dot(pi, np.diag(A[k]))) for k in range(K)])
    # enforce exact zero eigenvalue / unit-rate normalisation
    lam = np.clip(lam, -np.inf, 0.0)
    lam[0] = 0.0
    return Spectral(pi=pi, lam=lam, U=U, A=A, tau=tau, w=w)
