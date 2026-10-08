"""Reversible nucleotide substitution models, normalised to unit rate.

Normalisation:  -sum_j pi_j Q_jj = 1  (one expected substitution per unit time),
so branch lengths are in expected substitutions per site -- matching PAML /
IQ-TREE conventions and Part I assumption A2.

States are ordered A, C, G, T = 0, 1, 2, 3.  Transitions are A<->G and C<->T.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import expm
from typing import Literal

SYMBOLS = "ACGT"
TS = frozenset({(0, 2), (2, 0), (1, 3), (3, 1)})

ModelName = Literal["JC69", "K80", "K2P", "F81", "F84", "HKY85", "TN93", "GTR"]


def _normalise(Q: np.ndarray, pi: np.ndarray) -> np.ndarray:
    Q = np.asarray(Q, dtype=float).copy()
    np.fill_diagonal(Q, 0.0)
    for i in range(4):
        Q[i, i] = -Q[i].sum()
    rate = -float(np.sum(pi * np.diag(Q)))
    return Q / rate


def build_Q(model: ModelName, pi=None, kappa: float = 1.0, **kw) -> np.ndarray:
    """Return the normalised rate matrix Q.

    pi    : stationary frequencies.  Uniform for JC69/K80; required for the rest.
    kappa : transition/transversion rate ratio.  kappa = 1 gives F81 when
            pi is supplied (the maximally degenerate spectrum).
    """
    pi = np.ones(4) / 4 if pi is None else np.asarray(pi, dtype=float)
    pi = pi / pi.sum()
    m = model.upper()

    if m == "JC69":
        Q = np.ones((4, 4))
    elif m in ("K80", "K2P"):
        Q = np.ones((4, 4))
        for i, j in TS:
            Q[i, j] = kappa
    elif m in ("F81", "HKY85", "TN93", "F84"):
        Q = np.tile(pi, (4, 1))               # transversions: q_ij = pi_j
        for i, j in TS:
            Q[i, j] = pi[j] * kappa           # transitions: q_ij = kappa * pi_j
        if m == "TN93":                       # two transition rates
            kappa2 = kw.get("kappa2", kappa)
            Q[0, 2] = pi[2] * kappa
            Q[2, 0] = pi[0] * kappa
            Q[1, 3] = pi[3] * kappa2
            Q[3, 1] = pi[1] * kappa2
        if m == "F84":                        # Felsenstein 1984 parameterisation
            piY = pi[1] + pi[3]               # C + T
            piR = pi[0] + pi[2]               # A + G
            Q = np.tile(pi, (4, 1))
            Q[0, 2] = pi[2] * (1.0 + kappa / piY)
            Q[2, 0] = pi[0] * (1.0 + kappa / piY)
            Q[1, 3] = pi[3] * (1.0 + kappa / piR)
            Q[3, 1] = pi[1] * (1.0 + kappa / piR)
    elif m == "GTR":
        R = kw.get("R")
        if R is None:
            raise ValueError("GTR needs the 6 exchangeabilities via R=")
        R = np.asarray(R, dtype=float)
        Q = np.zeros((4, 4))
        k = 0
        for i in range(4):
            for j in range(i + 1, 4):
                Q[i, j] = R[k] * pi[j]
                Q[j, i] = R[k] * pi[i]
                k += 1
    else:
        raise ValueError(f"unknown model {model!r}")

    return _normalise(Q, pi)


def P_matrix(Q: np.ndarray, t) -> np.ndarray:
    """Transition matrix at time t, vectorised over t."""
    t = np.atleast_1d(np.asarray(t, dtype=float))
    return np.stack([expm(Q * ti) for ti in t]) if t.size > 1 else expm(Q * float(t[0]))


def spectrum(Q: np.ndarray) -> np.ndarray:
    """Eigenvalues of Q, descending: 0 = l1 > l2 >= ... >= lK (reversible Q)."""
    w = np.linalg.eigvals(Q).real
    return np.sort(w)[::-1]


def stationary_from_Q(Q: np.ndarray) -> np.ndarray:
    w, V = np.linalg.eig(Q.T)
    idx = int(np.argmin(np.abs(w)))
    pi = np.abs(V[:, idx].real)
    return pi / pi.sum()
