"""Entropy primitives and the two H_FSS conventions.

Convention A is Xia et al. (2003) Eq. (3)  -- the *expected* entropy of a
multinomial sample of size N drawn with the pooled empirical frequencies:

    H_FSS^A = E_{n ~ Multinomial(N, pi_hat)} [ h(n/N) ]

Convention B is the plug-in composition entropy h(pi_hat), which is what the
Part I / Part II proofs use.

The two are related exactly, pointwise in the data, by

    I^A = c * I^B ,    c = h(pi_hat) / H_FSS^A(N, pi_hat)

so the affine-identity *criterion* is convention-independent; only the slope and
intercept scale by c.  See docs/Hfss口径核对_Xia2003.md.
"""
from __future__ import annotations

import numpy as np
from math import comb, log, log2
from typing import Literal

Convention = Literal["A", "B"]

LOG = {"nats": np.log, "bits": log2}


def h(p, base: str = "nats") -> float:
    """Shannon entropy of a probability vector; zero-probability states skipped."""
    p = np.asarray(p, dtype=float)
    p = p[p > 0]
    if base == "nats":
        return float(-np.sum(p * np.log(p)))
    if base == "bits":
        return float(-np.sum(p * log2(float(p))))
    raise ValueError(base)


def h_from_counts(n, base: str = "nats") -> float:
    """Site entropy from a count vector n_j (sum_j n_j = N)."""
    n = np.asarray(n, dtype=float)
    N = n.sum()
    if N <= 0:
        return 0.0
    p = n / N
    p = p[p > 0]
    if base == "nats":
        return float(-np.sum(p * np.log(p)))
    return float(-np.sum(p * log2(float(p))))


def multinomial_entropy(N: int, q, base: str = "nats") -> float:
    """Exact E[h(n/N)] for n ~ Multinomial(N, q).

    Lemma 5 of Part I:
        E[h] = ln N - (1/N) sum_j sum_{m=1}^{N} m ln m C(N,m) q_j^m (1-q_j)^{N-m}

    This is *exactly* Xia's H_FSS (convention A) when q = pi_hat.
    """
    q = np.asarray(q, dtype=float)
    lnN = log(N) if base == "nats" else log2(N)
    acc = 0.0
    for qj in q:
        if qj <= 0.0:
            continue
        for m in range(1, N + 1):
            w = m * log(m) * comb(N, m) * (qj ** m) * ((1.0 - qj) ** (N - m))
            acc += w
    val = lnN - acc / N
    return float(val) if base == "nats" else float(val / log(2.0))


def H_FSS(q, N: int, convention: Convention = "A", base: str = "nats") -> float:
    """Xia's saturation reference entropy under either convention.

    convention="A" -> E_{Mult(N,q)}[h]   (Xia et al. 2003, Eq. 3, base 2)
    convention="B" -> h(q)               (plug-in composition entropy)
    """
    if convention == "A":
        return multinomial_entropy(N, q, base)
    return h(q, base)


def miller_basharin_gap(N: int, q, base: str = "nats") -> float:
    """h(q) - H_FSS^A(N,q), i.e. the finite-N bias removed by convention A."""
    return h(q, base) - multinomial_entropy(N, q, base)


def convention_ratio(N: int, q, base: str = "nats") -> float:
    """c = h(q) / H_FSS^A(N,q)  s.t.  I^A = c * I^B."""
    return h(q, base) / multinomial_entropy(N, q, base)
