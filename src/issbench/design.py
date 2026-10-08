"""Uniform-Iss sampling design (Part II, Theorem 11).

Fix a branch-length *pattern* d^0 and a scale s, so that every branch length is
multiplied by s.  With mu(s) = E[Iss] monotone, drawing

    p_s(s) = mu'(s) / Z,   Z = mu(s_max) - mu(s_min),   F(s) = (mu(s)-mu_min)/Z

makes Iss uniform on [mu_min, mu_max].  The normalising constant is the
increment of mu itself, so no pilot run is needed.

The closure below is exact for N=4 (mu from the pattern enumeration in
exact.py); the second-order affine form is available as `mu_affine` for
comparison, and the difference between the two marginals is itself a check of
Theorem 1's validity window.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .exact import E_Iss
from .gamma import discrete_gamma_rates, gamma_weights
from .trees import Node, quartet, star


@dataclass
class Design:
    shape: Callable[[float], Node]      # s -> tree
    N: int
    spec_Q: np.ndarray
    pi: np.ndarray
    s_min: float
    s_max: float
    rates: tuple
    weights: tuple
    convention: str = "B"

    def mu(self, s: float) -> float:
        return E_Iss(self.N, self.shape(s), self.spec_Q, self.pi,
                     rates=self.rates, weights=self.weights,
                     convention=self.convention)

    def mu_grid(self, n: int = 200):
        ss = np.exp(np.linspace(np.log(self.s_min), np.log(self.s_max), n))
        return ss, np.array([self.mu(s) for s in ss])

    def sample(self, n_draw: int, seed: int = 0) -> np.ndarray:
        """Draw s values whose induced Iss are uniform on [mu_min, mu_max]."""
        ss, mus = self.mu_grid(400)
        order = np.argsort(mus)
        mus, ss = mus[order], ss[order]
        mus, idx = np.unique(mus, return_index=True)
        ss = ss[idx]
        rng = np.random.default_rng(seed)
        u = rng.random(n_draw)
        target = mus[0] + u * (mus[-1] - mus[0])
        return np.interp(target, mus, ss)

    def density(self, s) -> np.ndarray:
        """p_s(s) = mu'(s)/Z, by central differences."""
        s = np.atleast_1d(np.asarray(s, dtype=float))
        h = np.maximum(1e-6, 1e-4 * s)
        d = np.array([(self.mu(a + b) - self.mu(a - b)) / (2 * b)
                      for a, b in zip(s, h)])
        Z = self.mu(self.s_max) - self.mu(self.s_min)
        return d / Z


def quartet_design(spec_Q, pi, internal_fraction: float, N: int = 4,
                   alpha: float | np.inf = np.inf, s_min: float = 0.05,
                   s_max: float = 60.0, convention: str = "B") -> Design:
    """Quartet with internal branch = internal_fraction * total tree length."""
    fi = float(internal_fraction)

    def shape(s: float) -> Node:
        r = fi * s
        l = (s - r) / 4.0
        return quartet(l, l, l, l, r)

    rates, weights = _rates(alpha)
    return Design(shape, N, spec_Q, np.asarray(pi, float), s_min, s_max,
                  rates, weights, convention)


def star_design(spec_Q, pi, N: int, alpha: float | np.inf = np.inf,
                s_min: float = 0.02, s_max: float = 40.0,
                convention: str = "B") -> Design:
    """Star tree; s is the root-to-leaf length (total tree length = N*s)."""
    rates, weights = _rates(alpha)
    return Design(lambda s: star(N, s / N), N, spec_Q, np.asarray(pi, float),
                  s_min, s_max, rates, weights, convention)


def _rates(alpha):
    if alpha is None or np.isinf(alpha):
        return (1.0,), (1.0,)
    return tuple(discrete_gamma_rates(alpha, 4)), tuple(gamma_weights(4))


def mu_affine(design: Design, s: float) -> float:
    """Second-order prediction (Part I Thm 1), for comparison with the exact mu.

        E[Iss] = 1 - [ (K-1)/N + (N-1)/N * Psi ] / (2 h(pi)),   Psi = 3 xi

    valid for JC69 with uniform pi (the case used in the designs here).
    """
    from math import log
    from .entropy import h
    from .trees import pairwise_distances
    D = pairwise_distances(design.shape(s))
    iu = np.triu_indices(design.N, 1)
    xi = float(np.mean([np.mean(np.exp(-4 * r * D[iu] / 3.0))
                        for r in design.rates]))
    K = len(design.pi)
    return 1.0 - ((K - 1) / design.N + (design.N - 1) / design.N * 3 * xi) \
        / (2 * h(design.pi, "nats"))
