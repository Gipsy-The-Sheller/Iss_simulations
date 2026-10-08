"""Discrete-Gamma rate categories and the continuous Gamma Laplace transform.

Yang (1994) J Mol Evol 39:306-314, doi:10.1007/BF00160154 defines both an
equal-probability *mean* discretisation (the standard one, used by PAML,
IQ-TREE, PhyML, RAxML, DAMBE) and a *median* variant.  They are NOT
interchangeable: for E[e^{-4tr/3}] they differ by ~2-5% (mean) vs ~3-11%
(median) relative to the continuous limit.  See docs/Hfss口径核对 and the
audit report, item A1.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import gamma as _gd
from scipy.special import gammainc


def discrete_gamma_rates(alpha: float, k: int = 4, method: str = "mean"):
    """Equal-probability discrete Gamma, normalised to mean 1.

    method="mean"   : category means (Yang 1994; the standard)
    method="median" : category quantile midpoints (non-standard)
    """
    if np.isinf(alpha):
        return np.arange(1, k + 1, dtype=float)
    scale = 1.0 / alpha
    if method == "median":
        r = np.array([_gd.ppf((2 * i + 1) / (2 * k), alpha, scale=scale)
                      for i in range(k)])
    elif method == "mean":
        bounds = [0.0] + [_gd.ppf(i / k, alpha, scale=scale) for i in range(1, k)] + [np.inf]
        f = lambda x: gammainc(alpha + 1.0, alpha * x)          # noqa: E731
        r = np.array([scale * (f(bounds[i + 1]) - f(bounds[i])) * k for i in range(k)])
    else:
        raise ValueError(method)
    return r / r.mean()


def gamma_weights(k: int = 4) -> np.ndarray:
    return np.full(k, 1.0 / k)


def laplace_gamma(alpha: float, lam: float, t: float) -> float:
    """E_r[e^{lam * r * t}] for r ~ Gamma(alpha, scale=1/alpha).

    Part I, Lemma 6: (1 - lam*t/alpha)^{-alpha}, valid for lam < 0.
    """
    if np.isinf(alpha):
        return float(np.exp(lam * t))
    return float((1.0 - lam * t / alpha) ** (-alpha))


def laplace_discrete(alpha: float, lam: float, t: float, k: int = 4,
                     method: str = "mean") -> float:
    r = discrete_gamma_rates(alpha, k, method)
    return float(np.mean(np.exp(lam * r * t)))
