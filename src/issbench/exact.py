"""Exact per-site law of H_i, hence exact E[I_ss], sigma^2 and higher moments.

Two routes:
  * star_exact_*   -- closed form (Part I, Lemma 5 / Theorem 2), any N
  * tree_patterns  -- Felsenstein pruning + full pattern enumeration, N <= 8
                      (4^N patterns: 256 for N=4), arbitrary topology

This module is the ground truth against which AliSim output and the
second-order form (Theorem 1) are checked.
"""
from __future__ import annotations

import numpy as np
from itertools import product
from math import log
from typing import Dict, Tuple, Sequence

from .entropy import h, h_from_counts, multinomial_entropy
from .models import P_matrix
from .trees import Node


# ----------------------------------------------------------------- star trees
def star_site_law(N: int, q) -> Dict[Tuple[int, ...], float]:
    """Not needed: the star law is used only through its multinomial moment."""
    raise NotImplementedError("use multinomial_entropy / star_E_H instead")


def star_E_H(N: int, d: float, Q: np.ndarray, pi, rates, weights,
             base: str = "nats") -> float:
    """Exact E[h(n/N)] on a star tree with root-to-leaf distance d (Theorem 2)."""
    pi = np.asarray(pi, dtype=float)
    acc = 0.0
    for rc, wc in zip(rates, weights):
        Pm = P_matrix(Q, rc * d)
        for rho in range(len(pi)):
            acc += pi[rho] * wc * multinomial_entropy(N, Pm[rho], base)
    return float(acc)


def star_E_H2(N: int, d: float, Q: np.ndarray, pi, rates, weights) -> float:
    """Exact E[h^2] on a star tree (nats) -- for Var(H) and the design cost table."""
    pi = np.asarray(pi, dtype=float)
    acc = 0.0
    for rc, wc in zip(rates, weights):
        Pm = P_matrix(Q, rc * d)
        for rho in range(len(pi)):
            acc += pi[rho] * wc * _multinomial_h2(N, Pm[rho])
    return float(acc)


def _multinomial_h2(N: int, q) -> float:
    """E[h^2] for n ~ Multinomial(N,q), h = -sum p ln p, p = n/N."""
    from math import comb
    lnN = log(N)
    tot = 0.0
    K = len(q)
    rng = range(N + 1)
    # E[h^2] = ln^2 N - (2 lnN / N) sum_j E[n_j ln n_j]
    #          + (1/N^2) sum_j E[(n_j ln n_j)^2]
    #          + (1/N^2) sum_{j!=k} E[n_j ln n_j] E[n_k ln n_k]
    first = 0.0
    second = 0.0
    means = []
    for qj in q:
        s1 = 0.0
        s2 = 0.0
        for m in rng:
            if m == 0:
                continue
            pmf = comb(N, m) * (qj ** m) * ((1.0 - qj) ** (N - m))
            s1 += m * log(m) * pmf
            s2 += (m * log(m)) ** 2 * pmf
        first += s1
        second += s2
        means.append(s1)
    cross = sum(means) ** 2 - sum(m * m for m in means)
    return lnN ** 2 - 2.0 * lnN * first / N + (second + cross) / N ** 2


# ------------------------------------------------------- general tree topology
def tree_patterns(node: Node, Q: np.ndarray, pi, rates, weights) -> Dict[Tuple[int, ...], float]:
    """Full pattern distribution for a rooted binary tree with N <= 8 leaves.

    Felsenstein pruning, one pass per pattern, over all K^N patterns.
    Returns {state tuple in leaf order: probability}, marginalised over rates.
    """
    leaves = node.leaves()
    N = len(leaves)
    if N > 8:
        raise ValueError("pattern enumeration is exponential; use a star tree for large N")
    K = len(pi)
    pi = np.asarray(pi, dtype=float)
    pos = {id(l): i for i, l in enumerate(leaves)}

    def recurse(v: Node, pattern, Pm):
        if v.is_leaf:
            L = np.zeros(K)
            L[pattern[pos[id(v)]]] = 1.0
            return L
        L = np.ones(K)
        for c in v.children:
            L *= Pm[id(c)] @ recurse(c, pattern, Pm)
        return L

    out: Dict[Tuple[int, ...], float] = {}
    for rc, wc in zip(rates, weights):
        Pm = {id(c): P_matrix(Q, rc * c.length) for c in _all_children(node)}
        for pattern in product(range(K), repeat=N):
            p = float(np.dot(pi, recurse(node, pattern, Pm)))
            if p:
                out[pattern] = out.get(pattern, 0.0) + wc * p
    return out


def _all_children(node: Node):
    for c in node.children:
        yield c
        yield from _all_children(c)


# ------------------------------------------------------------------------ Iss
def equidistant_star(node):
    """Return the common root-to-leaf length if `node` is an equidistant star,
    else None.  A star has an exact closed form for any N, so it must never be
    routed to the (exponential) pattern enumeration."""
    if node.is_leaf or not node.children:
        return None
    if not all(c.is_leaf for c in node.children):
        return None
    d = node.children[0].length
    if all(abs(c.length - d) <= 1e-12 * max(1.0, abs(d)) for c in node.children):
        return float(d)
    return None


def E_Iss(N: int, tree_or_d, Q, pi, rates=(1.0,), weights=(1.0,),
          convention: str = "A", base: str = "nats") -> float:
    """Exact E[I_ss] -- star trees analytically (any N), other trees with N<=8
    by pattern enumeration."""
    pi = np.asarray(pi, dtype=float)
    if not isinstance(tree_or_d, (int, float)):
        n_leaf = len(tree_or_d.leaves())
        if n_leaf != N:
            raise ValueError(
                f"N={N} disagrees with the tree, which has {n_leaf} leaves; "
                "the reference entropy would be computed at the wrong N")
    if isinstance(tree_or_d, (int, float)):
        EH = star_E_H(N, float(tree_or_d), Q, pi, rates, weights, base)
    elif (d_star := equidistant_star(tree_or_d)) is not None:
        EH = star_E_H(N, d_star, Q, pi, rates, weights, base)
    else:
        pats = tree_patterns(tree_or_d, Q, pi, rates, weights)
        EH = 0.0
        for pat, pr in pats.items():
            cnt = np.bincount(list(pat), minlength=len(pi))
            EH += pr * h_from_counts(cnt, base)
    Href = multinomial_entropy(N, pi, base) if convention == "A" else h(pi, base)
    return EH / Href
