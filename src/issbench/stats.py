"""Statistics computed on alignments: Iss, mean p-distance, C-value, lambda."""
from __future__ import annotations

import numpy as np
from pathlib import Path
from typing import Sequence

from .entropy import h, h_from_counts, multinomial_entropy

CODES = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3,
         "a": 0, "c": 1, "g": 2, "t": 3, "u": 3}


def read_fasta(path: str | Path) -> tuple[list[str], np.ndarray]:
    """Read a FASTA alignment -> (names, matrix of shape (n_seq, n_sites))."""
    names, seqs, cur = [], [], []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
                cur = []
            names.append(line[1:].split()[0])
        else:
            cur.append(line)
    if cur:
        seqs.append("".join(cur))
    n = min(len(s) for s in seqs)
    M = np.full((len(seqs), n), -1, dtype=np.int8)
    for i, s in enumerate(seqs):
        for j, ch in enumerate(s[:n]):
            M[i, j] = CODES.get(ch, -1)
    return names, M


def site_counts(M: np.ndarray, K: int = 4) -> np.ndarray:
    """(n_sites, K) count matrix; sites containing ambiguity codes are dropped."""
    bad = (M < 0).any(axis=0)
    M = M[:, ~bad]
    return np.stack([(M == k).sum(axis=0) for k in range(K)], axis=1).astype(np.int32)


def _entropy_matrix(C: np.ndarray, base: str = "nats") -> np.ndarray:
    """Per-row Shannon entropy of a count matrix, vectorised."""
    P = C / C.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        T = np.where(P > 0, P * np.log(P), 0.0)
    if base == "bits":
        T = T / np.log(2.0)
    return -T.sum(axis=1)


def h_bar(M: np.ndarray, K: int = 4, base: str = "nats") -> float:
    """Mean per-site entropy H_bar (vectorised over sites)."""
    C = site_counts(M, K)
    if len(C) == 0:
        return float("nan")
    return float(_entropy_matrix(C, base).mean())


def pi_hat(M: np.ndarray, K: int = 4) -> np.ndarray:
    """Pooled empirical base frequencies over all (sequence, site) observations."""
    C = site_counts(M, K)
    tot = C.sum(axis=0).astype(float)
    return tot / tot.sum()


def iss(M: np.ndarray, convention: str = "A", base: str = "nats",
        K: int = 4) -> float:
    """Xia's index of substitution saturation.

    convention "A": H_FSS = E_Mult(N, pi_hat)[h]   (Xia et al. 2003, Eq. 3)
    convention "B": H_FSS = h(pi_hat)
    """
    N = M.shape[0]
    q = pi_hat(M, K)
    Href = (multinomial_entropy(N, q, base) if convention == "A" else h(q, base))
    return h_bar(M, K, base) / Href


def mean_p_distance(M: np.ndarray) -> float:
    """Mean over unordered sequence pairs of the uncorrected p-distance."""
    n, L = M.shape
    if n < 2:
        return float("nan")
    iu = np.triu_indices(n, k=1)
    diffs = (M[iu[0]] != M[iu[1]]).mean(axis=1)
    return float(diffs.mean())


def c_value(M: np.ndarray) -> float:
    """Struck et al. (2008) C factor.  DEFINITION TO BE PINNED DOWN.

    Ref: Struck TH et al. (2008) Mol Phylogenet Evol 48:628-645,
         doi:10.1016/j.ympev.2008.05.015
    Implemented per the original: SD of the transition/transversion ratio
    across sequence pairs, divided by the SD of the pairwise distance.
    Verify against the paper before use.
    """
    raise NotImplementedError("pin the definition against Struck et al. 2008 first")


def lam(M: np.ndarray) -> float:
    """Jermiin & Misof (2020) lambda (loss of historical signal).

    Ref: Jermiin LS & Misof B (2020) bioRxiv 2020.01.03.894097 (preprint only)
    DEFINITION TO BE PINNED DOWN against the preprint before use.
    """
    raise NotImplementedError("pin the definition against Jermiin & Misof 2020 first")


def parsimony_informative_fraction(M: np.ndarray, K: int = 4) -> float:
    """Fraction of sites that are parsimony-informative.

    For K=4 the informative pattern is exactly 2:2 (two states, two sequences
    each): 2:1:1 is informative for some tree shapes but not unrooted-quartet
    informative, and 3:1 / 4:0 never are.  Used as the target Y in the
    conditional-independence test.
    """
    C = site_counts(M, K)
    if len(C) == 0:
        return float("nan")
    n_seq = M.shape[0]
    two_two = (C.max(axis=1) == n_seq // 2) & ((C > 0).sum(axis=1) == 2)
    return float(np.mean(two_two))
