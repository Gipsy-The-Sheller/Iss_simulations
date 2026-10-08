"""Statistics computed on alignments: Iss, mean p-distance, C-value, lambda."""
from __future__ import annotations

import numpy as np
from pathlib import Path
from typing import Sequence

from .entropy import h, h_from_counts, multinomial_entropy

CODES = {"A": 0, "C": 1, "G": 2, "T": 3,
         "a": 0, "c": 1, "g": 2, "t": 3}


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
    out = np.zeros((M.shape[1], K), dtype=np.int32)
    for j in range(M.shape[1]):
        out[j] = np.bincount(M[:, j], minlength=K)
    return out


def h_bar(M: np.ndarray, K: int = 4, base: str = "nats") -> float:
    """Mean per-site entropy H_bar."""
    C = site_counts(M, K)
    return float(np.mean([h_from_counts(c, base) for c in C])) if len(C) else float("nan")


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
    tot, cnt = 0.0, 0
    for i in range(n):
        for j in range(i + 1, n):
            tot += float(np.mean(M[i] != M[j]))
            cnt += 1
    return tot / cnt


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
