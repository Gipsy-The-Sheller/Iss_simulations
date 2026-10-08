"""Topology recovery for the four-taxon experiments.

Three criteria, all of which reduce to choosing one of three splits:
  ML  - delegated to IQ-TREE (so no part of the loop is our own implementation)
  MP  - Fitch parsimony on the unrooted quartet, via a precomputed 4^4 table
  NJ  - for four taxa, neighbour-joining == minimise the within-pair distance
        sum; run both with JC-corrected and with uncorrected p-distances
        (the latter is the classical four-point rule)
"""
from __future__ import annotations

import re
import subprocess
from itertools import product
from pathlib import Path

import numpy as np

from .alisim import iqtree_path

LABELS = ("A", "B", "C", "D")
SPLITS = {"AB|CD": (0, 1, 2, 3), "AC|BD": (0, 2, 1, 3), "AD|BC": (0, 3, 1, 2)}
TRUE_SPLIT = "AB|CD"


# --------------------------------------------------------------------- ML
def iqtree_ml(msa: Path, model: str = "JC69", seed: int = 1,
              workdir: Path | None = None, timeout: int = 900) -> str | None:
    """Return the ML tree as Newick, or None on failure."""
    msa = Path(msa)
    workdir = Path(workdir or msa.parent)
    prefix = workdir / (msa.stem + ".ml")
    cmd = [iqtree_path(), "-s", str(msa), "-m", model, "-pre", str(prefix),
           "--seed", str(seed), "-T", "1", "--quiet", "-redo"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=timeout, cwd=str(workdir))
    except subprocess.TimeoutExpired:
        return None
    tf = Path(str(prefix) + ".treefile")
    if r.returncode != 0 or not tf.exists():
        return None
    return tf.read_text(encoding="utf-8").strip()


def split_of_newick(nwk: str | None) -> str | None:
    """Which quartet split does an unrooted 4-taxon Newick imply?"""
    if not nwk:
        return None
    t = nwk.split("\n")[0].strip()
    t = re.sub(r":[\d.eE+-]+", "", t)               # drop branch lengths
    t = re.sub(r"\)[\d.eE+-]+(?=[,);])", ")", t)    # drop internal support
    for name, (i, j, k, l) in SPLITS.items():
        for pair in (((i, j), (k, l)), ((k, l), (i, j))):
            (a, b), (c, d) = pair
            pat = (r"\(\s*\(?[^()]*\b%s\b[^()]*\b%s\b[^()]*\)?\s*,\s*"
                   r"\(?[^()]*\b%s\b[^()]*\b%s\b[^()]*\)?\s*\)"
                   % (LABELS[a], LABELS[b], LABELS[c], LABELS[d]))
            if re.search(pat, t):
                return name
    return None


# --------------------------------------------------------------------- MP
def _build_fitch_table() -> np.ndarray:
    """cost[(a,b,c,d)] = min substitutions on an unrooted quartet."""
    tbl = np.zeros((4, 4, 4, 4), dtype=np.int16)
    for a, b, c, d in product(range(4), repeat=4):
        best = 99
        for x, y in product(range(4), repeat=2):
            best = min(best, (a != x) + (b != x) + (c != y) + (d != y) + (x != y))
        tbl[a, b, c, d] = best
    return tbl


_FITCH = _build_fitch_table()


def mp_split(M: np.ndarray) -> str:
    """Fitch parsimony; for 4 taxa the tree cost is the sum over the other split
    of the site-wise quartet costs, minimised over the three topologies."""
    cost = {}
    for name, (i, j, k, l) in SPLITS.items():
        cost[name] = int(_FITCH[M[i], M[j], M[k], M[l]].sum())
    return min(cost, key=cost.get)


# --------------------------------------------------------------------- NJ
def jc_distance(p):
    p = np.asarray(p, dtype=float)
    out = np.where(p >= 0.75, np.inf, -0.75 * np.log(np.maximum(1.0 - 4.0 * p / 3.0, 1e-300)))
    return out


def nx_split(M: np.ndarray, corrected: bool = True) -> str:
    """Minimise the sum of within-pair distances (== NJ for four taxa)."""
    n = M.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            p = float(np.mean(M[i] != M[j]))
            D[i, j] = D[j, i] = float(jc_distance(p)) if corrected else p
    sums = {name: D[o[0], o[1]] + D[o[2], o[3]] for name, o in SPLITS.items()}
    return min(sums, key=sums.get)


def nj_split(M: np.ndarray) -> str:
    return nx_split(M, corrected=True)


def four_point_split(M: np.ndarray) -> str:
    return nx_split(M, corrected=False)


def ptrue(paths, method: str, model: str = "JC69", seed0: int = 1,
          workdir: Path | None = None, cache: dict | None = None) -> float:
    """Fraction of replicates whose inferred split is the true one."""
    from .stats import read_fasta
    cache = cache if cache is not None else {}
    hits = 0
    for k, f in enumerate(paths):
        key = (str(f), method)
        if key in cache:
            s = cache[key]
        else:
            if method == "ML":
                s = split_of_newick(iqtree_ml(f, model, seed0 + k, workdir))
            else:
                M = read_fasta(f)[1]
                s = {"MP": mp_split, "NJ": nj_split,
                     "4P": four_point_split}[method](M)
            cache[key] = s
        hits += int(s == TRUE_SPLIT)
    return hits / len(paths)
