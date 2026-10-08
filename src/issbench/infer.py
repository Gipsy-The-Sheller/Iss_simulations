"""Topology recovery for the four-taxon experiments.

Three criteria, all of which reduce to choosing one of three splits:
  ML  - delegated to IQ-TREE (so no part of the loop is our own implementation)
  MP  - Fitch parsimony on the unrooted quartet, via a precomputed 4^4 table
  NJ  - for four taxa, neighbour-joining == minimise the within-pair distance
        sum; run both with JC-corrected and with uncorrected p-distances
        (the latter is the classical four-point rule)
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
from itertools import product
from pathlib import Path

import numpy as np

from .alisim import iqtree_path

MPBOOT_DEFAULT = os.environ.get(
    "MPBOOT",
    "E:/Codes/YRTools/plugins/YR-MPE/YR_MPE/bin/mpboot/mpboot.exe")

LABELS = ("A", "B", "C", "D")
SPLITS = {"AB|CD": (0, 1, 2, 3), "AC|BD": (0, 2, 1, 3), "AD|BC": (0, 3, 1, 2)}
TRUE_SPLIT = "AB|CD"


# --------------------------------------------------------------------- ML
def iqtree_ml(msa: Path, model: str = "JC69", seed: int = 1,
              workdir: Path | None = None, timeout: int = 900) -> str | None:
    """Return the ML tree as Newick, or None on failure."""
    msa = Path(msa).resolve()
    workdir = Path(workdir or msa.parent).resolve()
    stem = msa.name
    prefix = msa.stem + ".ml"
    cmd = [iqtree_path(), "-s", stem, "-m", model, "-pre", prefix,
           "--seed", str(seed), "-T", "1", "--quiet", "-redo"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=timeout, cwd=str(workdir))
    except subprocess.TimeoutExpired:
        return None
    tf = workdir / (prefix + ".treefile")
    if r.returncode != 0 or not tf.exists():
        return None
    return tf.read_text(encoding="utf-8").strip()


def mpboot_path() -> str:
    import os
    p = os.environ.get("MPBOOT", MPBOOT_DEFAULT)
    if Path(p).exists():
        return p
    found = shutil.which("mpboot")
    if not found:
        raise FileNotFoundError(f"MPBoot not found; set $MPBOOT (tried {p!r})")
    return found


def mpboot_mp(msa: Path, seed: int = 1, workdir: Path | None = None,
              timeout: int = 900) -> str | None:
    """Maximum-parsimony tree from MPBoot (Hoang et al. 2018).

    Preferred over the in-process Fitch implementation for anything reported:
    MPBoot is the standard tool and its ratchet/SPR search explores more than a
    four-taxon enumeration can.  The in-process version is kept as `mp_split`
    for a free consistency check.
    """
    msa = Path(msa).resolve()
    workdir = Path(workdir or msa.parent).resolve()
    # MPBoot resolves -s and -pre relative to the CWD, and we already set cwd to
    # workdir; passing a path that also contains workdir doubles it up and the
    # run dies with "Cannot read file".  Bare names only, and cwd absolute.
    stem = msa.name
    prefix = msa.stem + ".mpb"
    cmd = [mpboot_path(), "-s", stem, "-pre", prefix,
           "-seed", str(seed), "-st", "DNA"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=timeout, cwd=str(workdir))
    except subprocess.TimeoutExpired:
        return None
    if r.returncode != 0:
        return None
    for cand in (workdir / (prefix + ".treefile"), workdir / (prefix + ".tree")):
        if cand.exists():
            return cand.read_text(encoding="utf-8").strip()
    return None


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
JC_P_CAP = 0.74          # p >= 3/4 is undefined for JC69; cap rather than inf


def jc_distance(p):
    """Jukes-Cantor corrected distance, capped at p = JC_P_CAP.

    Returning inf for saturated pairs (the previous behaviour) makes every
    split sum infinite and the minimisation arbitrary, which silently destroys
    the NJ arm near saturation.  Capping at p = 0.74 is the usual convention;
    it is a real limitation of distance methods in this regime and is reported
    as such rather than hidden.
    """
    p = np.asarray(p, dtype=float)
    pc = np.minimum(p, JC_P_CAP)
    return -0.75 * np.log(1.0 - 4.0 * pc / 3.0)


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
            elif method == "MPB":
                s = split_of_newick(mpboot_mp(f, seed0 + k, workdir))
            else:
                M = read_fasta(f)[1]
                s = {"MP": mp_split, "NJ": nj_split,
                     "4P": four_point_split}[method](M)
            cache[key] = s
        hits += int(s == TRUE_SPLIT)
    return hits / len(paths)
