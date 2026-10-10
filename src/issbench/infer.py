"""Topology recovery for the four-taxon experiments.

Three criteria, all of which reduce to choosing one of three splits:
  ML   - delegated to IQ-TREE
  MPB  - maximum parsimony, delegated to MPBoot
  NJ   - for four taxa, neighbour joining == minimise the within-pair distance
         sum, evaluated on JC-corrected distances

The two searching criteria are run through the published programs rather than
reimplemented, so that no reported topology comes from code of ours that a
reader would have to take on trust.  NJ is the one in-process criterion: on a
quartet it is a deterministic function of six pairwise distances, with no
search and no heuristic.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
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
    """Maximum-parsimony tree from MPBoot (Hoang et al. 2018)."""
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


def nj_split(M: np.ndarray) -> str:
    """Neighbour joining. For four taxa this is the exact criterion: minimise
    the sum of within-pair JC-corrected distances over the three splits."""
    n = M.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = float(jc_distance(np.mean(M[i] != M[j])))
    sums = {name: D[o[0], o[1]] + D[o[2], o[3]] for name, o in SPLITS.items()}
    return min(sums, key=sums.get)


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
                s = {"NJ": nj_split}[method](M)
            cache[key] = s
        hits += int(s == TRUE_SPLIT)
    return hits / len(paths)
