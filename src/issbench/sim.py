"""Cell-level simulation and scoring.

One 'cell' = (tree, model, sequence length, number of replicates, seed).
AliSim is the only simulator; everything here is book-keeping around it.
Cells are addressed by a stable key so a killed run restarts where it stopped.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np

from .alisim import run_alisim
from .stats import (h_bar, iss, mean_p_distance, pi_hat, read_fasta,
                    site_counts)


def cell_key(**kw) -> str:
    """Deterministic, filesystem-safe key for a parameter set."""
    blob = json.dumps(kw, sort_keys=True, default=str)
    return hashlib.sha1(blob.encode()).hexdigest()[:16]


def simulate(tree_file: Path, model: str, length: int, n_reps: int,
             seed: int, workdir: Path, tag: str = "sim",
             threads: int = 1) -> list[Path]:
    """Run AliSim for one cell; return the alignment paths."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    prefix = workdir / tag
    outs = run_alisim(prefix, tree_file, model, length=int(length),
                      n_alignments=int(n_reps), seed=int(seed),
                      threads=threads, log=workdir / f"{tag}.log")
    # AliSim writes one file per replicate unless --single-output is honoured,
    # in which case the single file is a concatenation of all n_reps
    if len(outs) == 1 and n_reps > 1:
        outs = _split_single(outs[0], n_reps, workdir, tag)
    if len(outs) != n_reps:
        # never return a short list: the caller zips it against replicate
        # indices, and a short list silently mislabels every replicate after
        # the missing one
        raise RuntimeError(
            f"AliSim returned {len(outs)} file(s) for n_reps={n_reps} "
            f"(tag={tag}, dir={workdir}); refusing to record a partial cell")
    return outs


def _split_single(path: Path, n: int, workdir: Path, tag: str) -> list[Path]:
    """Split a single concatenated FASTA of n alignments into n files."""
    txt = Path(path).read_text(encoding="utf-8", errors="replace")
    blocks = [b for b in re.split(r"\n(?=>)", txt.strip()) if b.strip()]
    seqs = {}
    for b in blocks:
        lines = b.splitlines()
        name = lines[0][1:].split()[0]
        seq = "".join(lines[1:]).strip()
        seqs.setdefault(name, []).append(seq)
    names = list(seqs)
    out = []
    for i in range(n):
        f = workdir / f"{tag}_r{i}.fa"
        with f.open("w", encoding="utf-8") as fh:
            for nm in names:
                fh.write(f">{nm}\n{seqs[nm][i]}\n")
        out.append(f)
    # the concatenation is superseded by the per-replicate files; leaving it
    # around only creates something for a later glob to trip over
    Path(path).unlink(missing_ok=True)
    return out


def score(path: Path) -> dict:
    """All per-alignment statistics needed downstream."""
    names, M = read_fasta(path)
    q = pi_hat(M)
    return dict(n_seq=len(names), n_sites=M.shape[1],
                pbar=mean_p_distance(M),
                hbar=h_bar(M),
                iss_A=iss(M, "A"), iss_B=iss(M, "B"),
                pi_A=float(q[0]), pi_C=float(q[1]),
                pi_G=float(q[2]), pi_T=float(q[3]))


def score_cell(paths, extra: dict | None = None) -> list[dict]:
    rows = []
    for p in paths:
        r = score(p)
        r["file"] = str(p)
        if extra:
            r.update(extra)
        rows.append(r)
    return rows


def true_topology_split(a: float, b: float, c: float, d: float) -> str:
    """Canonical split string for the true quartet ((A,B),(C,D))."""
    return f"((A:{a:g},B:{b:g}),(C:{c:g},D:{d:g}))"
