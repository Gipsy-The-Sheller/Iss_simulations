"""The two other tree-free saturation indices, on the same alignments.

Both are pairwise constructions and neither is a function of the per-site
count vector, so unlike the index itself they need the alignment retained
rather than summarised:

  C factor (Struck et al. 2008)   C = sigma(ti/tv) / sigma(p), the standard
      deviations taken over all pairwise comparisons.  As saturation grows the
      ti/tv values converge on 1, so their spread falls, while the spread of
      the pairwise distances rises; C therefore falls with saturation.  There
      is no intrinsic cut-off.

  lambda (Jermiin & Misof 2020)   lambda(ij) = delta_obs(ij) / delta_ran(ij),
      with delta_obs the p-distance between i and j and delta_ran the value
      expected if the two had evolved independently to equilibrium.  lambda = 0
      is a fully retained historical signal and lambda = 1 a fully eroded one.
      It is defined per pair; we report the mean over pairs.

The avian alignment is read twice: once to accumulate the per-site counts and
so determine which columns are complete, and once to extract those columns.
Storing the whole 51 x 322M alignment would be 16 GB; the complete columns are
51 x 32M = 1.6 GB.

    python scripts/s11_alt_indices.py --avian
    python scripts/s11_alt_indices.py --coi
"""
from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
OUT = ROOT / "results" / "empirical"
GZ = ROOT / "work" / "empirical" / "ALL.aln.gz"

_LUT = np.full(256, 4, dtype=np.uint8)
for _ws in (9, 10, 13, 32):
    _LUT[_ws] = 5
for _ch, _k in (("A", 0), ("C", 1), ("G", 2), ("T", 3)):
    _LUT[ord(_ch)] = _k
    _LUT[ord(_ch.lower())] = _k

# transitions: A<->G, C<->T
_TI = {(0, 2), (2, 0), (1, 3), (3, 1)}


def pairwise(M: np.ndarray):
    """Per-pair p-distance, ti/tv and lambda for an (n_tax, n_sites) uint8 matrix."""
    n, L = M.shape
    comp = np.stack([(M == k).sum(axis=1) for k in range(4)]) / L   # (4, n)
    joint = np.zeros((16,), dtype=np.int64)
    rows = []
    for i in range(n):
        for j in range(i + 1, n):
            key = M[i].astype(np.uint16) * 4 + M[j]
            c = np.bincount(key, minlength=16).astype(np.float64)
            same = c[0] + c[5] + c[10] + c[15]
            p = (L - same) / L
            ti = sum(c[a * 4 + b] for a, b in _TI)
            tv = (L - same) - ti
            d_ran = 1.0 - float((comp[:, i] * comp[:, j]).sum())
            rows.append(dict(i=i, j=j, p=p,
                             ti_tv=(ti / tv) if tv > 0 else np.nan,
                             delta_ran=d_ran,
                             lam=(p / d_ran) if d_ran > 0 else np.nan))
    return pd.DataFrame(rows)


def summarise(df: pd.DataFrame, label: str, n_tax: int, n_sites: int) -> dict:
    return dict(dataset=label, n_tax=n_tax, n_sites=n_sites, n_pairs=len(df),
                C=float(df.ti_tv.std() / df.p.std()),
                lam_mean=float(df.lam.mean()), lam_median=float(df.lam.median()),
                lam_lo=float(df.lam.quantile(0.025)),
                lam_hi=float(df.lam.quantile(0.975)),
                ti_tv_mean=float(df.ti_tv.mean()),
                p_mean=float(df.p.mean()))


def avian(limit_sites: int = 0) -> None:
    t0 = time.time()
    counts = None
    n_tax = n_sites = 0
    M_blocks = []
    # pass 1: counts -> complete columns
    with gzip.open(GZ, "rb") as fh:
        head = fh.readline().decode("ascii", "replace").split()
        n_tax, n_sites = int(head[0]), int(head[1])
        counts = np.zeros((5, n_sites), dtype=np.uint8)
        for i in range(n_tax):
            seq = fh.readline().strip().split(None, 1)[1]
            cat = _LUT[np.frombuffer(seq, dtype=np.uint8)]
            for k in range(5):
                counts[k] += (cat == k)
            if (i + 1) % 10 == 0:
                print(f"  pass1 {i + 1}/{n_tax}  {time.time() - t0:.0f}s",
                      flush=True)
    valid = counts[4] == 0
    n_valid = int(valid.sum())
    print(f"  complete columns: {n_valid:,} of {n_sites:,}", flush=True)
    del counts

    # pass 2: extract the complete columns
    M = np.empty((n_tax, n_valid), dtype=np.uint8)
    with gzip.open(GZ, "rb") as fh:
        fh.readline()
        for i in range(n_tax):
            seq = fh.readline().strip().split(None, 1)[1]
            cat = _LUT[np.frombuffer(seq, dtype=np.uint8)]
            M[i] = cat[valid]
            if (i + 1) % 10 == 0:
                print(f"  pass2 {i + 1}/{n_tax}  {time.time() - t0:.0f}s",
                      flush=True)
    del valid
    print(f"  matrix {M.shape} ({M.nbytes / 2**30:.2f} GB), "
          f"{time.time() - t0:.0f}s", flush=True)

    if limit_sites:
        M = M[:, :limit_sites]

    def hbar_of(X: np.ndarray) -> tuple[float, np.ndarray]:
        """Mean per-site entropy, and the pooled composition, in blocks."""
        acc, seen = 0.0, 0
        for s in range(0, X.shape[1], 4_000_000):
            # axis=0: sum over TAXA, giving per-SITE counts.  Summing over
            # axis=1 would give per-taxon compositions, whose mean entropy is
            # a different (and near-constant) quantity.
            c = np.stack([(X[:, s:s + 4_000_000] == k).sum(axis=0)
                          for k in range(4)]).astype(np.float32)
            P = c / np.maximum(c.sum(axis=0, keepdims=True), 1)
            acc += float(-(np.where(P > 0, P * np.log(P), 0.0)).sum())
            seen += P.shape[1]
        c = np.stack([(X == k).sum(axis=0) for k in range(4)]).astype(float)
        return acc / seen, c.sum(axis=1) / c.sum()

    rows, pairs = [], []
    subsets = {"all complete columns": M}
    # Duchene et al.'s modification: restrict the entropy to parsimony-
    # informative sites.  A site qualifies when two or more states each occur
    # in two or more taxa.
    pi_mask = np.zeros(M.shape[1], dtype=bool)
    for s in range(0, M.shape[1], 4_000_000):
        blk = M[:, s:s + 4_000_000]
        cc = np.stack([(blk == k).sum(axis=0) for k in range(4)])
        pi_mask[s:s + 4_000_000] = (cc >= 2).sum(axis=0) >= 2
    print(f"  parsimony-informative: {int(pi_mask.sum()):,} sites "
          f"({100 * pi_mask.mean():.1f}%)", flush=True)
    subsets["parsimony-informative sites"] = M[:, pi_mask]

    from issbench.entropy import h, multinomial_entropy
    for label, X in subsets.items():
        df = pairwise(X)
        row = summarise(df, f"avian, {label}", n_tax, X.shape[1])
        Eh, q = hbar_of(X)
        row.update(hbar=Eh, hfss_A=multinomial_entropy(n_tax, q),
                   hfss_B=h(q), iss_A=Eh / multinomial_entropy(n_tax, q),
                   iss_B=Eh / h(q))
        rows.append(row)
        df["subset"] = label
        pairs.append(df)
        print(f"  {label}: {X.shape[1]:,} sites  Iss(A)={row['iss_A']:.4f}  "
              f"C={row['C']:.1f}  lambda={row['lam_mean']:.4f}", flush=True)

    pd.concat(pairs, ignore_index=True).to_csv(OUT / "avian_pairs.csv",
                                               index=False)
    pd.DataFrame(rows).to_csv(OUT / "avian_alt_indices.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    print(f"total {time.time() - t0:.0f}s")


def coi() -> None:
    from issbench.stats import read_fasta
    names, A = read_fasta(ROOT / "data" / "VertCOI.FAS")
    PARTS = {"all": (0, 1, 2), "pos1": (0,), "pos2": (1,), "pos3": (2,),
             "pos1+2": (0, 1)}
    rows = []
    for tag, pos in PARTS.items():
        M = np.concatenate([A[:, p::3] for p in pos], axis=1).astype(np.uint8)
        df = pairwise(M)
        row = summarise(df, tag, M.shape[0], M.shape[1])
        c = np.stack([(M == k).sum(axis=0) for k in range(4)]).astype(float)
        q = c.sum(axis=1) / c.sum()
        from issbench.entropy import h, multinomial_entropy
        P = c / c.sum(axis=0, keepdims=True)
        Eh = float(-(np.where(P > 0, P * np.log(P), 0.0)).sum(axis=0).mean())
        row.update(hbar=Eh, hfss_A=multinomial_entropy(M.shape[0], q),
                   hfss_B=h(q), iss_A=Eh / multinomial_entropy(M.shape[0], q),
                   iss_B=Eh / h(q))
        rows.append(row)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "coi_alt_indices.csv", index=False)
    print(d.to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--avian", action="store_true")
    ap.add_argument("--coi", action="store_true")
    ap.add_argument("--sites", type=int, default=0)
    a = ap.parse_args()
    if a.avian:
        avian(a.sites)
    if a.coi:
        coi()
    if not (a.avian or a.coi):
        ap.error("pass --avian and/or --coi")
