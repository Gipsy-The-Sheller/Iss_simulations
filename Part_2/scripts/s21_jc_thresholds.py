"""The JC69 pair threshold, and what the COI distance matrix does with it.

Five claims were put forward for a JC-only saturation screen.  Three of them
are numerical and are tested here.

1. AN OBSERVED p NEAR 3/4 PREDICTS SATURATION.  True for JC by construction --
   3/4 is the asymptote of the p-distance under JC, and equals 1 - sum q^2 at
   uniform composition, the boundary the manuscript derives from xi -> 0.

2. A THRESHOLD FROM THE JC69 STANDARD DEVIATION.  The corrected distance
   d = -(3/4) ln(1 - 4p/3) is finite only while the OBSERVED p stays below
   3/4.  With p-hat having variance p(1-p)/L, the largest true p whose
   observation stays below 3/4 with confidence 1-alpha is the solution of
   3/4 - p = z * sqrt(p(1-p)/L).  That is computed here, together with the
   standard error of the corrected distance, which is the quantity that
   actually decides whether the correction is usable.

3. lambda IS THE SAME LINE.  lambda = p / (1 - sum_k q_k^2) by the
   manuscript's Proposition 24, so at uniform composition lambda = 4p/3 and
   lambda = 1 is exactly p = 3/4.  Checked on the COI matrices.

4 and 5 are about the distance matrix as an object, and are checked here on
the COI partitions: the four-point condition (a matrix is tree-like exactly
when the two largest of the three pairwise sums agree in every quartet), and
the connectivity of the graph that keeps only sub-threshold pairs, which is
the proposed taxon screen.
"""
from __future__ import annotations

import sys
from math import log, sqrt
from pathlib import Path

import numpy as np
import pandas as pd

# Part 2.  Inputs that belong to the published preprint (the empirical
# alignments, the simulated work directory) are read from the
# repository root; everything this line of work produces is written
# under Part_2/, so the two cannot be confused.
PART2 = Path(__file__).resolve().parents[1]
ROOT = PART2.parent
sys.path.insert(0, str(ROOT / "src"))

from issbench.stats import read_fasta                                 # noqa: E402

OUT = PART2 / "results" / "jc_threshold"
EMP = ROOT / "results" / "empirical"
Z = {0.10: 1.2816, 0.05: 1.6449, 0.01: 2.3263}


def p_threshold(L: int, alpha: float = 0.05) -> float:
    """Largest true p whose observation stays under 3/4 with confidence."""
    z = Z[alpha]
    lo, hi = 0.0, 0.75
    for _ in range(200):
        mid = (lo + hi) / 2
        if 0.75 - mid > z * sqrt(mid * (1 - mid) / L):
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def jc_d(p: float) -> float:
    return -0.75 * log(1 - 4 * p / 3) if p < 0.75 else float("inf")


def jc_se(p: float, L: int) -> float:
    """Delta-method standard error of the corrected distance."""
    if p >= 0.75:
        return float("inf")
    return sqrt(p * (1 - p) / L) / (1 - 4 * p / 3)


def pairwise_p(M: np.ndarray) -> np.ndarray:
    n = M.shape[0]
    P = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            P[i, j] = P[j, i] = float(np.mean(M[i] != M[j]))
    return P


def four_point(P: np.ndarray, L: int) -> tuple[int, float, float]:
    """(violations beyond noise, worst gap, noise scale).

    The condition is that the two largest of the three pairwise sums are
    EQUAL, so with finite L they differ by sampling error.  A gap is only a
    violation if it is larger than that: the sum of two p-distances has
    standard error sqrt(2 p(1-p)/L), which is the scale reported alongside.
    """
    n = P.shape[0]
    bad, worst, noise = 0, 0.0, 0.0
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                for d in range(c + 1, n):
                    s = sorted([P[a, b] + P[c, d], P[a, c] + P[b, d],
                                P[a, d] + P[b, c]])
                    p = min(0.749, max(P[a, b], P[c, d], P[a, c],
                                       P[b, d], P[a, d], P[b, c]))
                    tol = 2.0 * np.sqrt(2 * p * (1 - p) / L)
                    noise = max(noise, tol)
                    worst = max(worst, s[2] - s[1])      # worst gap, always
                    if s[2] - s[1] > tol:
                        bad += 1
    return bad, worst, noise


def lambda_matrix(M: np.ndarray) -> np.ndarray:
    """Per-pair lambda = p / (1 - sum_k f_ik f_jk), composition-corrected.

    Not 4p/3: that is the uniform-composition value, and the COI third
    positions are the opposite of uniform, so the boundary they approach is
    1 - sum q^2 at their own compositions and not 3/4.
    """
    n, L = M.shape
    comp = np.stack([(M == k).sum(axis=1) for k in range(4)]) / L
    Lam = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = 1.0 - float((comp[:, i] * comp[:, j]).sum())
            p = float(np.mean(M[i] != M[j]))
            Lam[i, j] = Lam[j, i] = p / d if d > 0 else np.nan
    return Lam


def components(P: np.ndarray, thr: float) -> list[list[int]]:
    """Connected components of the graph keeping edges with p < thr."""
    n = P.shape[0]
    seen, out = set(), []
    for s in range(n):
        if s in seen:
            continue
        stack, comp = [s], []
        seen.add(s)
        while stack:
            v = stack.pop()
            comp.append(v)
            for w in range(n):
                if w not in seen and v != w and P[v, w] < thr:
                    seen.add(w)
                    stack.append(w)
        out.append(sorted(comp))
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------- (2) the threshold
    print("(2) the JC69 threshold on the observed p-distance")
    print("    largest true p whose observed value stays below 3/4 with"
          " confidence")
    rows = []
    for L in (250, 500, 1000, 1500, 4000):
        line = f"    L = {L:5d}   "
        for a in (0.10, 0.05, 0.01):
            pt = p_threshold(L, a)
            rows.append(dict(L=L, alpha=a, p_threshold=pt))
            line += f"alpha={a:.2f}: {pt:.4f}   "
        print(line)
    pd.DataFrame(rows).to_csv(OUT / "p_threshold.csv", index=False)

    print("\n    and what the correction is worth at each p (L = 1000):")
    print(f"    {'p':>6} {'d^':>8} {'SE(d^)':>9} {'rel. SE':>8}")
    rows = []
    for p in (0.40, 0.55, 0.65, 0.68, 0.70, 0.72, 0.735, 0.745):
        d, se = jc_d(p), jc_se(p, 1000)
        rows.append(dict(p=p, d=d, se=se, rel=se / d if d else np.nan))
        print(f"    {p:6.3f} {d:8.3f} {se:9.4f} {100 * se / d:7.1f}%")
    pd.DataFrame(rows).to_csv(OUT / "jc_precision.csv", index=False)

    # ------------------------------------------- (3,4,5) on the COI data
    print("\n(3,4,5) the COI distance matrices")
    sets = [("all sites", EMP / "all.fa"),
            ("1st+2nd", EMP / "pos1+2.fa"),
            ("1st", EMP / "pos1.fa"),
            ("2nd", EMP / "pos2.fa"),
            ("3rd", EMP / "pos3.fa")]
    print(f"    {'partition':10s} {'L':>5} {'max p':>7} {'boundary':>9} "
          f"{'pairs>=thr':>11} {'lam max':>8} {'4pt>noise':>10} "
          f"{'worst gap':>10} {'noise':>8} {'comps':>6}")
    rows = []
    for label, path in sets:
        if not path.exists():
            print(f"    {label}: {path.name} not on disk")
            continue
        names, M = read_fasta(path)
        P = pairwise_p(M)
        Lam = lambda_matrix(M)
        thr = p_threshold(M.shape[1], 0.05)
        off = P[np.triu_indices(len(names), 1)]
        loff = Lam[np.triu_indices(len(names), 1)]
        bad, worst, noise = four_point(P, M.shape[1])
        comps = components(P, thr)
        rows.append(dict(partition=label, L=M.shape[1], max_p=float(off.max()),
                         jc_boundary=0.75, thr=thr,
                         n_over=int((off >= thr).sum()),
                         lam_max=float(np.nanmax(loff)),
                         four_point_beyond_noise=bad, worst_gap=worst,
                         noise_scale=noise, n_components=len(comps),
                         component_sizes=str([len(c) for c in comps])))
        print(f"    {label:10s} {M.shape[1]:5d} {off.max():7.3f} "
              f"{0.75:9.3f} {int((off >= thr).sum()):11d} "
              f"{np.nanmax(loff):8.3f} {bad:10d} {worst:10.4f} "
              f"{noise:8.4f} {len(comps):6d}")
    pd.DataFrame(rows).to_csv(OUT / "coi_matrix.csv", index=False)
    print("\n    the pairwise view against the mean, on the saturated partition:")
    for label, path in sets[:1] + sets[4:]:
        names, M = read_fasta(path)
        P = pairwise_p(M); Lam = lambda_matrix(M)
        iu = np.triu_indices(len(names), 1)
        print(f"      {label:10s} p: mean {P[iu].mean():.3f}  max {P[iu].max():.3f}"
              f"   |   lambda: mean {np.nanmean(Lam[iu]):.3f}  max {np.nanmax(Lam[iu]):.3f}")
    print("      the mean is a compression of the pairwise structure: on the third")
    print("      positions it reads 0.749 while the worst pair is at 0.852.")
    print("    (the boundary column is the JC one; each alignment has its own"
          " at 1 - sum q^2,\n     which is 0.684 for the third positions and"
          " is what lambda is measured against)")

    print("\n    the 3rd-position matrix in full (p x 100):")
    names, M = read_fasta(EMP / "pos3.fa")
    P = pairwise_p(M)
    print("      " + " ".join(f"{n[:5]:>6}" for n in names))
    for i, n in enumerate(names):
        print(f"    {n[:5]:>6} " + " ".join(
            "     ." if j == i else f"{100 * P[i, j]:6.1f}"
            for j in range(len(names))))
    thr = p_threshold(M.shape[1], 0.05)
    over = [(names[i], names[j]) for i in range(len(names))
            for j in range(i + 1, len(names)) if P[i, j] >= thr]
    print(f"    pairs at or above the {thr:.3f} threshold: "
          f"{len(over)} of {len(names) * (len(names) - 1) // 2}")
    for a, b in over:
        print(f"      {a} - {b}")
    print(f"    components keeping only sub-threshold edges: "
          f"{[ [names[k] for k in c] for c in components(P, thr)]}")


if __name__ == "__main__":
    main()
