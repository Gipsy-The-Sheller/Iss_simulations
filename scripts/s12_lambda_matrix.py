"""What the pairwise lambda matrix sees that a single global index cannot.

The claim being tested is that lambda is not merely the mean p-distance in
other units.  It is not, once the mean is not what you look at: lambda is
defined per pair, and the author's own tooling plots the matrix as a heat map
precisely so that individual long branches can be located.  A single global
number -- Iss, mean p-distance, or mean lambda -- has no way to say *which*
sequences are the problem.

The design is the one from the review notes: six taxa, one of them on a long
pendant (1.6 expected substitutions against 0.05 for everything else), gamma
rates with alpha = 1, 4000 sites.  The global summaries are unremarkable; the
matrix is not.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench import trees                                            # noqa: E402
from issbench.models import build_Q                                   # noqa: E402
from issbench.sim import simulate                                     # noqa: E402
from issbench.stats import (h_bar, mean_p_distance, pi_hat, read_fasta)  # noqa: E402
from issbench.entropy import h, multinomial_entropy                   # noqa: E402

OUT = ROOT / "results" / "empirical"
LONG, SHORT = 1.6, 0.05
N_REPS, L = 3, 4000
LABELS = ["A", "B", "C", "D", "E", "F"]


def tree():
    """Unrooted six-taxon tree with F on a long pendant.

    Rooted representation with three children, so the root is a degree-3 node
    of the unrooted tree: (Z, E, F) with Z = ((A,B),(C,D)).
    """
    S = SHORT
    ab = trees.Node("", S, [trees.Node("A", S, []), trees.Node("B", S, [])])
    cd = trees.Node("", S, [trees.Node("C", S, []), trees.Node("D", S, [])])
    z = trees.Node("", S, [ab, cd])
    return trees.Node("", 0.0, [z, trees.Node("E", S, []),
                                trees.Node("F", LONG, [])])


def lam_matrix(M: np.ndarray) -> pd.DataFrame:
    """lambda_ij = delta_obs / delta_ran, both computed on the same sites."""
    n, L = M.shape
    comp = np.stack([(M == k).mean(axis=1) for k in range(4)])     # (4, n)
    out = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(i + 1, n):
            p = float((M[i] != M[j]).mean())
            d_ran = 1.0 - float((comp[:, i] * comp[:, j]).sum())
            out[i, j] = out[j, i] = p / d_ran
    return pd.DataFrame(out, index=LABELS, columns=LABELS)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t = tree()
    nwk = OUT / "lambda_tree.nwk"
    trees.write_newick(t, nwk)
    wd = ROOT / "work" / "empirical" / "lambda_matrix"
    files = simulate(nwk, "JC69+G4{1.0}", L, N_REPS, 7, wd, tag="lam")
    print(f"simulated {len(files)} replicates, {L} sites, 6 taxa, "
          f"F pendant {LONG} vs {SHORT} elsewhere")

    rows, mats = [], []
    for k, f in enumerate(files):
        _names, M = read_fasta(f)
        q = pi_hat(M)
        iss = h_bar(M) / multinomial_entropy(M.shape[0], q)
        lm = lam_matrix(M)
        off = lm.to_numpy()[np.triu_indices(6, 1)]
        f_pairs = [v for i in range(6) for j in range(i + 1, 6)
                   if "F" in (LABELS[i], LABELS[j]) for v in [lm.iloc[i, j]]]
        rest = [v for i in range(6) for j in range(i + 1, 6)
                if "F" not in (LABELS[i], LABELS[j]) for v in [lm.iloc[i, j]]]
        rows.append(dict(rep=k, iss_A=iss, pbar=mean_p_distance(M),
                         lam_mean=float(np.mean(off)),
                         lam_F=float(np.mean(f_pairs)),
                         lam_rest=float(np.mean(rest)),
                         bound=1 - float((q ** 2).sum())))
        mats.append(lm)
        print(f"  rep{k}: Iss={iss:.3f}  mean lam={np.mean(off):.3f}  "
              f"F-pairs={np.mean(f_pairs):.3f}  other={np.mean(rest):.3f}")

    d = pd.DataFrame(rows)
    d.to_csv(OUT / "lambda_matrix_summary.csv", index=False)
    mean_mat = sum(mats) / len(mats)
    mean_mat.to_csv(OUT / "lambda_matrix.csv")
    print("\nmean lambda over replicates:")
    print(mean_mat.round(3).to_string())
    print(f"\nbound 1-Sum q^2 = {d.bound.mean():.3f}   (0.75 only under "
          f"uniform composition)")


if __name__ == "__main__":
    main()
