"""Item 1 -- (tree length x sequence length) maps of Iss location and scale.

No presupposed conclusion: the run produces the full 2-D grid so that the four
candidate invariances (median vs L, mean vs L, SD vs L, SD vs divergence) can
each be read off where they hold.

Two quartets per tree length: equal pendants ("symmetric") and a graded ladder
("asymmetric"), so the topology axis is a branch-length arrangement, which is
the only sense in which a 4-taxon unrooted tree can be more or less balanced.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench import trees                                             # noqa: E402
from issbench.models import build_Q                                    # noqa: E402
from issbench.runner import dump_meta, run_cells                       # noqa: E402
from issbench.trees import Node, write_newick                          # noqa: E402

EXP = "exp1_Lmap"
# TL = 2 already gives mean pairwise ~1.12 (pbar ~ 0.58, Iss ~ 0.71), so a
# grid starting at 2 misses the whole low-divergence regime; and TL = 200 is
# deep in the plateau.  Span 0.2 to 100 instead.
TL_GRID = np.exp(np.linspace(np.log(0.2), np.log(100.0), 28))
L_GRID = [100, 200, 500, 1000, 2000, 5000, 10000]
N_REPS = 200            # per (TL, L, topology); raise for the final figure
ALPHA = np.inf          # homogeneity first; add +G4 arms later
MODEL = "JC69"
TREES = ROOT / "results" / EXP / "trees"


# Both shapes use the SAME internal-branch fraction (10% of TL); they differ
# only in the pendant arrangement, so the topology axis is clean.  (An earlier
# version gave the symmetric arm a zero internal branch, i.e. a star, which
# confounded pendant equality with loss of the internal branch.)
FI = 0.10


def sym(TL: float) -> Node:
    """Equal pendants."""
    r = FI * TL
    l = (TL - r) / 4.0
    return trees.quartet(l, l, l, l, r)


def asym(TL: float) -> Node:
    """Graded ladder: pendants in ratio 1 : 2 : 3 : 4."""
    r = FI * TL
    rest = TL - r
    w = np.array([1, 2, 3, 4], float)
    p = rest * w / w.sum()
    return trees.quartet(p[0], p[1], p[2], p[3], r)


def cells():
    TREES.mkdir(parents=True, exist_ok=True)
    for shape, fn in (("sym", sym), ("asym", asym)):
        for TL in TL_GRID:
            f = TREES / f"{shape}_TL{TL:.3f}.nwk"
            write_newick(fn(float(TL)), f)
            for L in L_GRID:
                yield dict(tag=f"{shape}_TL{TL:.2f}_L{L}", tree_file=f,
                           tree_nwk=trees.to_newick(fn(float(TL))),
                           model=MODEL, L=L, n_reps=N_REPS,
                           shape=shape, TL=float(TL))


if __name__ == "__main__":
    build_Q("JC69")
    dump_meta(EXP, dict(item=1, model=MODEL, alpha="inf",
                        tl_grid=[float(x) for x in TL_GRID], L_grid=L_GRID,
                        n_reps=N_REPS, topology=["sym", "asym"]))
    n = run_cells(EXP, cells(), ROOT)
    print(f"\n{EXP}: {n} new replicate rows")
