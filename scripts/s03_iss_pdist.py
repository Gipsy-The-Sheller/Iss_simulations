"""Item 2 -- Iss vs mean p-distance under JC69, with the closed forms overlaid.

Design: the uniform-Iss sampler of Item 0, so the points are spread evenly in
Iss rather than evenly in tree length (the latter wastes the budget in the
saturated plateau).  Both conventions are recorded so the c-ratio can be
checked against Eq. (3) of the manuscript.

N=4 is the main arm; N=8 and N=16 spot points are included to show the relation
survives the number of taxa.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench import design as D                                       # noqa: E402
from issbench.models import build_Q                                    # noqa: E402
from issbench.runner import dump_meta, run_cells                       # noqa: E402
from issbench.trees import star, write_newick                          # noqa: E402

EXP = "exp2_iss_pdist"
TREES = ROOT / "results" / EXP / "trees"
N_REPS = 200
MODEL = "JC69"
N_SPOT = (4, 8, 16)   # N=4 is a star, so the N-comparison is shape-controlled
N_SPOT_POINTS = 8


def quartet_grid():
    """60 s-values uniform in Iss on [10%, 90%] of the attainable range."""
    d = D.quartet_design(build_Q("JC69"), np.ones(4) / 4, internal_fraction=0.35,
                         N=4, alpha=np.inf, s_min=0.05, s_max=60.0)
    ss, mus = d.mu_grid(300)
    lo, hi = mus.min(), mus.max()
    tgt = lo + np.linspace(0.10, 0.90, 60) * (hi - lo)
    o = np.argsort(mus)
    return np.interp(tgt, mus[o], ss[o]), d


def star_grid(N: int):
    d = D.star_design(build_Q("JC69"), np.ones(4) / 4, N=N, alpha=np.inf,
                      s_min=0.1, s_max=200.0)
    ss, mus = d.mu_grid(300)
    lo, hi = mus.min(), mus.max()
    tgt = lo + np.linspace(0.15, 0.90, N_SPOT_POINTS) * (hi - lo)
    o = np.argsort(mus)
    return np.interp(tgt, mus[o], ss[o])


def cells():
    TREES.mkdir(parents=True, exist_ok=True)
    # ---- N = 4 quartet, uniform-Iss design
    grid, _ = quartet_grid()
    for k, s in enumerate(grid):
        r = 0.35 * s
        l = (s - r) / 4.0
        import issbench.trees as T
        node = T.quartet(l, l, l, l, r)
        f = TREES / f"q4_s{k:03d}.nwk"
        write_newick(node, f)
        for L in (1000, 5000):
            yield dict(tag=f"q4_s{k:03d}_L{L}", tree_file=f,
                       tree_nwk=T.to_newick(node), model=MODEL, L=L,
                       n_reps=N_REPS, design="uniform_iss", s=float(s))
    # ---- N = 8, 16 spot points (star)
    for N in N_SPOT:
        for k, s in enumerate(star_grid(N)):
            node = star(N, s / N)
            f = TREES / f"star{N}_s{k:02d}.nwk"
            write_newick(node, f)
            yield dict(tag=f"star{N}_s{k:02d}_L2000", tree_file=f,
                       tree_nwk=f"star{N}", model=MODEL, L=2000,
                       n_reps=N_REPS, design="spot", s=float(s))


if __name__ == "__main__":
    dump_meta(EXP, dict(item=2, model=MODEL, alpha="inf",
                        n_reps=N_REPS, N_spot=list(N_SPOT)))
    n = run_cells(EXP, cells(), ROOT)
    print(f"\n{EXP}: {n} new replicate rows")
