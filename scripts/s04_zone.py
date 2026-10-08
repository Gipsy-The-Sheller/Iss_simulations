"""Item 3 -- Felsenstein -> Farris transfer.

    (a, b, c, d) = (L - t, S, L, S + t),   t in [0, L - S],  internal r fixed

t = 0      : (L,S,L,S) -- long branches A and C separated by r  (Felsenstein)
t = L - S  : (S,S,L,L) -- long branches C and D adjacent        (Farris)

The true topology ((A,B),(C,D)) never changes; the sum of pairwise distances is
constant along the path (3(a+b+c+d)+4r); and at the endpoints the distance
multisets coincide exactly as r -> 0, so Iss is identical there.

Recovery is measured three ways -- ML (IQ-TREE), MP (Fitch) and NJ (both
JC-corrected and the uncorrected four-point rule) -- because the point is
precisely that recovery depends on the branch-length arrangement while Iss does
not.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench import trees                                             # noqa: E402
from issbench.infer import ptrue                                       # noqa: E402
from issbench.models import build_Q                                    # noqa: E402
from issbench.runner import (dump_meta, out_path, run_cells,          # noqa: E402
                            work_dir)
from issbench.stats import read_fasta                                  # noqa: E402
from issbench.trees import write_newick                                # noqa: E402

EXP = "exp3_zone"
TREES = ROOT / "results" / EXP / "trees"
P = 1.0          # long pendant
Q = 0.05         # short pendant
R = 0.02         # internal branch
N_T = 17         # points along the path
N_REPS = 200
L = 1000
MODEL = "JC69"


def cells():
    TREES.mkdir(parents=True, exist_ok=True)
    for k, t in enumerate(np.linspace(0.0, P - Q, N_T)):
        node = trees.quartet(P - t, Q, P, Q + t, R)
        f = TREES / f"t{k:02d}.nwk"
        write_newick(node, f)
        yield dict(tag=f"t{k:02d}", tree_file=f, tree_nwk=trees.to_newick(node),
                   model=MODEL, L=L, n_reps=N_REPS, t=float(t),
                   a=float(P - t), b=float(Q), c=float(P), d=float(Q + t))


def recovery():
    """P_true per path point for each method; written next to cells.csv.

    ML shells out to IQ-TREE once per replicate, so this is the expensive half:
    at N_T=17 and n_reps=200 that is 3400 IQ-TREE runs of 4 taxa, minutes each
    in aggregate on one core.  MP / NJ / 4P are computed in-process.
    """
    rows = []
    for k, t in enumerate(np.linspace(0.0, P - Q, N_T)):
        wd = work_dir(EXP, f"t{k:02d}", ROOT)
        files = sorted(wd.glob(f"t{k:02d}*.fa"))
        if not files:
            print(f"  [t{k:02d}] no alignments found, skipping", flush=True)
            continue
        rec = dict(t=t)
        for meth in ("ML", "MP", "NJ", "4P"):
            rec[f"Ptrue_{meth}"] = ptrue(files, meth, MODEL, workdir=wd)
        rows.append(rec)
        print(f"  t={t:.3f}  " + "  ".join(
            f"{m}={rec[f'Ptrue_{m}']:.3f}" for m in ("ML", "MP", "NJ", "4P")),
            flush=True)
    if rows:
        pd.DataFrame(rows).to_csv(out_path(EXP).parent / "recovery.csv",
                                  index=False)
        print(f"wrote {out_path(EXP).parent / 'recovery.csv'}")


if __name__ == "__main__":
    build_Q("JC69")
    dump_meta(EXP, dict(item=3, P=P, Q=Q, R=R, n_t=N_T, n_reps=N_REPS,
                        L=L, model=MODEL, methods=["ML", "MP", "NJ", "4P"]))
    n = run_cells(EXP, cells(), ROOT)
    print(f"\n{EXP}: {n} new replicate rows; now measuring recovery")
    recovery()
