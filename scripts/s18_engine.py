"""The critical value is a property of the engine, not of the alignment.

The 2003 calibration is the joint product of a data generator and an inference
engine, and the source says so: reconstruction used neighbour-joining with the
F84 distance and maximum likelihood under the same F84 model for NOTU = 4, 8
and 16, and "only the neighbour-joining method is used for reconstructing
topologies with NOTU > 16" (Xia et al. 2003, Material and Methods).  The
tabulated Iss.c is therefore an Iss.c for that pair, and the engine is not a
coordinate of the table.

This script measures how large that omission is, using the alignments already
simulated for the alpha sweep: the same 150 replicates per tree length, the
same generator, the same N and L and the same definition of the critical tree
length, scored twice -- by neighbour joining in process and by maximum
likelihood in IQ-TREE.  The exact E[Iss] at each crossing comes from the same
enumeration as Theorem 8, so the comparison is between engines only.

    python scripts/s18_engine.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench import exact, trees                                    # noqa: E402
from issbench.infer import TRUE_SPLIT, nj_split, split_of_newick, iqtree_ml  # noqa: E402
from issbench.models import build_Q                                   # noqa: E402
from issbench.stats import read_fasta                                 # noqa: E402

WORK = ROOT / "work" / "exp7_tlc_alpha"
OUT = ROOT / "results" / "engine"
TAG = "inf"                      # the no-gamma arm, so the comparison is clean
FI = 0.35
L = 1000
N_REPS = 150


def tree_for(tl: float):
    r = FI * tl
    return trees.quartet((tl - r) / 4, (tl - r) / 4, (tl - r) / 4,
                         (tl - r) / 4, r)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    dirs = sorted(d for d in WORK.glob(f"{TAG}_*") if d.is_dir())
    if not dirs:
        sys.exit(f"no alignments under {WORK} -- run scripts/s14_tlc_alpha.py")
    rows = []
    for d in dirs:
        tl = float(d.name.split("_", 1)[1])
        fas = sorted(d.glob(f"{d.name}_r*.fa"))
        if len(fas) < N_REPS:
            print(f"  {d.name}: only {len(fas)} replicates, skipped")
            continue
        fas = fas[:N_REPS]
        nj = ml = 0
        for f in fas:
            _, M = read_fasta(f)
            nj += int(nj_split(M) == TRUE_SPLIT)
            ml += int(split_of_newick(iqtree_ml(f, "JC69", seed=1,
                                                workdir=f.parent)) == TRUE_SPLIT)
        rows.append(dict(tl=tl, n=len(fas), p_nj=nj / len(fas), p_ml=ml / len(fas)))
        print(f"  TL={tl:>7g}  NJ {rows[-1]['p_nj']:.3f}  ML "
              f"{rows[-1]['p_ml']:.3f}", flush=True)

    d = pd.DataFrame(rows).sort_values("tl")
    d.to_csv(OUT / "nj_vs_ml.csv", index=False)

    Q = build_Q("JC69")
    pi = np.ones(4) / 4

    def cross(col: str) -> tuple[float, float]:
        """Tree length at which P_true crosses 0.95, and the exact E[Iss] there."""
        tl, p = d.tl.to_numpy(), d[col].to_numpy()
        i = next((k for k in range(len(p) - 1)
                  if p[k] >= 0.95 > p[k + 1]), None)
        if i is None:
            return np.nan, np.nan
        x = np.log(tl[i:i + 2])
        tl_c = float(np.exp(np.interp(0.95, p[i:i + 2][::-1], x[::-1])))
        return tl_c, exact.E_Iss(4, tree_for(tl_c), Q, pi, convention="A")

    print(f"\n{'engine':28s} {'TL_c':>9} {'Iss.c':>8}")
    res = []
    for label, col in (("neighbour joining, JC distance", "p_nj"),
                       ("maximum likelihood, IQ-TREE", "p_ml")):
        tl_c, iss_c = cross(col)
        res.append(dict(engine=label, tl_c=tl_c, iss_c=iss_c))
        print(f"  {label:28s} {tl_c:9.3f} {iss_c:8.4f}")
    pd.DataFrame(res).to_csv(OUT / "critical_by_engine.csv", index=False)

    nj_c, ml_c = res[0]["iss_c"], res[1]["iss_c"]
    print(f"\n  gap {abs(ml_c - nj_c):.4f}, against the {0.0115:.4f} the alpha "
          f"sweep spans\n  over two orders of magnitude of rate heterogeneity "
          f"(ratio {abs(ml_c - nj_c) / 0.0115:.0f}x)")
    print("  same alignments, same generator, same N, L and definition of the\n"
          "  crossing; only the inference engine differs.")


if __name__ == "__main__":
    main()
