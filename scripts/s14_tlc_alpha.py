"""The critical tree length against rate heterogeneity, and the index at it.

Section 2.6 rests on a comparison that the manuscript currently reports as
four numbers in a sentence: as alpha goes from infinity to 0.5 the tree length
at which the true quartet stops being recoverable moves by a factor of 18,
while the index at that point moves by less than 0.01.  This script produces
that comparison as a curve.

Design: the N=4 quartet with the internal fraction used elsewhere in the
paper (0.35), L=1000, 150 replicates per point, and a geometric grid of tree
lengths wide enough to bracket the crossing for every alpha.  Recovery is
maximum likelihood by IQ-TREE, the criterion the critical value is defined
against, so nothing here is our own inference.

Two outputs per alpha: the tree length at which P_true crosses 0.95, and the
exact E[Iss] at that tree length (Theorem exact, evaluated by enumeration, not
read off the simulation).
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench import exact, trees                                    # noqa: E402
from issbench.gamma import discrete_gamma_rates                       # noqa: E402
from issbench.infer import TRUE_SPLIT, iqtree_ml, split_of_newick     # noqa: E402
from issbench.models import build_Q                                   # noqa: E402
from issbench.sim import simulate                                     # noqa: E402

EXP = "exp7_tlc_alpha"
OUT = ROOT / "results" / EXP
ALPHAS = [np.inf, 64.0, 32.0, 24.0, 16.0, 12.0, 8.0, 4.0, 2.0, 1.0, 0.5,
          0.4, 0.3, 0.25, 0.2, 0.15, 0.1]
TL_GRID = [2, 3, 4.5, 7, 10, 15, 22, 33, 48, 70, 100, 150, 220, 320]
FI = 0.35
L = 1000
N_REPS = 150


def grid_for(a: float) -> list[float]:
    """Tree lengths to sample for this alpha.

    For alpha >= 0.5 the fixed grid is used.  Below that the critical tree
    length grows steeply -- an order of magnitude per halving, roughly
    TL_c ~ alpha^-2.08 fitted to the upper row of the grid -- so a single grid
    cannot bracket every crossing without wasting most of its points in the
    plateau.  Each smaller alpha gets a geometric grid spanning the predicted
    crossing by a factor of four below and six above.
    """
    if np.isinf(a) or a >= 0.5:
        return list(TL_GRID)
    pred = 25.47 * a ** -2.08            # fitted to alpha = 0.5, 1 in the grid
    lo, hi = pred / 4.0, pred * 6.0
    return list(np.geomspace(lo, hi, 7))


def model_for(a: float) -> str:
    return "JC69" if np.isinf(a) else f"JC69+G4{{{a:g}}}"


def alpha_tag(a: float) -> str:
    return "inf" if np.isinf(a) else f"a{a:g}"


def tree_for(tl: float):
    r = FI * tl
    l = (tl - r) / 4.0
    return trees.quartet(l, l, l, l, r)


def run() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "trees").mkdir(exist_ok=True)
    Q = build_Q("JC69")
    rows = []
    done = set()
    f = OUT / "pilot.csv"
    if f.exists():
        prev = pd.read_csv(f)
        done = set(zip(prev.tag, prev.tl))
        rows = prev.to_dict("records")
    t0 = time.time()
    for a in ALPHAS:
        tag = alpha_tag(a)
        m = model_for(a)
        for tl in grid_for(a):
            if (tag, float(tl)) in done:
                continue
            node = tree_for(float(tl))
            nwk = OUT / "trees" / f"{tag}_{tl:.4g}.nwk"
            trees.write_newick(node, nwk)
            wd = ROOT / "work" / EXP / f"{tag}_{tl:.4g}"
            try:
                files = simulate(nwk, m, L, N_REPS, 11, wd,
                                 tag=f"{tag}_{tl:.4g}")
            except Exception as e:                                    # noqa: BLE001
                print(f"  [{tag} TL={tl}] SIM FAILED: {e}", flush=True)
                continue
            hits = 0
            for k, fa in enumerate(files):
                tr = iqtree_ml(fa, "JC69", seed=1 + k, workdir=fa.parent)
                hits += int(split_of_newick(tr) == TRUE_SPLIT)
            p = hits / len(files)
            rows.append(dict(alpha=float(a), tag=tag, tl=float(tl), ptrue=p,
                             n=len(files)))
            pd.DataFrame(rows).to_csv(f, index=False)
            print(f"  {tag:5s} TL={tl:>6g}  P_true={p:.3f}  "
                  f"({time.time() - t0:.0f}s)", flush=True)

    # critical tree length by linear interpolation on the log tree length,
    # which is how P_true falls; report the bracketing points too
    d = pd.DataFrame(rows).sort_values(["tag", "tl"])
    cres = []
    for tag, g in d.groupby("tag", sort=False):
        g = g.sort_values("tl")
        above = g[g.ptrue >= 0.95]
        below = g[g.ptrue < 0.95]
        rec = dict(tag=tag, alpha=float(g.alpha.iloc[0]))
        if len(above) and len(below):
            lo, hi = above.iloc[-1], below.iloc[0]
            x = np.log([lo.tl, hi.tl]); y = [lo.ptrue, hi.ptrue]
            tl_c = float(np.exp(np.interp(0.95, y[::-1], x[::-1])))
            rec.update(tl_c=tl_c, tl_lo=float(lo.tl), tl_hi=float(hi.tl),
                       p_lo=float(lo.ptrue), p_hi=float(hi.ptrue))
        else:
            rec.update(tl_c=np.nan)
        # exact E[Iss] at that tree length, by enumeration at N=4
        if not np.isnan(rec.get("tl_c", np.nan)):
            node = tree_for(rec["tl_c"])
            rates = (np.array([1.0]) if np.isinf(rec["alpha"])
                     else discrete_gamma_rates(rec["alpha"], 4))
            w = np.full(len(rates), 1.0 / len(rates))
            rec["iss_c"] = exact.E_Iss(4, node, Q, np.ones(4) / 4,
                                       tuple(rates), tuple(w), convention="A")
        cres.append(rec)
    c = pd.DataFrame(cres)
    c.to_csv(OUT / "critical.csv", index=False)
    print()
    print(c.to_string(index=False))
    print(f"\ntotal {time.time() - t0:.0f}s")


if __name__ == "__main__":
    run()
