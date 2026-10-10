"""The critical tree length against rate heterogeneity, and the index at it.

Section 2.6 rests on a comparison that the manuscript reports as a row of
numbers: as alpha goes from infinity (no gamma) down to 0.5 the tree length at
which the true quartet stops being recoverable runs over a factor of 15, while
the index at that point moves by less than 0.01.  This script produces that
comparison as a curve, sampled densely at the no-gamma end so that the
position of the departure is visible.

Design: the N=4 quartet with the internal fraction used elsewhere in the
paper (0.35), L=1000, 150 replicates per point, and a geometric grid of tree
lengths wide enough to bracket the crossing for every alpha.  Recovery is
maximum likelihood by IQ-TREE, the criterion the critical value is defined
against, so nothing here is our own inference.

The sweep stops at alpha = 0.5 deliberately.  Below it a discrete gamma with
four categories stops representing a gamma: the slowest class is drawn from a
tail that a four-point quadrature resolves as a pile of effectively invariant
sites, not as a continuum.  Down there the recovery curve turns over and rises
again with tree length, which is a property of the discretisation rather than
of the model -- an earlier version of this experiment ran alpha down to 0.1
and spent a paragraph on that turnover before we recognised it.  A finer
categorisation would settle it; the point this section makes does not need it.

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
ALPHAS = [np.inf, 64.0, 32.0, 24.0, 16.0, 12.0, 8.0, 4.0, 2.0, 1.0, 0.5]
TL_GRID = [2, 3, 4.5, 7, 10, 15, 22, 33, 48, 70, 100, 150, 220, 320]
FI = 0.35
L = 1000
N_REPS = 150


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
        # round before matching: the grid is rebuilt from alpha on every run
        # and the CSV round-trip is not bit-exact (1100.861740811548 came back
        # as ...5479), which silently re-runs a point and appends it twice
        done = set(zip(prev.tag, prev.tl.round(6)))
        rows = prev.to_dict("records")
    t0 = time.time()
    for a in ALPHAS:
        tag = alpha_tag(a)
        m = model_for(a)
        for tl in TL_GRID:
            if (tag, round(float(tl), 6)) in done:
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

    analyze(rows)
    print(f"\ntotal {time.time() - t0:.0f}s")


def analyze(rows: list[dict] | None = None) -> pd.DataFrame:
    """Critical tree length per alpha, by interpolation in log TL.

    The crossing is taken as the FIRST downward one in tree length, not as
    "the last point above 0.95 against the first point below it".  The two
    coincide only if P_true is monotone, and read off a finite replicate
    count it never quite is: the old rule paired bracketing points from
    opposite sides of a rise and put alpha = 0.15 at TL_c = 4261 against a
    first crossing at 406.  The `rises` column reports how many upward steps
    each curve has, so the assumption is visible rather than implicit.
    """
    if rows is None:
        rows = pd.read_csv(OUT / "pilot.csv").to_dict("records")
    Q = build_Q("JC69")
    d = pd.DataFrame(rows).sort_values(["tag", "tl"])
    cres = []
    for tag, g in d.groupby("tag", sort=False):
        g = g.drop_duplicates(subset="tl").sort_values("tl")
        tl = g.tl.to_numpy(); p = g.ptrue.to_numpy()
        rec = dict(tag=tag, alpha=float(g.alpha.iloc[0]), n_pts=len(p),
                   rises=int((np.diff(p) > 0).sum()))
        i = next((k for k in range(len(p) - 1)
                  if p[k] >= 0.95 > p[k + 1]), None)
        if i is None:
            rec.update(tl_c=np.nan)
        else:
            x = np.log(tl[i:i + 2])              # increasing
            y = p[i:i + 2]                       # decreasing across the cross
            tl_c = float(np.exp(np.interp(0.95, y[::-1], x[::-1])))
            rec.update(tl_c=tl_c, tl_lo=float(tl[i]), tl_hi=float(tl[i + 1]),
                       p_lo=float(p[i]), p_hi=float(p[i + 1]))
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
    return c


if __name__ == "__main__":
    if "--analyze" in sys.argv:
        analyze()
    else:
        run()
