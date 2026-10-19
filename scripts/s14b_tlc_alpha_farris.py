"""Figure 2's experiment, run on a Farris arrangement instead of a symmetric one.

s14 produced Part I's Figure 2 and concluded that rate heterogeneity moves the
critical tree length by a factor of fifteen while moving the index at that
tree length by less than 0.01, attributing the near-invariance to the index
being close to its ceiling.  That attribution is right about the configuration
and wrong as a statement about the index.

s14's tree is perfectly symmetric: four equal pendants and an internal branch
carrying 0.35 of the tree length.  Such a tree is easy at every tree length,
its critical index sits at 0.97--0.98, and a quantity pinned just below its
ceiling has nowhere to move.  Change only the arrangement --- put the long
pendants on sisters, the Farris configuration --- and the same sweep gives an
index that runs from 0.77 to 0.88.

Nothing else is changed here.  Same model set, same alignment length, same
replicate count, same alpha grid, same tree-length grid, same recovery method,
same criterion for locating the crossing.  The only difference is
`tree_for`, and the comparison it licenses is therefore about the arrangement
and nothing else.

Outputs
    results/exp7b_tlc_alpha_farris/pilot.csv
    results/exp7b_tlc_alpha_farris/critical.csv
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

EXP = "exp7b_tlc_alpha_farris"
OUT = ROOT / "results" / EXP
ALPHAS = [np.inf, 64.0, 32.0, 24.0, 16.0, 12.0, 8.0, 4.0, 2.0, 1.0, 0.5]
TL_GRID = [1.5, 2, 3, 4.5, 7, 10, 15, 22, 33, 48, 70, 100, 150, 220, 320]
FI = 0.35          # internal branch as a fraction of tree length, as in s14
RATIO = 20.0       # long pendant / short pendant; 1.0 would be s14's tree
L = 1000
N_REPS = 150


def model_for(a: float) -> str:
    return "JC69" if np.isinf(a) else f"JC69+G4{{{a:g}}}"


def alpha_tag(a: float) -> str:
    return "inf" if np.isinf(a) else f"a{a:g}"


def tree_for(tl: float):
    """The same tree as s14 except that the long pendants are sisters.

    s14 used four equal pendants.  Here the pendant sum on each side of the
    internal branch is the same, so the tree length and the internal fraction
    are unchanged, and only the split of that sum between the two pendants
    differs.
    """
    r = FI * tl
    half = (tl - r) / 2.0
    p = half * RATIO / (1.0 + RATIO)          # long, on A and B
    q = half / (1.0 + RATIO)                  # short, on C and D
    return trees.quartet(p, p, q, q, r)


def run() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "trees").mkdir(exist_ok=True)
    Q = build_Q("JC69")
    rows = []
    done = set()
    f = OUT / "pilot.csv"
    if f.exists():
        prev = pd.read_csv(f)
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
    """Critical tree length per alpha, exactly as in s14."""
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
            x = np.log(tl[i:i + 2])
            y = p[i:i + 2]
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
