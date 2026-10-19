"""Is convergence slowest in the Farris zone?

Xia et al. report that the asymmetrical ultrametric tree's Iss.c converges
more slowly than the symmetrical one's.  If p = lam_Iss / lam_Delta is right,
the rate is set by how fast the topological signal decays with tree length,
and that is fastest where the long branches are sisters -- so the prediction
on the table is that p is smallest in the Farris zone.

The prediction is plausible but not automatic, and the distinction it turns
on has to be kept straight:

    difficulty   is Delta sqrt(L) / sqrt(v2) at a given tree;
    p            is the rate at which the critical tree moves with L.

A tree can be hard to resolve and still converge quickly, and vice versa.  So
this is a measurement, not a corollary.

One controlled family.  Hold the pendant sum fixed at 4 -- so the mean
pendant is 1 throughout -- and move only the arrangement:

    Farris(r)       (1+r, 1+r, 1-r, 1-r)   long branches are sisters
    Felsenstein(r)  (1+r, 1-r, 1+r, 1-r)   long branches are apart

with r = 0 giving the same all-equal tree on both sides, so any difference at
r > 0 is attributable to the arrangement alone and not to the scale.  The
internal branch is held at 0.2 throughout.

Reported for each shape: p, and the critical tree length at a common L, so
that difficulty and rate can be compared rather than conflated.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PART2 = HERE.parent
ROOT = PART2.parent
sys.path.insert(0, str(ROOT / "src"))

OUT = PART2 / "results" / "farris_slowest"

_spec = importlib.util.spec_from_file_location(
    "s36", str(HERE / "s36_factors.py"))
s36 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s36)

ISS_SAT = s36.ISS_SAT
LS = (30000, 100000, 300000, 1000000)
E_INT = 0.2


def family(r):
    return {"Farris": (1 + r, 1 + r, 1 - r, 1 - r, E_INT),
            "Felsenstein": (1 + r, 1 - r, 1 + r, 1 - r, E_INT)}


def measure(shape):
    """(p, r^2, critical tree length at LS[0], its Iss)."""
    ss, gs, used = [], [], []
    for L in LS:
        s = s36.boundary_fast(shape, L, target=0.95)
        if not np.isfinite(s):
            continue
        ss.append(s)
        gs.append(ISS_SAT - s36.e_iss(shape, s))
        used.append(L)
    if len(used) < 3:
        return float("nan"), float("nan"), float("nan"), float("nan")
    L_ = np.array(used, float)
    g_ = np.array(gs)
    A = np.column_stack([np.ones_like(L_), -np.log(L_)])
    c, *_ = np.linalg.lstsq(A, np.log(g_), rcond=None)
    pred = np.exp(c[0]) * L_ ** (-c[1])
    r2 = 1 - float(((g_ - pred) ** 2).sum()) / float(((g_ - g_.mean()) ** 2).sum())
    return float(c[1]), r2, float(ss[0]), float(s36.e_iss(shape, ss[0]))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print("Convergence rate against pendant arrangement, pendant sum fixed")
    print("at 4 and internal branch at 0.2.  r = 0 is the same all-equal")
    print("tree for both families, so any difference is the arrangement.\n")
    print(f"  {'r':>5} {'family':>12} {'p':>8} {'r2':>7} "
          f"{'TL* at 3e4':>11} {'Iss*':>8} {'difficulty z':>13}")
    rows = []
    for r in (0.0, 0.2, 0.4, 0.6, 0.8):
        for fam, shape in family(r).items():
            p, r2, s0, iss0 = measure(shape)
            if not np.isfinite(p):
                print(f"  {r:5.1f} {fam:>12}   no crossing")
                continue
            # difficulty at the SAME L, measured as the critical scale
            # relative to the r = 0 case at the same L: bigger TL* means
            # the topology stays resolvable further into saturation, i.e.
            # easier.
            rows.append(dict(r=r, family=fam, p=p, r2=r2, s_crit=s0,
                             TL_crit=sum(shape) * s0, iss_crit=iss0))
            print(f"  {r:5.1f} {fam:>12} {p:8.4f} {r2:7.4f} "
                  f"{sum(shape) * s0:11.3f} {iss0:8.4f} "
                  f"{'-' if not rows else ''}", flush=True)
        print()
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "farris_slowest.csv", index=False)

    if len(d):
        print("  p against r, by family:")
        for fam in d.family.unique():
            sub = d[d.family == fam].sort_values("r")
            print(f"    {fam:12s} r: {list(sub.r)}  p: "
                  f"{[round(v, 4) for v in sub.p]}")
        print("\n  Xia's observation is that the asymmetrical tree converges")
        print("  more slowly.  If p falls with r in BOTH families the effect")
        print("  is about asymmetry, not about which branches are sisters;")
        print("  if it falls only for Farris, it is about the Farris zone.")
        print("  Difficulty is the TL_crit column: larger means resolvable")
        print("  further into saturation, i.e. easier at the same L.")


if __name__ == "__main__":
    main()
