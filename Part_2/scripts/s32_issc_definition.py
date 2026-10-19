"""Stage 2, on our own definition.

Xia's tabulated Iss.c is not used.  It was produced with DNAML under F84,
with NJ mixed in, and its convention does not match any Iss we can compute
(the ceiling of our scale is 1.0 in convention A and 0.662 in convention B,
against roughly 0.81 in their table) -- so a level comparison is a
comparison of two different quantities wearing the same name.

Iss.c is therefore DEFINED here, the way it has to be defined to mean
anything: the Iss of the tree at which resampling recovers the topology
just 95% of the time.  That is a definition, not a fit, and it needs no
external data.  The question then becomes one of internal consistency, which
is a stronger position:

    Iss.c^MC(L)   -- found by bisecting the resampling recovery rate to 0.95
    Iss.c^pred(L) -- found by solving Delta(TL) sqrt(L) = z sqrt(v2(TL))

against each other.  If the analytic boundary is right, the two agree -- and
they are compared in Iss, not in probability, which is a harder test because
it goes through the map TL -> Iss.  Near saturation that map is very flat,
so a small error in the critical tree barely shows in Iss; away from it, the
map has real slope and an error in TL* shows up amplified.

Two recovery rules are used, and the difference between them is the Stage 1
correction rather than an unknown:

  * oracle -- topology 1 at its true parameters, the alternatives at their
    KL optima.  Exact sampling by multinomial, cheap, and it is what the
    analytic formula actually describes.
  * profiled -- every topology maximised over its five branch lengths.  What
    a program really does, and slower.

The oracle is what the analytic boundary can be expected to match.  Comparing
the two also measures how much of Iss.c is definition and how much is
estimator.
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

OUT = PART2 / "results" / "issc_definition"

_spec = importlib.util.spec_from_file_location(
    "s27", str(HERE / "s27_mc_boundary.py"))
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)

from issbench import exact, trees                                  # noqa: E402
from issbench.models import build_Q                                 # noqa: E402

Q = build_Q("JC69")
PI = np.ones(4) / 4
ISS_SAT = float(exact.E_Iss(4, trees.quartet(40.0, 40.0, 40.0, 40.0, 40.0),
                            Q, PI, convention="A"))

# The Farris zone, which is where this belongs.  There is no reason to
# reproduce Xia's ultrametric families: they are not Farris configurations
# at all (the balanced one has all four pendants equal, the asymmetrical one
# is a caterpillar), and the whole point here is the region where the long
# pendants are sisters and the tree is hardest to resolve.
#
# Leaves A and B are sisters, so long branches on A and B together is the
# aligned case; c = d keeps the (A B) symmetry that makes the two wrong
# topologies coincide exactly.
SHAPES = {
    "a=b=1.0, c=d=0.05, e=0.02": (1.0, 1.0, 0.05, 0.05, 0.02),
    "a=b=1.0, c=d=0.05, e=0.10": (1.0, 1.0, 0.05, 0.05, 0.10),
    "a=b=1.0, c=d=0.20, e=0.05": (1.0, 1.0, 0.20, 0.20, 0.05),
    "a=b=0.6, c=d=0.05, e=0.02": (0.6, 0.6, 0.05, 0.05, 0.02),
}
# These shapes are not resolvable at small L at ANY tree length -- s27 found
# the best achievable rate at L = 10^4 is still only 0.90 -- so the sweep
# has to start where recovery exists at all.
LS = (30000, 60000, 100000, 200000)
TARGET = 0.95


def e_iss(b1):
    return float(exact.E_Iss(4, trees.quartet(*b1), Q, PI, convention="A"))


def rate_at(shape, s, L, n_rep, kind, seed):
    """Recovery rate at tree scale s, by exact multinomial resampling."""
    b1 = tuple(v * s for v in shape)
    p_true, fits, es = m.three_fits(b1)
    rng = np.random.default_rng(seed)
    if kind == "oracle":
        return m.mc(p_true, fits, L, n_rep, rng), b1
    counts = rng.multinomial(L, p_true, size=n_rep)
    starts = [np.asarray(b1, float)]
    for lab in m.SPLITS:
        if lab == "12|34":
            continue
        _, bb = m.best_fit(p_true, m.SPLITS[lab], b1)
        starts.append(np.asarray(bb, float))
    labs = list(m.SPLITS)
    wins = 0
    for i in range(n_rep):
        vals = [s28.profile(counts[i], m.SPLITS[lab], starts[a])[0]
                for a, lab in enumerate(labs)]
        top = max(vals)
        tied = [j for j, v in enumerate(vals) if top - v <= 1e-6]
        if 0 in tied:
            wins += 1.0 / len(tied)
    return wins / n_rep, b1


def critical_scale(shape, L, n_rep, kind, grid=np.geomspace(0.05, 60.0, 13),
                   seed=17):
    """Largest scale at which the resampling recovery rate still reaches 0.95.

    Scanned for the maximum before bisecting, because recovery is not
    monotone in tree length: a tree shrunk to a point is a star where the
    answer is a coin toss, and a tree grown past saturation is noise where
    it is again, so the rate peaks in between and a bracket on the endpoints
    can straddle the peak and report a crossing that is not there.
    """
    cache = {}

    def g(s):
        key = round(float(s), 10)
        if key not in cache:
            cache[key], _ = rate_at(shape, float(s), L, n_rep, kind, seed)
        return cache[key]

    vals = np.array([g(s) for s in grid])
    j = int(np.argmax(vals))
    if vals[j] < TARGET:
        return float("nan")
    for i in range(j, len(grid) - 1):
        if not (vals[i] >= TARGET > vals[i + 1]):
            continue
        lo, hi = grid[i], grid[i + 1]
        for _ in range(12):
            mid = 0.5 * (lo + hi)
            if g(mid) >= TARGET:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1e-5 * hi:
                break
        return 0.5 * (lo + hi)
    return float("nan")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    global s28
    _s28 = importlib.util.spec_from_file_location(
        "s28", str(HERE / "s28_stage1_profile.py"))
    s28 = importlib.util.module_from_spec(_s28)
    _s28.loader.exec_module(s28)

    print("STAGE 2 (own definition).  Iss.c := the Iss of the tree at which")
    print("our own resampling recovers the topology just 95% of the time.")
    print(f"No external data.  Iss scale: convention A, ceiling {ISS_SAT:.4f}.\n")

    rows = []
    print(f"  {'shape':40s} {'L':>7} {'s*_MC':>8} {'TL*':>7} "
          f"{'Iss.c MC':>9} {'+/-':>7} {'Iss.c pred':>11} {'diff':>9} "
          f"{'dP/ds':>7} {'n':>6}")
    for name, shape in SHAPES.items():
        for L in LS:
            # The transition is FLAT in this zone: dP/ds is about -0.27
            # per unit scale, so a crossing location carries an uncertainty
            # of se/|dP/ds|.  At n = 1500 that is 0.021 in scale -- the
            # same size as the effect being measured, which is exactly what
            # the first run showed.  n is set from a fixed event budget and
            # the resulting scale uncertainty is reported per point rather
            # than assumed small.
            n = int(min(40000, max(6000, 8e8 / L)))
            s_mc = critical_scale(shape, L, n, "oracle")
            if not np.isfinite(s_mc):
                print(f"  {name:40s} {L:7d}   no crossing of 0.95")
                continue
            iss_mc = e_iss(tuple(v * s_mc for v in shape))
            s_pred = m.boundary(shape, L, target=TARGET)
            iss_pred = (e_iss(tuple(v * s_pred for v in shape))
                        if np.isfinite(s_pred) else float("nan"))
            se = float(np.sqrt(TARGET * (1 - TARGET) / n))
            # |dP/ds| at the crossing from the analytic curve (the two
            # curves share a slope to leading order), propagated to a scale
            # uncertainty and then to Iss.
            _h = 0.004 * s_mc
            _slope = abs(m.pr_at(shape, s_mc + _h, L)[0]
                         - m.pr_at(shape, s_mc - _h, L)[0]) / (2 * _h)
            dsc = se / _slope if _slope > 0 else float("inf")
            iss_err = 0.5 * abs(
                e_iss(tuple(v * (s_mc + dsc) for v in shape))
                - e_iss(tuple(v * (s_mc - dsc) for v in shape)))
            rows.append(dict(shape=name, L=L, s_mc=s_mc, TL=sum(shape) * s_mc,
                             iss_mc=iss_mc, s_pred=s_pred, iss_pred=iss_pred,
                             diff=iss_mc - iss_pred, n_rep=n,
                             dPds=_slope, scale_se=dsc, iss_err=iss_err))
            print(f"  {name:40s} {L:7d} {s_mc:8.4f} {sum(shape) * s_mc:7.3f} "
                  f"{iss_mc:9.5f} {iss_err:7.5f} {iss_pred:11.5f} "
                  f"{iss_mc - iss_pred:+9.5f} {_slope:7.3f} {n:6d}",
                  flush=True)
            # written after EVERY point, not at the end: this sweep runs for
            # over an hour and a crash at the last point should not cost the
            # first nineteen
            pd.DataFrame(rows).to_csv(OUT / "issc_definition.csv", index=False)
        print()
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "issc_definition.csv", index=False)

    if len(d):
        print(f"  worst |Iss.c^MC - Iss.c^pred| = {d['diff'].abs().max():.5f}")
        print(f"  mean  |Iss.c^MC - Iss.c^pred| = {d['diff'].abs().mean():.5f}")
        print("\n  For reference, the whole range Iss.c spans here is")
        print(f"  {d[['iss_mc', 'iss_pred']].min().min():.5f} to "
              f"{d[['iss_mc', 'iss_pred']].max().max():.5f}.")
    print("\n  The comparison is made in Iss, not in probability, so it goes")
    print("  through the map TL -> Iss.  That map is very flat near")
    print("  saturation, which is where most of these points sit, so the")
    print("  Iss agreement should be BETTER than the probability agreement")
    print("  -- and if it is not, the map or the boundary is wrong.")


if __name__ == "__main__":
    main()
