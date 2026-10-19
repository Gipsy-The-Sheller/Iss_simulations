"""The admissible window in the three-nucleotide weight, as r grows.

CASTER's supplementary figure S1 scans the weight w3 against a branch-length
ratio 1 <= r < 12 and finds that -1/2 is the only weight leaving both the
Farris and the Felsenstein zone margins positive.  This recomputes the same
object exactly and extends the scan past 12, which is where the interesting
behaviour is.

For each r there are two boundaries.  The Farris margin falls through zero at
some w3 and the Felsenstein margin at another; the admissible set is the
interval between them, and a weight works only if it lies inside it.  The
bounds are found by scanning w3 and locating the sign changes, which is
robust when the margin is flat; a root finder was tried first and walked off
the end of the bracket wherever the margin had no crossing at all.

The ratio r is the same one Fig. S1 uses: the long branches carry length
0.1 r, the short ones 0.1, and the internal branch 0.1, so r = 1 is a tree
with every branch equal and r is the long-to-short ratio.

Outputs
    results/caster_weights/window.csv
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

OUT = PART2 / "results" / "caster_weights"

_spec = importlib.util.spec_from_file_location("s49", str(HERE / "s49_caster_weights.py"))
s49 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s49)
m = s49.m

W_GRID = np.round(np.linspace(-2.0, 0.5, 1001), 6)
R_GRID = np.unique(np.round(np.concatenate([
    np.geomspace(1.0, 200.0, 60), [4.0, 8.0, 12.0, 16.0, 20.0, 30.0]]), 4))


def margin(w3, shape):
    V = s49.vote_table(w3)
    p = m.pattern_batch(shape, m.SPLITS["12|34"])[0]
    return float((V[:, 0] - V[:, 1]) @ p)


def admissible(shape):
    """Bounds of the set of w3 where this zone's margin is positive."""
    vals = np.array([margin(w, shape) for w in W_GRID])
    ok = vals > 0
    if not ok.any():
        return None, None
    return float(W_GRID[ok][0]), float(W_GRID[ok][-1])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for r in R_GRID:
        far = (0.1 * r, 0.1 * r, 0.1, 0.1, 0.1)     # long branches are sisters
        fel = (0.1 * r, 0.1, 0.1 * r, 0.1, 0.1)     # long branches are not
        fa, fb = admissible(far)
        la, lb = admissible(fel)
        if None in (fa, fb, la, lb):
            lo = hi = np.nan
        else:
            lo, hi = max(fa, la), min(fb, lb)       # both must hold
        if not np.isfinite(lo):
            kind = "undefined"
        elif hi < lo - 1e-9:
            kind = "empty"
        elif hi - lo < 1e-9:
            kind = "point"
        else:
            kind = "interval"
        rows.append(dict(r=float(r), farris_lo=fa, farris_hi=fb,
                         felsen_lo=la, felsen_hi=lb,
                         window_lo=lo, window_hi=hi, kind=kind,
                         width=(hi - lo) if kind == "interval" else 0.0,
                         margin_w0=margin(0.0, far),
                         margin_whalf=margin(-0.5, far)))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "window.csv", index=False)

    print("Admissible w3 window against the branch-length ratio r\n")
    print(f"  {'r':>7} {'Farris allows':>20} {'Felsenstein allows':>22} "
          f"{'window':>20}")
    for _, x in d.iterrows():
        if x.kind == "interval":
            w = f"({x.window_lo:+.4f}, {x.window_hi:+.4f})"
        elif x.kind == "point":
            w = f"{{{x.window_lo:+.4f}}}"
        else:
            w = x.kind
        print(f"  {x.r:7.2f} ({x.farris_lo:+.3f}, {x.farris_hi:+.3f})"
              f"{'':>4} ({x.felsen_lo:+.3f}, {x.felsen_hi:+.3f})"
              f"{'':>3} {w:>20}")

    fin = d[(d.kind == "interval") & np.isfinite(d.width)]
    if len(fin):
        print(f"\n  widest  w = {fin.width.max():.4f} at r = "
              f"{fin.loc[fin.width.idxmax(), 'r']:.2f}")
        thin = fin[fin.width < 1e-9]
        if len(thin):
            print(f"  width falls below the 0.0025 grid step by r = "
                  f"{thin.r.iloc[0]:.2f}")
    print("\n  Beyond r about 20 both margins are at the 1e-18 level, where the")
    print("  sign is floating-point noise and the bounds mean nothing.  The")
    print("  last row at which the window is a genuine interval is the")
    print("  statement worth quoting, not the closure point.")


if __name__ == "__main__":
    main()
