"""Is the 1.4% residual an averaging artefact, or something unmodelled?

The window comparison is between two differently-averaged objects:

    p_fit  = the exponent of a straight line through log(gap) against log L
    p_pred = lam_Iss / (2 lam_Delta - lam_v2), each rate fitted over a scale
             window

A ratio of window averages is not the window average of the ratio, and the
local exponents are known to drift (0.352 down to 0.180 on one Farris tree),
so the two cannot agree exactly even if the mechanism is right.  The residual
therefore says almost nothing about whether the mechanism is right.

This script removes the averaging from both sides.  At each L it compares

    p_local(L) = -d log gap / d log L        measured pointwise, and
    ratio(L)   = lam_Iss / lam_Delta         both rates measured locally at
                                             the critical scale s*(L)

If these agree far better than the window-averaged versions do, the residual
is an averaging artefact and the mechanism statement is exact locally, which
is a far stronger claim than "1.4% is good enough".  If they do not, then
something is genuinely missing and the honest answer is that multiple hits
are not the whole story.

Everything uses the star-projection fast path validated in s36.
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

OUT = PART2 / "results" / "factors"

_spec = importlib.util.spec_from_file_location(
    "s36", str(HERE / "s36_factors.py"))
s36 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s36)

m = s36.m
ISS_SAT = s36.ISS_SAT


def gap(shape, s):
    return ISS_SAT - s36.e_iss(shape, s)


def local_rates(shape, s, h=1e-3):
    """(lam_Delta, lam_Iss) at the scale s, by central differences."""
    def ld(x):
        b1 = tuple(v * x for v in shape)
        p = m.pattern_batch(b1, m.SPLITS["12|34"])[0]
        fits = s36.fits_fast(b1)
        mu, cov = m.moments(p, m.logs_of(fits))
        m2, m3, v2, v3, c23 = m.params(mu, cov)
        return max(m2, 1e-300)
    d = (np.log(ld(s + h)) - np.log(ld(s - h))) / (2 * h)
    g = (np.log(gap(shape, s + h)) - np.log(gap(shape, s - h))) / (2 * h)
    return -d, -g                                  # both decay, so negate


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    shapes = {"reference a=b=1.0 c=d=0.05 e=0.02": s36.REF,
              "e x5      a=b=1.0 c=d=0.05 e=0.10":
                  (1.0, 1.0, 0.05, 0.05, 0.10),
              "c=d x4    a=b=1.0 c=d=0.20 e=0.02":
                  (1.0, 1.0, 0.20, 0.20, 0.02)}

    Ls = np.geomspace(40000, 250000, 7)
    print("Locally measured exponent against the locally measured ratio.")
    print("Both sides pointwise, so no averaging can account for a gap.\n")
    rows = []
    for name, shape in shapes.items():
        ss, gs = [], []
        for L in Ls:
            s = s36.boundary_fast(shape, L, target=0.95)
            if not np.isfinite(s):
                continue
            ss.append(s)
            gs.append(gap(shape, s))
        if len(gs) < 5:
            print(f"  {name}: too few points")
            continue
        ss = np.array(ss)
        gs = np.array(gs)
        print(f"  {name}")
        print(f"    {'L':>8} {'s*':>8} {'p_local':>9} {'ratio':>9} "
              f"{'rel err':>9} {'lam_D':>8} {'lam_I':>8}")
        rels = []
        for i in range(2, len(ss) - 2):
            # 5-point stencil on log(gap) against log(s*)
            k = np.arange(i - 2, i + 3)
            c = np.polyfit(np.log(ss[k]), np.log(gs[k]), 2)
            dlogs = np.polyval(np.polyder(c), np.log(ss[i]))
            # d log s* / d log L at this point, from a local line
            kk = np.arange(max(0, i - 1), min(len(ss), i + 2))
            cs = np.polyfit(np.log(Ls[kk]), np.log(ss[kk]), 1)
            p_local = -dlogs * cs[0]
            lD, lI = local_rates(shape, ss[i])
            ratio = lI / lD
            rel = (ratio - p_local) / p_local if p_local else np.nan
            rels.append(rel)
            rows.append(dict(config=name, L=float(Ls[i]), s=float(ss[i]),
                             p_local=float(p_local), ratio=float(ratio),
                             rel=float(rel), lam_D=float(lD), lam_I=float(lI)))
            print(f"    {Ls[i]:8.0f} {ss[i]:8.4f} {p_local:9.4f} {ratio:9.4f} "
                  f"{rel:+9.4f} {lD:8.4f} {lI:8.4f}")
        print(f"    mean |rel err| = {np.abs(rels).mean():.4f}\n")

    d = pd.DataFrame(rows)
    d.to_csv(OUT / "local_ratio.csv", index=False)
    if len(d):
        win = pd.read_csv(OUT / "factors.csv")
        w = win[win.rel_err.abs() < 0.2]
        print(f"  pointwise mean |rel err| : {d.rel.abs().mean():.4f}")
        print(f"  window   mean |rel err|  : {w.rel_err.abs().mean():.4f}")
        print("\n  If the pointwise figure is much the smaller, the window")
        print("  residual was averaging, and the mechanism statement holds")
        print("  exactly in the local limit.")


if __name__ == "__main__":
    main()
