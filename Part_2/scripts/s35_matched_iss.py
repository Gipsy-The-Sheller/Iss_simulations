"""Same Iss, two trees, different power.

The objection this answers: the Farris zone is hard because of substitution
saturation, and Iss measures saturation, so nothing new is being said.

Take it literally and it is testable.  If Iss is the causally relevant
saturation variable, then two trees with the SAME Iss must have the same
recoverability, whatever their shapes.  So: for a range of Iss values, find
a Farris-zone tree and a balanced tree that share that Iss, and compare them.

Iss is a monotone function of the tree scale within a shape, so matching is a
scalar root-find, and everything below is exact pattern enumeration -- no
simulation, no fitted parameters.

The companion construction (s34) is stronger but blunter: it holds Iss
exactly fixed while power sweeps its whole range, using data that are not
trees at all.  This one stays inside the tree family, so it answers the
objection on its own ground -- both data sets are genuine trees growing
under the same model, and one is hard for reasons that Iss does not see.

The quantity compared is z = Delta sqrt(L) / sqrt(v2), the standardised
signal the decision rests on.  Pr[correct] follows from it and is reported
alongside.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
PART2 = HERE.parent
ROOT = PART2.parent
sys.path.insert(0, str(ROOT / "src"))

OUT = PART2 / "results" / "matched_iss"

_spec = importlib.util.spec_from_file_location(
    "s27", str(HERE / "s27_mc_boundary.py"))
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)

from issbench import exact, trees                                  # noqa: E402
from issbench.models import build_Q                                 # noqa: E402

Q = build_Q("JC69")
PI = np.ones(4) / 4

FARRIS = (1.0, 1.0, 0.05, 0.05, 0.02)
BALANCED = (0.1, 0.1, 0.1, 0.1, 0.2)


def e_iss(shape, s):
    return float(exact.E_Iss(4, trees.quartet(*[v * s for v in shape]),
                             Q, PI, convention="A"))


def z_at(shape, s, L):
    """Standardised signal z, with all three topologies fitted at KL."""
    b1 = tuple(v * s for v in shape)
    p = m.pattern_batch(b1, m.SPLITS["12|34"])[0]
    fits = []
    for lab in m.SPLITS:
        _, bb = m.best_fit(p, m.SPLITS[lab], b1)
        fits.append(m.pattern_batch(bb, m.SPLITS[lab])[0])
    mu, cov = m.moments(p, m.logs_of(fits))
    m2, m3, v2, v3, c23 = m.params(mu, cov)
    if v2 <= 0:
        return float("nan"), float("nan")
    return float(m2 * np.sqrt(L) / np.sqrt(v2)), m.analytic(p, fits, L)


def match(s_bal_target, lo=0.01, hi=40.0):
    """Scale the balanced tree so its Iss equals the target."""
    f = lambda s: e_iss(BALANCED, s) - s_bal_target
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return float("nan")
    return float(brentq(f, lo, hi, xtol=1e-10))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    L = 30000

    print("Same Iss, two trees, different power.  L = 30000.")
    print("Both trees are genuine; only the shape differs.  If Iss were the")
    print("causal saturation variable, z would have to agree.\n")
    print(f"  {'Iss':>8} {'Farris TL':>10} {'z_Farris':>10} "
          f"{'bal TL':>9} {'z_bal':>9} {'ratio':>8} {'Pr_F':>7} {'Pr_B':>7}")

    rows = []
    for s_f in (0.15, 0.25, 0.40, 0.60, 0.90, 1.30):
        iss = e_iss(FARRIS, s_f)
        s_b = match(iss)
        if not np.isfinite(s_b):
            continue
        zf, pf = z_at(FARRIS, s_f, L)
        zb, pb = z_at(BALANCED, s_b, L)
        ratio = zf / zb if zb and np.isfinite(zb) and zb != 0 else float("nan")
        rows.append(dict(iss=iss, s_farris=s_f, TL_farris=sum(FARRIS) * s_f,
                         z_farris=zf, p_farris=pf, s_balanced=s_b,
                         TL_balanced=sum(BALANCED) * s_b, z_balanced=zb,
                         p_balanced=pb, ratio=ratio))
        print(f"  {iss:8.5f} {sum(FARRIS) * s_f:10.3f} {zf:10.4f} "
              f"{sum(BALANCED) * s_b:9.3f} {zb:9.4f} {ratio:8.3f} "
              f"{pf:7.4f} {pb:7.4f}", flush=True)

    d = pd.DataFrame(rows)
    d.to_csv(OUT / "matched_iss.csv", index=False)

    if len(d):
        print(f"\n  at matched Iss, z differs by a factor of up to "
              f"{d.ratio.max():.2f}")
        print(f"  Pr[correct] ranges {d.p_farris.min():.4f}-{d.p_farris.max():.4f}"
              f" on the Farris tree and {d.p_balanced.min():.4f}-"
              f"{d.p_balanced.max():.4f} on the balanced one")
        print("\n  The two data sets are real trees under the same model and")
        print("  the same Iss.  Whatever makes one harder is not what Iss")
        print("  measures, so 'the Farris zone is hard because of saturation'")
        print("  cannot be read as 'because of Iss'.")
    print("\n  Caveat carried deliberately: z is not monotone in the tree")
    print("  scale, so a ratio of z is not a ratio of difficulty at every")
    print("  point.  The claim made here is the weaker, sufficient one --")
    print("  matched Iss does not imply matched z.")


if __name__ == "__main__":
    main()
