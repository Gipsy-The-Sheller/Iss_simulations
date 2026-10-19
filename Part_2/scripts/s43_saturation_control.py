"""Is the Farris effect just geometry-mediated saturation?

The objection: fixing the pendant sum and the internal branch fixes the TOTAL
amount of multiple substitution, but not its distribution across branches.
The Farris arrangement concentrates the length on two sister branches, so
local multiple substitution is heavier there, and the reviewer concludes that
Delta is only measuring saturation under another name.

The objection has a testable consequence, and it fails that test.  At every r
the two arrangements have the SAME MULTISET of branch lengths:

    Farris(r)       (1+r, 1+r, 1-r, 1-r, 0.2)
    Felsenstein(r)  (1+r, 1-r, 1+r, 1-r, 0.2)

so every per-branch substitution statistic is identical between them.
Substitutions on different branches are conditionally independent given the
tree, so anything that is a function of the per-branch counts alone -- the
expected number of substitutions on each branch, the probability that a
branch carries two or more, the probability that a site carries a multiple
hit somewhere -- is identical.  Not approximately: identically, as a
multiset of branch lengths.

If Delta differs by a factor of 36 between two configurations whose local
saturation profile is the same multiset, then local saturation cannot be the
variable that sets Delta.  What differs is the wiring, and the wiring is what
Delta measures.

This script computes the saturation statistics explicitly so the claim is a
number rather than an argument, and reports them beside Delta.
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

OUT = PART2 / "results" / "saturation_control"

_spec = importlib.util.spec_from_file_location(
    "s36", str(HERE / "s36_factors.py"))
s36 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s36)

m = s36.m
E_INT = 0.2
RS = (0.0, 0.2, 0.4, 0.6, 0.8)


def branch_stats(t):
    """(P(>=2 hits), P(no hit)) on a branch of expected length t, JC69.

    The branch length is the expected number of substitutions per site, so
    the count is Poisson(t) and P(0) = e^-t, P(>=2) = 1 - e^-t (1 + t).
    """
    return 1.0 - np.exp(-t) * (1.0 + t), float(np.exp(-t))


def site_multiple(ts):
    """P(at least one branch of the quartet carries a multiple hit)."""
    return 1.0 - np.prod([branch_stats(t)[1] + np.exp(-t) * t for t in ts])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    print("The saturation profile is a multisets property; Delta is not.")
    print("Pendant sum fixed at 4, internal branch at 0.2.\n")
    print(f"  {'r':>5} {'arrangement':>13} {'sorted branches':>34} "
          f"{'P(>=2) per branch':>26} {'P(site has a mult. hit)':>24} "
          f"{'Delta':>12}")
    rows = []
    for r in RS:
        for fam in ("Farris", "Felsenstein"):
            b1 = ((1 + r, 1 + r, 1 - r, 1 - r, E_INT) if fam == "Farris"
                  else (1 + r, 1 - r, 1 + r, 1 - r, E_INT))
            p = m.pattern_batch(b1, m.SPLITS["12|34"])[0]
            fits = s36.fits_fast(b1)
            mu, cov = m.moments(p, m.logs_of(fits))
            m2, m3, v2, v3, c23 = m.params(mu, cov)
            srt = np.sort(np.array(b1))[::-1]
            p2 = np.sort([branch_stats(t)[0] for t in b1])[::-1]
            psite = site_multiple(b1)
            rows.append(dict(r=r, family=fam, TL=float(sum(b1)),
                             branches=",".join(f"{v:.1f}" for v in srt),
                             p_ge2=",".join(f"{v:.4f}" for v in p2),
                             p_site_multiple=psite, delta=float(m2)))
            print(f"  {r:5.1f} {fam:>13} "
                  f"{','.join(f'{v:.1f}' for v in srt):>34} "
                  f"{','.join(f'{v:.4f}' for v in p2):>26} "
                  f"{psite:24.6f} {m2:12.6e}", flush=True)
        print()
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "saturation_control.csv", index=False)

    print("  At every r the two arrangements have the same sorted branch")
    print("  multiset, hence the same per-branch multiple-hit probabilities")
    print("  and the same probability that a site carries a multiple hit")
    print("  somewhere.  Delta differs nonetheless:\n")
    print(f"    {'r':>5} {'P(site mult. hit)':>19} {'ratio of Delta':>16}")
    for r in RS:
        a = d[(d.r == r) & (d.family == "Farris")].iloc[0]
        b = d[(d.r == r) & (d.family == "Felsenstein")].iloc[0]
        same = abs(a.p_site_multiple - b.p_site_multiple) < 1e-12
        print(f"    {r:5.1f} {a.p_site_multiple:19.6f} "
              f"{b.delta / a.delta:16.2f}   {'identical' if same else 'DIFFER'}")


if __name__ == "__main__":
    main()
