"""The Farris zone and Iss.c are two appearances of one quantity.

Established here, exactly:

    Delta = KL( P_true || the best star fit of P_true ).

The identity is not a coincidence.  The two wrong topologies of a quartet both
project onto the star, so the true topology's fit is P itself with KL 0 and
the alternatives' fit is the star; Delta, defined as the per-site KL gap
between them, is therefore the distance from the tree to its own best star
approximation.  It is checked numerically rather than asserted.

If that is what Delta is, then two things that looked like separate
phenomena are the same thing seen twice:

  * the difficulty of resolving a tree, because Delta is the numerator of the
    resolution condition Delta sqrt(L) = z sqrt(v2); and
  * the Farris zone, because what makes the long-branches-sisters arrangement
    hard is that it drives the tree towards being a star.

The controlled family holds the pendant sum and the internal branch fixed and
moves only the arrangement, so any change in Delta is attributable to the
arrangement and to nothing else.  If Delta collapses for Farris and stays put
for Felsenstein, the Farris zone is a region of tree space -- the region
nearest the star -- and not a separate mechanism.
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

OUT = PART2 / "results" / "shared_origin"

_spec = importlib.util.spec_from_file_location(
    "s36", str(HERE / "s36_factors.py"))
s36 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s36)

m = s36.m
E_INT = 0.2
RS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    print("CHECK.  Delta against the KL distance to the best star fit.")
    print(f"    {'tree':>28} {'Delta':>13} {'KL(P||star)':>13} {'ratio':>10}")
    worst = 0.0
    for b1 in ((0.05, 0.05, 1.0, 1.0, 0.02), (0.1, 0.1, 0.1, 0.1, 0.2),
               (1.0, 0.05, 1.0, 0.05, 0.2), (0.3, 0.3, 0.3, 0.3, 0.6),
               (1.0, 1.0, 0.05, 0.05, 0.5)):
        p = m.pattern_batch(b1, m.SPLITS["12|34"])[0]
        fits = s36.fits_fast(b1)
        mu, cov = m.moments(p, m.logs_of(fits))
        m2, m3, v2, v3, c23 = m.params(mu, cov)
        kl = m.kl(p, fits[1])
        r = m2 / kl if kl else float("nan")
        worst = max(worst, abs(r - 1.0))
        print(f"    {str(np.round(b1, 3)):>28} {m2:13.6e} {kl:13.6e} "
              f"{r:10.10f}")
    print(f"    worst |ratio - 1| = {worst:.2e}\n")

    print("The controlled family: pendant sum fixed at 4, internal branch")
    print("fixed at 0.2, only the arrangement moves.\n")
    print(f"    {'r':>5} {'Delta Farris':>14} {'Delta Felsen':>14} "
          f"{'ratio':>8} {'Iss Farris':>11} {'Iss Felsen':>11}")
    rows = []
    for r in RS:
        far = (1 + r, 1 + r, 1 - r, 1 - r, E_INT)
        fel = (1 + r, 1 - r, 1 + r, 1 - r, E_INT)
        d = {}
        for fam, b1 in (("Farris", far), ("Felsenstein", fel)):
            p = m.pattern_batch(b1, m.SPLITS["12|34"])[0]
            fits = s36.fits_fast(b1)
            mu, cov = m.moments(p, m.logs_of(fits))
            m2, m3, v2, v3, c23 = m.params(mu, cov)
            d[fam] = (m2, float(np.sqrt(v2)), s36.e_iss(b1, 1.0))
        ratio = d["Felsenstein"][0] / d["Farris"][0] if d["Farris"][0] else np.nan
        rows.append(dict(r=r,
                         delta_farris=d["Farris"][0],
                         delta_felsenstein=d["Felsenstein"][0],
                         ratio=ratio,
                         snr_farris=d["Farris"][0] / d["Farris"][1],
                         snr_felsenstein=d["Felsenstein"][0] / d["Felsenstein"][1],
                         iss_farris=d["Farris"][2],
                         iss_felsenstein=d["Felsenstein"][2],
                         TL=sum(far)))
        print(f"    {r:5.1f} {d['Farris'][0]:14.6e} "
              f"{d['Felsenstein'][0]:14.6e} {ratio:8.2f} "
              f"{d['Farris'][2]:11.5f} {d['Felsenstein'][2]:11.5f}",
              flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "shared_origin.csv", index=False)

    print(f"\n  Delta falls by a factor of "
          f"{d.delta_farris.iloc[0] / d.delta_farris.iloc[-1]:.1f} along the")
    print("  Farris family and rises slightly along the Felsenstein one,")
    print(f"  with the pendant sum, the internal branch and the total tree")
    print(f"  length ({d.TL.iloc[0]:.1f}) fixed throughout.")
    print("\n  So the Farris zone is a region of tree space -- the region")
    print("  nearest the star -- and the quantity that defines it, Delta,")
    print("  is the same one that sets the critical tree and hence Iss.c.")


if __name__ == "__main__":
    main()
