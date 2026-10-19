"""Re-run s26's scan: were its counterexamples real?

Commit c0894c1 records that the star-projection claim is false in general,
on the strength of s26's scan, which flagged 11 of 230 projections as
landing off the star.  Two of those configurations were re-examined by hand
and both turned out to have a KL identical to the star fit's to 1e-14 --
so the "interior projections" may be the optimiser wandering along a flat
direction rather than finding a genuinely better point.

That distinction decides whether the claim is false or exact, so the scan is
repeated here with three changes, none of which is a matter of taste:

  * the objective is compared, not the parameter.  A projection is interior
    only if it beats the star in KL.  Where the two agree the configuration
    is recorded as flat, whatever branch length the optimiser returned.
  * the optimiser is warm-started from the star fit, which is where the
    answer is when the claim holds, so a spurious interior point has to be
    found *from* the star rather than merely being reachable.
  * the domain is widened: the old scan took e_true <= 3, and if the star
    projection fails anywhere it fails at saturation, so e_true is taken out
    to 12 here.

The second point matters most.  Starting from a distant point and descending
into a flat valley lands anywhere; starting from the star and failing to
improve on it is evidence.
"""
from __future__ import annotations

import importlib.util
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "results" / "star_rescan"

_spec = importlib.util.spec_from_file_location(
    "s27", str(HERE / "s27_mc_boundary.py"))
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)

ALT = m.SPLITS["13|24"]
TRUE = m.SPLITS["12|34"]


def star_fit(p, b1):
    """The best STAR on ALT, with e held at exactly zero.

    It is not enough to take the unconstrained fit and set its internal
    branch to zero: the pendants of that fit were chosen jointly with a
    nonzero internal branch and are not the star's pendants.  Doing it that
    way inflates KL_star and makes ordinary configurations look like
    counterexamples.  The four pendants are fitted here on their own.
    """
    x0 = np.log([b1[0], b1[2], b1[1], b1[3]])

    def neg(par):
        th = np.exp(par)
        q = m.pattern_batch((th[0], th[1], th[2], th[3], 0.0), ALT)[0]
        return float(np.sum(p * np.log(np.maximum(p, 1e-300)
                                       / np.maximum(q, 1e-300))))

    best = (np.inf, x0)
    for t in (x0, x0 * 0.3, np.log([0.5, 0.5, 0.5, 0.5]),
              np.log([1.0, 0.1, 1.0, 0.1])):
        r = minimize(neg, np.asarray(t, float), method="L-BFGS-B",
                     bounds=[(-20.0, 5.0)] * 4,
                     options=dict(ftol=1e-20, gtol=1e-14, maxiter=800))
        if r.fun < best[0]:
            best = (float(r.fun), np.asarray(r.x, float))
    th = np.exp(best[1])
    return np.array([th[0], th[1], th[2], th[3], 0.0]), best[0]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261015)
    t0 = time.time()

    cases = [(*rng.uniform(0.02, 1.5, 4), float(rng.uniform(1e-4, 12.0)))
             for _ in range(120)]
    for e in (1e-4, 0.01, 0.1, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0):
        for pend in ((1.0, 0.05, 1.0, 0.05), (1.0, 0.05, 0.7, 0.2),
                     (0.8, 0.3, 0.5, 0.1), (1.5, 1.2, 0.08, 0.05),
                     (0.5, 0.5, 0.5, 0.5), (1.427, 1.266, 1.109, 0.423),
                     (0.05, 0.05, 1.0, 1.0), (0.26, 0.335, 1.292, 0.997),
                     (0.708, 0.975, 0.07, 0.22), (0.373, 0.178, 0.369, 1.443),
                     (1.5, 1.2, 0.08, 0.05)):
            cases.append((*pend, e))

    print(f"re-scanning {len(cases)} configurations, e_true out to 12")
    print(f"    {'case':>5} {'e_true':>8} {'KL_star':>12} {'KL_best':>12} "
          f"{'KL_star-KL_best':>16} {'e_fit':>8} verdict")
    rows = []
    for n, b1 in enumerate(cases):
        p = m.pattern_batch(b1, TRUE)[0]
        star, kl_star = star_fit(p, b1)
        # the unconstrained fit, started FROM the star: a counterexample has
        # to be reached from the star, not merely be reachable
        _, bb = m.best_fit(p, ALT, star)
        kl_best = m.kl(p, m.pattern_batch(bb, ALT)[0])
        drop = kl_star - kl_best
        tol = 1e-6 * max(1e-12, abs(kl_star))
        if drop > tol and bb[4] > 1e-4:
            verdict = "GENUINE INTERIOR"
        elif bb[4] > 0.02:
            verdict = "flat, wandered"
        else:
            verdict = "star"
        rows.append(dict(case=n, a=b1[0], b=b1[1], c=b1[2], d=b1[3],
                         e_true=b1[4], kl_star=kl_star, kl_best=kl_best,
                         drop=drop, e_fit=float(bb[4]), verdict=verdict))
        if verdict != "star" or (n % 25 == 0):
            print(f"    {n:5d} {b1[4]:8.4f} {kl_star:12.5e} "
                  f"{kl_best:12.5e} {drop:+16.3e} {bb[4]:8.4f} {verdict}",
                  flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "rescan.csv", index=False)

    print(f"\n  {time.time() - t0:.0f}s")
    print("\n  verdict counts")
    print(d.verdict.value_counts().to_string())
    gen = d[d.verdict == "GENUINE INTERIOR"]
    if len(gen):
        print(f"\n  genuine interior projections: {len(gen)} of {len(d)}")
        print(gen.nlargest(12, "drop")[
            ["a", "b", "c", "d", "e_true", "drop", "e_fit"]].to_string(index=False))
    else:
        print("\n  no configuration in this scan has a projection that beats")
        print("  the star fit in KL.  Every interior branch length the earlier")
        print("  scan reported is an optimiser artefact on a flat objective.")
    fl = d[d.verdict == "flat, wandered"]
    if len(fl):
        print(f"\n  flat directions with a nonzero fitted branch: {len(fl)}, "
              f"largest {fl.e_fit.max():.4f} at e_true = "
              f"{fl.loc[fl.e_fit.idxmax(), 'e_true']:.3f}")
        print(f"  their KL gap runs to {fl['drop'].abs().max():.2e}")


if __name__ == "__main__":
    main()
