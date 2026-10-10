"""Theorem 13: is Psi affine in p-bar, and how badly is it not?

The theorem is a statement about Psi and the mean p-distance, both of which
are exact linear functionals of the spectral coordinates.  Testing it on
E[Iss] instead would confound the criterion with the O(3) remainder of the
entropy expansion, and the first version of this script did exactly that: the
criterion-holding model came out with the LARGEST residual, which is the
signature of a test that is measuring the wrong thing.

So the primary test is on Psi against p-bar, computed from the spectral
decomposition with no entropy and no simulation.  The residual from a
best-fit line is reported in Psi units and, after scaling by
(N-1)/(2 N h(pi)), in units of the index itself, which is what a reader can
compare with Remark 18.

The secondary trace gives the exact E[Iss] at N = 8 and 16 on star trees,
where the closed form makes it available, so that the criterion's magnitude is
visible where the paper actually makes its claim.

Three models, chosen to exercise all three corollaries:

  HKY, uniform composition   criterion holds  (Corollary uniform)
  F81, skewed composition    criterion holds  (Corollary degenerate)
  HKY, skewed composition    criterion FAILS  (Corollary failure)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench import exact, trees                                    # noqa: E402
from issbench.entropy import h                                       # noqa: E402
from issbench.models import build_Q, P_matrix                        # noqa: E402
from issbench.spectral import factorize                              # noqa: E402

OUT = ROOT / "results" / "theorem13"
UNIFORM = np.ones(4) / 4
SKEW = np.array([0.40, 0.30, 0.20, 0.10])
CASES = [
    ("HKY, uniform", build_Q("HKY85", kappa=5.0, pi=UNIFORM), UNIFORM),
    ("F81, skewed", build_Q("F81", pi=SKEW), SKEW),
    ("HKY, skewed", build_Q("HKY85", kappa=5.0, pi=SKEW), SKEW),
]
D_GRID = np.exp(np.linspace(np.log(0.02), np.log(1.5), 120))


def star_trace(Q, pi, N):
    """Psi, p-bar and exact E[Iss] along a star of root-to-leaf length d."""
    sp = factorize(Q, pi)
    rows = []
    for d in D_GRID:
        rows.append(dict(
            d=float(d), psi=float(sp.Psi(d)), pbar=float(sp.w_p(d)),
            iss=exact.E_Iss(N, trees.star(N, float(d)), Q, pi, convention="A")))
    return pd.DataFrame(rows)


def depart(t, ycol: str) -> dict:
    """Departure from the best-fit straight line in p-bar."""
    x, y = t.pbar.to_numpy(), t[ycol].to_numpy()
    b, a = np.polyfit(x, y, 1)
    r = y - (a + b * x)
    span = float(y.max() - y.min())
    return dict(slope=float(b), intercept=float(a),
                max_abs=float(np.abs(r).max()),
                as_pct_of_span=float(100 * np.abs(r).max() / span))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    crit, prim, sec = [], [], []
    for name, Q, pi in CASES:
        sp = factorize(Q, pi)
        ok = sp.affine_ok()                       # it is a method, not a property
        ratios = [round(float(r), 3) for r in sp.ratios]
        crit.append(dict(model=name, ratios=str(ratios), affine=ok))
        print(f"{name:14s} tau/w = {ratios}  ->  "
              f"{'criterion HOLDS' if ok else 'criterion FAILS'}")

        t = star_trace(Q, pi, 8)
        dep = depart(t, "psi")
        # a residual in Psi moves the index by (N-1)/(2 N h(pi)); report at N=16
        scale16 = (16 - 1) / (2 * 16 * h(pi))
        prim.append(dict(model=name, affine=ok, **dep,
                         max_in_index_units=dep["max_abs"] * scale16))
        print(f"   Psi vs p-bar: slope {dep['slope']:.4f}  "
              f"max |resid| {dep['max_abs']:.2e} Psi "
              f"= {dep['max_abs'] * scale16:.4f} in index units "
              f"({dep['as_pct_of_span']:.3f}% of the range of Psi)")
        t.to_csv(OUT / f"trace_{name.split(',')[0]}_{'uni' if 'uniform' in name else 'skew'}.csv",
                 index=False)

    print()
    for N in (8, 16):
        for name, Q, pi in CASES:
            t = star_trace(Q, pi, N)
            dep = depart(t, "iss")
            sec.append(dict(N=N, model=name, **dep))
            print(f"N={N:2d}  {name:14s} E[Iss] vs p-bar: slope {dep['slope']:.4f}"
                  f"  max |resid| {dep['max_abs']:.2e} "
                  f"({dep['as_pct_of_span']:.2f}% of range)")
            t.to_csv(OUT / f"iss_N{N}_{name.replace(', ', '_').replace(' ', '')}.csv",
                     index=False)

    pd.DataFrame(crit).to_csv(OUT / "criterion.csv", index=False)
    pd.DataFrame(prim).to_csv(OUT / "psi_affinity.csv", index=False)
    pd.DataFrame(sec).to_csv(OUT / "iss_affinity.csv", index=False)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
