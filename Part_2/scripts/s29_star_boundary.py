"""Where does the star projection stop holding?

s26 found numerically that the KL projection of a quartet onto a wrong
topology lands on the star, and that it stops doing so once the true internal
branch gets large.  That is a phase boundary in a five-dimensional space and
scanning it is not an explanation.  There is an exact local criterion.

Write F(theta) = KL(P || Q_theta) on the alternative topology, theta =
(pendants, e).  At the star point theta* the pendant gradient vanishes, so by
the envelope theorem the derivative of F along the *optimal* pendant path in
e is just the partial derivative at theta*:

    d/de F(theta(e))|_0 = dF/de|_theta*

and since e >= 0 is a boundary of the parameter set, e = 0 is a local optimum
exactly when that derivative is non-negative.  So there is one number per
configuration, computable in closed form:

    dQ/dx = cross - ssab (x) sscd / 4,   dQ/de = -(4/3) x dQ/dx,
    dKL/de = -sum_x P(x) (1/Q(x)) dQ/de,   at x = 1.

Two things then need separating, and the first version of this script
conflated them.

  * LOCAL.  Is the star a local minimum of the alternative's objective?
    That is the criterion above, and its zero is a curve e_c(pendants).
  * GLOBAL.  Is the star the *global* minimum?  A local optimum elsewhere
    can beat it, and nothing about a derivative at the star sees that.

So the script reports both on the same grid: the criterion, and the internal
branch of the full projection as found by s26's optimiser.  Where they
disagree the boundary is not the criterion's zero, and "e <= 2" from the
scan is not the criterion either.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "results" / "star_boundary"

_spec = importlib.util.spec_from_file_location(
    "s27", str(HERE / "s27_mc_boundary.py"))
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)

SPLITS = m.SPLITS
PI = m.PI
ALT = SPLITS["13|24"]          # labels ("A", "C", "B", "D")
TRUE = SPLITS["12|34"]


def pieces(lengths, labels):
    """(Q, dQ/dx_e) as 256-vectors, by the same algebra as pattern_batch."""
    th = np.atleast_2d(np.asarray(lengths, float))
    x = np.exp(-4.0 * th / 3.0)
    xa, xb, xc, xd, xe = x.T
    alpha, beta = (1.0 - xe) / 4.0, xe
    st = np.arange(4)

    def M(xi):
        return ((1.0 - xi)[:, None, None] / 4.0
                + xi[:, None, None] * (st[None, None, :] == st[None, :, None]))

    Ma, Mb, Mc, Md = M(xa), M(xb), M(xc), M(xd)
    ssab = np.einsum("bis,bjs,s->bij", Ma, Mb, PI)
    sscd = np.einsum("bks,bls->bkl", Mc, Md)
    cross = np.einsum("bis,bjs,bks,bls,s->bijkl", Ma, Mb, Mc, Md, PI)
    outer = ssab[:, :, :, None, None] * sscd[:, None, None, :, :]
    q = (alpha.reshape(-1, 1, 1, 1, 1) * outer
         + beta.reshape(-1, 1, 1, 1, 1) * cross)
    dq = cross - outer / 4.0
    pack = m._PACK[labels]
    qq = np.empty((q.shape[0], 256))
    dd = np.empty((dq.shape[0], 256))
    qq[:, pack] = q.reshape(-1, 256)
    dd[:, pack] = dq.reshape(-1, 256)
    return qq, dd


def pattern_true(a, b, c, d, e, labels=TRUE):
    return m.pattern_batch((a, b, c, d, e), labels)[0]


def star_fit(p, b1):
    """The best star on ALT: four pendants, no internal branch.

    General rather than symmetry-reduced.  Reducing the fit to two parameters
    was the first version's mistake: the reduction is only valid when the
    true quartet has the matching pendant symmetries, and imposing it where
    it does not hold makes the criterion meaningless rather than merely
    imprecise.
    """
    a, b, c, d = b1[:4]

    def neg(par):
        th = np.exp(par)
        q = m.pattern_batch((th[0], th[1], th[2], th[3], 0.0), ALT)[0]
        return float(np.sum(p * np.log(np.maximum(p, 1e-300)
                                       / np.maximum(q, 1e-300))))

    best = (np.inf, np.log([a, c, b, d]))
    for t in (np.log([a, c, b, d]), np.log([a, c, b, d]) * 0.3,
              np.log([0.5, 0.5, 0.5, 0.5]), np.log([1.0, 0.1, 1.0, 0.1])):
        r = minimize(neg, np.asarray(t, float), method="L-BFGS-B",
                     bounds=[(-20.0, 5.0)] * 4,
                     options=dict(ftol=1e-20, gtol=1e-14, maxiter=800))
        r = minimize(neg, r.x, method="Nelder-Mead",
                     options=dict(maxiter=6000, xatol=1e-14, fatol=1e-20))
        if r.fun < best[0]:
            best = (float(r.fun), np.asarray(r.x, float))
    return np.exp(best[1]), best[0]


def criterion(p, pars):
    """dKL/de of the alternative at the star point, e = 0.

    Non-negative means e = 0 is a LOCAL optimum of the alternative's
    objective.  It says nothing about whether a better point exists further
    out, which is the second thing the grid below measures.
    """
    q, dq = pieces((pars[0], pars[1], pars[2], pars[3], 0.0), ALT)
    q, dq = q[0], dq[0]
    dqde = -(4.0 / 3.0) * dq                     # x_e = 1 at e = 0
    return float(-np.sum(p * dqde / np.maximum(q, 1e-300)))


SHAPES = {"Farris (long sisters)":   (0.05, 0.05, 1.00, 1.00),
          "Felsenstein (apart)":     (1.00, 0.05, 1.00, 0.05),
          "generic":                 (0.80, 0.30, 0.50, 0.10),
          "near-equal pendants":     (0.50, 0.50, 0.50, 0.50),
          "one long pendant":        (1.00, 0.10, 0.10, 0.10)}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    print("CHECK.  at e = 0 the true tree IS the star, the star fit is exact")
    print("and the criterion must vanish.")
    for name, (a, b, c, d) in SHAPES.items():
        p = pattern_true(a, b, c, d, 0.0)
        pars, kl = star_fit(p, (a, b, c, d, 0.0))
        print(f"    {name:22s} KL={kl:+.2e}  dKL/de={criterion(p, pars):+.2e}")

    print("\nLOCAL criterion against the GLOBAL projection, on a grid.")
    print("  criterion >= 0 means the star is a local optimum; e_fit is what")
    print("  the full optimiser finds.  They agree only where e_fit = 0.")
    rows = []
    for name, (a, b, c, d) in SHAPES.items():
        print(f"\n  {name}  a,b,c,d = {a}, {b}, {c}, {d}")
        print(f"    {'e_true':>7} {'criterion':>12} {'star KL':>11} "
              f"{'e_fit':>8} {'KL(fit)':>11} {'agree':>6}")
        for e_t in (0.1, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0):
            p = pattern_true(a, b, c, d, e_t)
            pars, kl_star = star_fit(p, (a, b, c, d, e_t))
            crit = criterion(p, pars)
            _, bb = m.best_fit(p, ALT, (a, c, b, d, 0.1))
            kl_full = m.kl(p, m.pattern_batch(bb, ALT)[0])
            ok = (crit >= 0) == (bb[4] < 1e-4)
            rows.append(dict(shape=name, a=a, b=b, c=c, d=d, e_true=e_t,
                             criterion=crit, star_kl=kl_star, e_fit=bb[4],
                             kl_full=kl_full, agree=bool(ok)))
            print(f"    {e_t:7.2f} {crit:+12.4e} {kl_star:11.3e} "
                  f"{bb[4]:8.4f} {kl_full:11.3e} {'yes' if ok else 'NO':>6}",
                  flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "criterion_vs_projection.csv", index=False)

    n_ok = int(d.agree.sum())
    print(f"\n  the criterion and the full projection agree on {n_ok} of "
          f"{len(d)} grid points")
    if n_ok < len(d):
        print("\n  disagreements, where a non-star optimum beats the star")
        print("  without the local test at the star noticing:")
        print(d[~d.agree][["shape", "e_true", "criterion", "e_fit",
                           "star_kl", "kl_full"]].to_string(index=False))
    else:
        print("  so on this grid the boundary is exactly the criterion's")
        print("  zero, and no better point off the star was found anywhere.")


if __name__ == "__main__":
    main()
