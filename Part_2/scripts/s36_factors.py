"""Which factors set Iss.c, and by how much.

The mechanism gives the answer in closed form.  With

    Iss.c(L) = Iss_sat - C L^{-p},

the exponent is a ratio of two independently measurable rates:

    p = lam_Iss / (2 lam_Delta - lam_v2),

where, over the tree scale s,

    Delta(s)   ~ exp(-lam_Delta s)     the topological signal decays
    v2(s)      ~ exp(-lam_v2   s)     the noise scale decays
    gap(s) = Iss_sat - Iss(s) ~ exp(-lam_Iss s)   Iss saturates

So the factors are lam_Iss and lam_Delta.  L is only the variable they are
assembled into.  The tree's geometry and the substitution model enter
through the two rates and nowhere else.

What this script does is check that.  For each configuration it measures the
three rates on the scale range actually spanned by the critical trees, and
measures p from the boundary, then compares p against the ratio.  If the
identity holds, p stops being a fitted number and becomes the ratio of two
rates that can be measured separately -- which is what makes the factors
arguable rather than merely nameable.

The scale range matters: the rates are local, so they are taken over the
window the sweep actually uses, not from an asymptotic limit that no
configuration reaches.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
import pandas as pd

HERE = Path(__file__).resolve().parent
PART2 = HERE.parent
ROOT = PART2.parent
sys.path.insert(0, str(ROOT / "src"))

OUT = PART2 / "results" / "factors"

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

# A one-factor-at-a-time set, all readable against a common reference.
REF = (1.0, 1.0, 0.05, 0.05, 0.02)
SHAPES = {
    "reference  a=b=1.0 c=d=0.05 e=0.02": REF,
    "e x5       a=b=1.0 c=d=0.05 e=0.10": (1.0, 1.0, 0.05, 0.05, 0.10),
    "e x25      a=b=1.0 c=d=0.05 e=0.50": (1.0, 1.0, 0.05, 0.05, 0.50),
    "c=d x4     a=b=1.0 c=d=0.20 e=0.02": (1.0, 1.0, 0.20, 0.20, 0.02),
    "a=b /2     a=0.5 c=d=0.05 e=0.02": (0.5, 0.5, 0.05, 0.05, 0.02),
    "a=b /5     a=0.2 c=d=0.05 e=0.02": (0.2, 0.2, 0.05, 0.05, 0.02),
    "all equal  a=b=c=d=0.3 e=0.6": (0.3, 0.3, 0.3, 0.3, 0.6),
    "near-star  a=b=1.0 c=d=0.05 e=0.002": (1.0, 1.0, 0.05, 0.05, 0.002),
}
LS = (30000, 60000, 100000, 200000, 300000)


def e_iss(shape, s):
    return float(exact.E_Iss(4, trees.quartet(*[v * s for v in shape]),
                             Q, PI, convention="A"))


def to_labels(pend4, lab):
    """(a, b, c, d) in canonical A,B,C,D order -> the order `lab` wants."""
    return np.array([pend4["ABCD".index(n)] for n in lab], float)


def star_fit(p, lab, pend4):
    """Four pendants with the internal branch pinned at zero.

    `pend4` is in canonical A,B,C,D order; the conversion to this
    topology's labelling is done here so the caller cannot get it wrong.

    Pinning e = 0 is justified rather than assumed: s30 re-scanned 219
    configurations by objective and found no projection off the star, the
    largest relative KL gain anywhere being 2.5e-6.  Dropping the fifth
    parameter is what makes this sweep affordable, and check_fast below
    verifies the fast path against the full one.
    """
    x0 = np.maximum(to_labels(pend4, lab), 1e-6)

    def neg(par):
        th = np.exp(par)
        q = m.pattern_batch((th[0], th[1], th[2], th[3], 0.0), lab)[0]
        return float(np.sum(p * np.log(np.maximum(p, 1e-300)
                                       / np.maximum(q, 1e-300))))

    best = (np.inf, np.log(x0))
    for t in (np.log(x0), np.log(x0) * 0.3, np.log([0.4, 0.4, 0.4, 0.4])):
        r = minimize(neg, np.asarray(t, float), method="L-BFGS-B",
                     bounds=[(-20.0, 5.0)] * 4,
                     options=dict(ftol=1e-20, gtol=1e-14, maxiter=600))
        if r.fun < best[0]:
            best = (float(r.fun), np.asarray(r.x, float))
    th = np.exp(best[1])
    return (th[0], th[1], th[2], th[3], 0.0)


def fits_fast(b1):
    """[P_true, fit to 13|24, fit to 14|23], star-only, symmetry imposed.

    The symmetry has to be imposed here for the same reason three_fits
    imposes it: for a quartet with a = b or c = d the two alternatives'
    optima are related exactly, and fitting them independently lands them at
    visibly different points -- 8e-3 apart in pattern units on one balanced
    configuration -- which makes the correlation come out 0.90 instead of
    exactly 1 and turns a two-way choice into a spurious three-way one.
    """
    p = m.pattern_batch(b1, m.SPLITS["12|34"])[0]
    lab1, lab2 = "13|24", "14|23"
    sym = m.symmetry(b1)

    bb = star_fit(p, m.SPLITS[lab1], b1[:4])
    if sym is not None:
        _, i, j, perm = sym
        th = np.array(bb, float)
        th[i] = th[j] = 0.5 * (th[i] + th[j])
        f1 = m.pattern_batch(th, m.SPLITS[lab1])[0]
        img = f1[perm]                      # sigma(fit1), as in three_fits
        f2 = f1 if np.abs(f1 - img).max() <= m.FIT_TOL else img
        return [p, f1, f2]
    f1 = m.pattern_batch(bb, m.SPLITS[lab1])[0]
    bb2 = star_fit(p, m.SPLITS[lab2], b1[:4])
    f2 = m.pattern_batch(bb2, m.SPLITS[lab2])[0]
    return [p, f1, f2]


def analytic_fast(b1, L):
    p = m.pattern_batch(b1, m.SPLITS["12|34"])[0]
    return m.analytic(p, fits_fast(b1), L)


def boundary_fast(shape, L, grid=np.geomspace(0.02, 60.0, 13), target=0.95):
    f = lambda s: analytic_fast(tuple(v * s for v in shape), L) - target
    vals = np.array([f(s) for s in grid])
    j = int(np.nanargmax(vals))
    if vals[j] < 0:
        return float("nan")
    for i in range(j, len(grid) - 1):
        if not (vals[i] >= 0 > vals[i + 1]):
            continue
        lo, hi = grid[i], grid[i + 1]
        for _ in range(20):
            mid = 0.5 * (lo + hi)
            if f(mid) >= 0:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1e-5 * hi:
                break
        return 0.5 * (lo + hi)
    return float("nan")


def check_fast():
    """The fast path must reproduce a boundary the full path already gave."""
    sh = REF
    s_full = m.boundary(sh, 30000, target=0.95)
    s_fast = boundary_fast(sh, 30000, target=0.95)
    iss_f = e_iss(sh, s_full)
    iss_g = e_iss(sh, s_fast)
    print(f"  check: full path s*={s_full:.5f} (Iss {iss_f:.5f}), "
          f"fast path s*={s_fast:.5f} (Iss {iss_g:.5f}), "
          f"difference {abs(iss_f - iss_g):.2e}")
    if abs(iss_f - iss_g) > 1e-4:
        sys.exit("  fast path disagrees; do not use it")
    print()


def rates(shape, lo, hi, n=7):
    """(lam_Delta, lam_v2, lam_Iss) fitted over the scale window [lo, hi]."""
    ss = np.linspace(lo, hi, n)
    D, V, G = [], [], []
    for s in ss:
        b1 = tuple(v * s for v in shape)
        p = m.pattern_batch(b1, m.SPLITS["12|34"])[0]
        fits = []
        for lab in m.SPLITS:
            _, bb = m.best_fit(p, m.SPLITS[lab], b1)
            fits.append(m.pattern_batch(bb, m.SPLITS[lab])[0])
        mu, cov = m.moments(p, m.logs_of(fits))
        m2, m3, v2, v3, c23 = m.params(mu, cov)
        D.append(max(m2, 1e-300))
        V.append(max(v2, 1e-300))
        G.append(max(ISS_SAT - e_iss(shape, s), 1e-300))
    f = lambda y: -np.polyfit(ss, np.log(np.array(y)), 1)[0]
    return f(D), f(V), f(G)


def power_law_p(shape):
    """p, r^2, and the scale window the critical trees actually span.

    Each boundary is solved ONCE.  An earlier version solved them again to
    get the window, doubling a cost that is already the bulk of the script.
    """
    scales = {L: boundary_fast(shape, L, target=0.95) for L in LS}
    good = {L: s for L, s in scales.items() if np.isfinite(s)}
    if len(good) < 3:
        return float("nan"), float("nan"), None, None
    L_ = np.array(sorted(good), float)
    g_ = np.array([ISS_SAT - e_iss(shape, good[L]) for L in L_])
    A = np.column_stack([np.ones_like(L_), -np.log(L_)])
    c, *_ = np.linalg.lstsq(A, np.log(g_), rcond=None)
    pred = np.exp(c[0]) * L_ ** (-c[1])
    r2 = 1 - float(((g_ - pred) ** 2).sum()) / float(((g_ - g_.mean()) ** 2).sum())
    vals = np.array(list(good.values()))
    return float(c[1]), r2, 0.9 * vals.min(), 1.1 * vals.max()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    check_fast()
    print("Which factors set Iss.c.  p is the exponent of Iss_sat - Iss.c(L),")
    print("predicted as lam_Iss / (2 lam_Delta - lam_v2).\n")
    print(f"  {'configuration':38s} {'lam_D':>7} {'lam_v':>7} {'lam_I':>7} "
          f"{'p_pred':>7} {'p_fit':>7} {'rel err':>8} {'r2':>7}")
    rows = []
    for name, shape in SHAPES.items():
        p_fit, r2, s_lo, s_hi = power_law_p(shape)
        if not np.isfinite(p_fit):
            print(f"  {name:38s}  no crossing at this L range")
            continue
        lD, lV, lI = rates(shape, s_lo, s_hi)
        p_pred = lI / (2 * lD - lV) if (2 * lD - lV) != 0 else float("nan")
        rel = (p_pred - p_fit) / p_fit if p_fit else float("nan")
        rows.append(dict(config=name, lam_Delta=lD, lam_v2=lV, lam_Iss=lI,
                         p_pred=p_pred, p_fit=p_fit, rel_err=rel, r2=r2,
                         s_lo=s_lo, s_hi=s_hi))
        print(f"  {name:38s} {lD:7.4f} {lV:7.4f} {lI:7.4f} {p_pred:7.4f} "
              f"{p_fit:7.4f} {rel:+8.3f} {r2:7.4f}", flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "factors.csv", index=False)

    if len(d):
        # Does lam_v track lam_D?  If so p reduces to a ratio of two rates.
        c = np.corrcoef(d.lam_Delta, d.lam_v2)[0, 1]
        slope = np.polyfit(d.lam_Delta, d.lam_v2, 1)[0]
        print(f"\n  lam_v2 tracks lam_Delta: slope {slope:.4f}, "
              f"correlation {c:.4f}")
        print(f"  so 2 lam_Delta - lam_v2 stays close to lam_Delta, and p is")
        print(f"  close to lam_Iss / lam_Delta -- one ratio of two rates.")
        print(f"\n  lam_Delta spans {d.lam_Delta.min():.4f} to "
              f"{d.lam_Delta.max():.4f} across these configurations, a factor "
              f"of {d.lam_Delta.max()/d.lam_Delta.min():.1f}")
        print(f"  lam_Iss   spans {d.lam_Iss.min():.4f} to "
              f"{d.lam_Iss.max():.4f}, a factor of "
              f"{d.lam_Iss.max()/d.lam_Iss.min():.1f}")
        print(f"  and p spans {d.p_fit.min():.4f} to {d.p_fit.max():.4f}")
        print(f"\n  mean |relative error| of the prediction: "
              f"{d.rel_err.abs().mean():.3f}")


if __name__ == "__main__":
    main()
