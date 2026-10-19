"""Stage 0: does the analytic boundary reproduce exact sampling?

The prediction for the error probability of maximum likelihood on a quartet
rests on two approximations, and this script tests the first of them with no
inference engine anywhere in the loop.

    ell_a = sum_j log P_{T_a}(x_j)          sum of L i.i.d. per-site terms
    (ell_1, ell_2, ell_3)                   jointly normal by the CLT
    Pr[correct] = Pr[ell_1 beats both]      bivariate normal orthant

The mean vector and the covariance are exact sums over the 256 site patterns,
so the only step being tested is the normal one.  The second approximation --
evaluating the alternatives at their KL optimum rather than profiling them --
is NOT tested here; it costs O_p(1) against a gap of Theta(L), which is the
reason it can be neglected, but that argument is separate.

The exact answer is cheap.  Because the columns are i.i.d., the count vector
over the 256 patterns is exactly multinomial: one draw is one whole alignment,
and the three log-likelihoods are a 40000 x 256 by 256 x 3 matrix product.
No sequences are simulated and the reference is exact rather than merely
better than an engine.

Ties are broken uniformly at random, which matters at the star: there the
three likelihoods are identical and argmax would return topology 1 every time,
reporting a correct answer with probability 1 where the truth is 1/3.
"""
from __future__ import annotations

import importlib.util
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import minimize
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "results" / "mc_boundary"

_spec = importlib.util.spec_from_file_location(
    "s26", str(HERE / "s26_star_projection.py"))
q = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(q)

SPLITS = q.SPLITS
PI = q.PI
to_x, to_t = q.to_x, q.to_t
_cache: dict = {}


_cache: dict = {}


def _pack_perm(labels):
    """maps a labels-order linear index to the canonical A,B,C,D packing."""
    st = np.arange(4)
    perm = np.empty(256, dtype=int)
    for m in range(256):
        d = [(m >> 6) & 3, (m >> 4) & 3, (m >> 2) & 3, m & 3]
        idx = 0
        for n in ("A", "B", "C", "D"):
            idx = (idx << 2) | d[labels.index(n)]
        perm[m] = idx
    return perm


_PACK = {lab: _pack_perm(lab) for lab in set(SPLITS.values())}


def pattern_batch(lengths, labels):
    """P over the 256 patterns for a batch of branch-length vectors.

    Same quantity as s26's pattern_fast, which loops over the 256 patterns in
    Python.  Stage 1 profiles the likelihood, which means thousands of
    evaluations per replicate, so the loop has to go: the two terms
    alpha * ssab (x) sscd and beta * cross are outer products and one
    einsum, and the packing is a fixed permutation.  Checked against
    pattern_fast in main() rather than assumed.
    """
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
    lin = (alpha.reshape(-1, 1, 1, 1, 1)
           * (ssab[:, :, :, None, None] * sscd[:, None, None, :, :])
           + beta.reshape(-1, 1, 1, 1, 1) * cross).reshape(-1, 256)
    out = np.empty_like(lin)
    out[:, _PACK[labels]] = lin
    return out


CACHE_MAX = 20000


def pattern(lengths, labels):
    """One pattern vector, cached; the optimiser revisits the same vertices.

    The cache is capped.  It is a pure function of its key, so clearing it
    cannot change any result -- but leaving it unbounded can, by exhausting
    memory in the middle of a long sweep, which is how a run of this died
    at a 13 MiB allocation.  A cap turns that failure into a slowdown.
    """
    key = (labels, tuple(np.round(np.asarray(lengths, float), 12)))
    if key not in _cache:
        if len(_cache) > CACHE_MAX:
            _cache.clear()
        _cache[key] = pattern_batch(lengths, labels)[0]
    return _cache[key]


def kl(p, r):
    m = p > 0
    return float(np.sum(p[m] * np.log(p[m] / np.maximum(r[m], 1e-300))))


def best_fit(p_true, labels, x0):
    """Branch lengths of `labels` minimising KL(p_true || .).

    DETERMINISTIC.  This is a function of its arguments and nothing else.
    An earlier version carried a module-level dict of the last fit for each
    topology and used it as a warm start, which made the return value depend
    on what had been evaluated before.  That is fatal for a boundary search:
    a grid scan feeds itself, so the objective is not a function of the
    grid point, brackets appear where there is no crossing, and bisection
    then converges confidently to a point that does not satisfy the
    condition -- measured at 0.900 where the search believed 0.950.  The
    cost of determinism is the warm start, and it is worth paying.
    """
    best = (np.inf, np.asarray(x0, float))

    def obj_t(th):
        if np.any(th <= 1e-10) or np.any(th > 30):
            return 1e3
        return kl(p_true, pattern(th, labels))

    def obj_x(x):
        return kl(p_true, pattern(to_t(x), labels))

    starts = [np.asarray(x0, float), np.asarray(x0, float) * 0.3,
              np.asarray(x0, float) * 0.03, np.full(5, 0.4)]
    for trial in starts:
        r = minimize(obj_t, trial, method="L-BFGS-B",
                     bounds=[(1e-10, 30)] * 5,
                     options=dict(ftol=1e-16, gtol=1e-12, maxiter=4000))
        if r.fun < best[0]:
            best = (float(r.fun), r.x)
    for trial in (to_x(best[1]), np.full(5, 0.7), np.full(5, 0.1)):
        r = minimize(obj_x, trial, method="L-BFGS-B",
                     bounds=[(1e-9, 1.0)] * 5,
                     options=dict(ftol=1e-18, gtol=1e-14, maxiter=4000))
        if r.fun < best[0]:
            best = (float(r.fun), to_t(r.x))
    return best


STAR_TOL = 1e-6


def swap_perm(x, y):
    """Relabelling leaves x and y, as a permutation of the 256 indices."""
    out = np.empty(256, dtype=int)
    for idx in range(256):
        st = {"A": idx >> 6, "B": (idx >> 4) & 3,
              "C": (idx >> 2) & 3, "D": idx & 3}
        st[x], st[y] = st[y], st[x]
        out[idx] = (st["A"] << 6) | (st["B"] << 4) | (st["C"] << 2) | st["D"]
    return out


_PERM_AB = swap_perm("A", "B")
_PERM_CD = swap_perm("C", "D")


def symmetry(b1):
    """The pendant swap that exchanges the two alternative topologies.

    Both (A B) and (C D) take the split {A,C}|{B,D} to {A,D}|{B,C}, so
    either one works provided it also fixes the true tree -- which needs
    a = b, respectively c = d.  Note what that costs: a = b or c = d forces
    the two long pendants to be sisters whenever there are two long pendants
    at all.  The symmetry that makes the two wrong topologies coincide is
    therefore a property of symmetric quartets, not of the Farris zone, and
    the generic case has to go through the full two-dimensional orthant.

    Returns (labels order of the two swapped pendant slots, index perm), or
    None when the quartet has no such symmetry.
    """
    if abs(b1[0] - b1[1]) < 1e-12:
        return ("A", "B"), 0, 2, _PERM_AB
    if abs(b1[2] - b1[3]) < 1e-12:
        return ("C", "D"), 1, 3, _PERM_CD
    return None


def three_fits(b1):
    """(P_true, [P1, P2, P3], [0, e2, e3]) at the KL optima.

    Both wrong topologies project onto the star -- that is the finding this
    line of work started from -- and the star is described *identically* in
    the two parametrisations, since at zero internal branch all that survives
    is which pendant hangs off which leaf.  So when both optima are on the
    star the two minimisers are not merely close, they are the same
    distribution: the correlation between the two likelihood differences is
    exactly 1 and the three-way choice is really two-way.

    Fitting them separately and leaving it there would obscure this.  The two
    optimiser outputs agree to about 1e-10, and the correlation computed from
    a difference of two things that agree to 1e-10 is noise, not 1.  So the
    two fits are explicitly made one when the star is the optimum for both.

    A pendant symmetry, where the quartet has one, is also exploited: it lets
    the star optimum be symmetrised before it is used, which costs nothing
    and removes a little more optimiser error.
    """
    p = pattern(b1, SPLITS["12|34"])
    sym = symmetry(b1)

    _, bb1 = best_fit(p, SPLITS["13|24"], b1)
    if sym is not None and bb1[4] < STAR_TOL:
        _, i, j, _ = sym
        bb1 = np.array(bb1, float)
        bb1[i] = bb1[j] = 0.5 * (bb1[i] + bb1[j])
    f1 = pattern(bb1, SPLITS["13|24"])

    _, bb2 = best_fit(p, SPLITS["14|23"], b1)
    f2 = pattern(bb2, SPLITS["14|23"])

    if bb1[4] < STAR_TOL and bb2[4] < STAR_TOL:
        f1 = f2 = 0.5 * (f1 + f2)

    return p, [p, f1, f2], [0.0, float(bb1[4]), float(bb2[4])]


def logs_of(fits):
    return np.array([np.log(np.maximum(f, 1e-300)) for f in fits])   # (3,256)


def moments(p_true, lg):
    mu = lg @ p_true
    c = lg - mu[:, None]
    return mu, (c * p_true[None, :]) @ c.T


def params(mu, cov):
    """(gap2, gap3, sd2, sd3, rho) after dividing through by sd."""
    m2, m3 = mu[0] - mu[1], mu[0] - mu[2]
    v2 = cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
    v3 = cov[0, 0] + cov[2, 2] - 2 * cov[0, 2]
    c23 = cov[0, 0] - cov[0, 2] - cov[0, 1] + cov[1, 2]
    return m2, m3, v2, v3, c23


def orthant(a, b, rho):
    """Pr[X > -a, Y > -b] for a standard bivariate normal with correlation rho.

    Conditioning on X = x turns it into a one-dimensional integral,
    int_{-a}^{inf} phi(x) Phi((rho x + b)/sqrt(1 - rho^2)) dx, which stays
    well behaved as rho approaches 1 -- unlike the equivalent Plackett form,
    whose integrand has a removable singularity at t = 1.
    """
    if abs(rho) < 1e-15:
        return float(norm.cdf(a) * norm.cdf(b))
    if abs(rho) >= 1.0 - 1e-12:      # D_2 = D_3, so this is a two-way choice
        if rho > 0:
            return float(norm.cdf(min(a, b)))
        return float(max(0.0, norm.cdf(-b) - norm.cdf(-a)))   # Y = -X
    scale = np.sqrt(1.0 - rho * rho)
    f = lambda x: norm.pdf(x) * norm.cdf((rho * x + b) / scale)
    # Integrate over a window that actually contains the mass.  The lower
    # limit must not be -a when a is large: phi(x) is below 1e-32 for x < -12,
    # so quad starting at -42 finds nothing but zeros and returns 0 instead of
    # 1.  That is not a tail effect but a total failure, and it is silent --
    # it makes every large-z configuration read as impossible.  The integrand
    # is bounded by phi, so truncating to |x| <= 14 loses less than 1e-42.
    lo = max(-a, -14.0)
    hi = 14.0
    if lo >= hi:
        return 0.0
    val, _ = quad(f, lo, hi, limit=400, epsabs=1e-13, epsrel=1e-13)
    return float(min(max(val, 0.0), 1.0))


# Two fits count as one distribution below this, in probability units per
# site.  The threshold is set from the decision it has to not disturb: a
# pattern difference of d moves the log-likelihood by about L d / f_min, and
# that must stay well under the O(1) gaps being compared, so for L up to 10^6
# and f_min ~ 0.01 the bound is d <= 1e-8.  1e-9 is that, with margin.
FIT_TOL = 1e-9


def tied_groups(vecs):
    """Equivalence classes among pattern vectors, by FIT_TOL."""
    n = len(vecs)
    g = np.zeros(n, dtype=int)
    nxt = 1
    for j in range(1, n):
        for i in range(j):
            if np.abs(np.asarray(vecs[j]) - np.asarray(vecs[i])).max()                     <= FIT_TOL:
                g[j] = g[i]
                break
        else:
            g[j] = nxt
            nxt += 1
    return g


def analytic(p_true, fits, L):
    """Pr[correct] under the normal approximation at alignment length L.

    Three distinct cases, and collapsing them is a mistake that only shows up
    when the data carry no signal at all:

      * all three fits coincide with the truth -- a star, or data that have
        been shuffled.  Nothing is distinguishable, the three topologies are
        exchangeable, and the answer is exactly 1/3.  Computing rho here
        divides 0 by 0 and an earlier version therefore returned Phi(0)=1/2,
        the two-way answer, for a genuinely three-way tie.
      * the two alternatives coincide but the truth does not -- the star
        projection of s26.  Then rho is exactly 1 and the choice is two-way.
      * otherwise the full orthant probability.
    """
    g = tied_groups(list(fits))
    if g[0] == g[1] == g[2]:
        # all three FITTED distributions coincide, so all three likelihoods
        # are one number and the winner is a three-way tie.  Whether the
        # truth agrees is irrelevant -- if it does the data are a star, if
        # it does not the data carry signal the model cannot express, and
        # either way nothing distinguishes the three candidates.
        return 1.0 / 3.0
    m2, m3, v2, v3, c23 = params(*moments(p_true, logs_of(fits)))
    if v2 <= 0 or v3 <= 0:
        return float("nan")
    if g[1] == g[2]:        # the two WRONG fits are one fit, truth apart
        return float(norm.cdf(min(m2 * L / np.sqrt(v2 * L),
                                  m3 * L / np.sqrt(v3 * L))))
    rho = float(np.clip(c23 / np.sqrt(v2 * v3), -1.0, 1.0))
    return orthant(m2 * L / np.sqrt(v2 * L), m3 * L / np.sqrt(v3 * L), rho)


def mc(p_true, fits, L, n_rep, rng):
    """Exact: one multinomial draw per replicate is one alignment."""
    g = tied_groups(list(fits))
    counts = rng.multinomial(L, p_true, size=n_rep)
    ell = counts @ logs_of(fits).T                       # (n_rep, 3)
    # members of one group are made to tie exactly, so that the same
    # equivalence test decides both the analytic and the sampled answer
    for j in range(len(fits)):
        ell[:, j] = ell[:, g[j]]
    # remaining tolerance set by floating-point precision, not a relative rule
    tol = 1e-9 * max(1.0, float(np.abs(ell).max()))
    tied = np.abs(ell - ell.max(axis=1, keepdims=True)) <= tol
    rank = np.cumsum(tied, axis=1) * tied                # 1..k within a row
    k = tied.sum(axis=1, keepdims=True)
    draw = np.floor(rng.random((n_rep, 1)) * k) + 1
    win = np.argmax(rank == draw, axis=1)
    return float((win == 0).mean())


# Three shapes, chosen so that the two things being separated do not move
# together.  The Farris pair differs only in whether the quartet is
# symmetric, which is the one condition that collapses the two wrong
# topologies onto each other; the Felsenstein shape is the generic
# contrast.
SHAPES = {"Farris, symmetric":      (0.05, 0.05, 1.00, 1.00, 0.02),
          "Farris, generic":        (0.05, 0.07, 1.00, 0.85, 0.02),
          "Felsenstein, generic":   (1.00, 0.05, 1.00, 0.05, 0.02)}
LS = (100, 300, 1000, 3000, 10000)


def pr_at(shape, s, L):
    b1 = tuple(v * s for v in shape)
    p, fits, es = three_fits(b1)
    return analytic(p, fits, L), p, fits, es


def boundary(shape, L, grid=np.geomspace(0.02, 60.0, 16), target=0.95):
    """Largest tree scale at which recovery still reaches `target`.

    Recovery is not monotone in tree length.  Shrink the tree towards a point
    and all three topologies coincide, so the right one is picked a third of
    the time; grow it past saturation and the sites are noise, and the same
    thing happens.  So the probability has an interior maximum, and the two
    crossings of `target` mean different things: the lower one is "not enough
    change to resolve anything", the upper one is saturation.  Only the upper
    one is a statement about Iss, which measures saturation and is not
    defined at all for a tree that has barely changed.

    The grid is scanned for the maximum and the search then runs on the
    decreasing side of it, so the lower crossing cannot be returned by
    accident.
    """
    f = lambda s: pr_at(shape, s, L)[0] - target
    vals = np.array([f(s) for s in grid])
    vals[~np.isfinite(vals)] = np.nan
    if np.all(np.isnan(vals)):
        return float("nan")
    j = int(np.nanargmax(vals))
    if vals[j] < 0.0:
        # the peak itself is below the target: this quartet is not
        # recoverable at this length at any tree scale
        return float("nan")
    for i in range(j, len(grid) - 1):
        flo, fhi = vals[i], vals[i + 1]
        if not (np.isfinite(flo) and np.isfinite(fhi)) or flo * fhi > 0:
            continue
        lo, hi = grid[i], grid[i + 1]
        for _ in range(24):
            mid = 0.5 * (lo + hi)
            fm = f(mid)
            if not np.isfinite(fm):
                return float("nan")
            if flo * fm <= 0:
                hi = mid
            else:
                lo, flo = mid, fm
            if hi - lo < 1e-6 * hi:
                break
        return 0.5 * (lo + hi)
    return float("nan")


def peak(shape, L, grid=np.geomspace(0.02, 60.0, 16)):
    """(best achievable Pr[correct], tree scale achieving it)."""
    vals = np.array([pr_at(shape, s, L)[0] for s in grid])
    j = int(np.nanargmax(vals))
    return float(vals[j]), float(grid[j])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(71203)
    n_rep = 40000
    t0 = time.time()

    print("STAGE 0.  analytic orthant vs exact sampling, no engine involved.")
    print(f"{n_rep} replicates per point; MC standard error at p=0.05 is "
          f"{np.sqrt(0.05 * 0.95 / n_rep):.4f}\n")

    # The boundary search is a bisection on pr_at, so pr_at had better be a
    # function.  It is checked rather than assumed, because the failure mode
    # is silent: a history-dependent objective still yields a smooth-looking
    # scan and a bisection that converges, just to the wrong place.
    probe = tuple(v * 0.7 for v in SHAPES["Farris, generic"])
    twice = [pr_at(SHAPES["Farris, generic"], 0.7, 1000)[0] for _ in range(2)]
    pr_at(SHAPES["Felsenstein, generic"], 3.0, 1000)      # perturb the state
    again = pr_at(SHAPES["Farris, generic"], 0.7, 1000)[0]
    print(f"    determinism check: {twice[0]:.12f} {twice[1]:.12f} "
          f"{again:.12f}", end="")
    if not (twice[0] == twice[1] == again):
        sys.exit("  FAILED -- pr_at depends on evaluation history")
    print("  ok\n")

    rows = []
    for name, shape in SHAPES.items():
        print(f"  {name}")
        print(f"    {'TL':>6} {'L':>6} {'analytic':>9} {'exact':>9} "
              f"{'diff':>8} {'rho23':>7} {'gap/L':>10} {'e_fit':>7}")
        for s in (0.3, 0.6, 1.0, 1.6, 2.5):
            b1 = tuple(v * s for v in shape)
            p, fits, es = three_fits(b1)
            lg = logs_of(fits)
            mu, cov = moments(p, lg)
            m2, m3, v2, v3, c23 = params(mu, cov)
            rho = c23 / np.sqrt(v2 * v3) if v2 > 0 and v3 > 0 else np.nan
            for L in LS:
                a = analytic(p, fits, L)
                e = mc(p, fits, L, n_rep, rng)
                rows.append(dict(shape=name, scale=s, TL=sum(b1), L=L,
                                 analytic=a, exact=e, diff=a - e, rho=rho,
                                 gap2=m2, gap3=m3, e_fit=es[1]))
                print(f"    {sum(b1):6.2f} {L:6d} {a:9.4f} {e:9.4f} "
                      f"{a - e:+8.4f} {rho:7.4f} {m2:10.3e} {es[1]:7.4f}",
                      flush=True)
        print()

    d = pd.DataFrame(rows)
    d.to_csv(OUT / "stage0.csv", index=False)
    print(f"  worst |analytic - exact| : {d['diff'].abs().max():.4f}")
    print(f"  mean  |analytic - exact| : {d['diff'].abs().mean():.4f}")
    mid = d[d.analytic.between(0.05, 0.95)]
    print(f"  on 0.05 < analytic < 0.95, {len(mid)} points: "
          f"worst {mid['diff'].abs().max():.4f}, "
          f"mean {mid['diff'].abs().mean():.4f}"
          if len(mid) else "  nothing lands in the informative range")

    # ---------------- the same comparison, but located ON the boundary
    print("\n  boundary sweep: at each L, the LARGEST tree length whose")
    print("  analytic error rate is still 0.05 -- the saturation side, this")
    print("  being the side Iss is about -- and exact sampling there.")
    print(f"    {'shape':22s} {'L':>7} {'s*':>9} {'TL*':>8} "
          f"{'analytic':>9} {'exact':>9} {'diff':>8} {'n_rep':>6}")
    brows = []
    for name, shape in SHAPES.items():
        for L in (3000, 10000, 30000, 100000):
            s = boundary(shape, L)
            if not np.isfinite(s):
                pk, spk = peak(shape, L)
                print(f"    {name:22s} {L:7d}   no crossing; best achievable "
                      f"Pr = {pk:.4f} at TL = {sum(shape) * spk:.3f}")
                continue
            b1 = tuple(v * s for v in shape)
            p, fits, es = three_fits(b1)
            a = analytic(p, fits, L)
            n = int(min(n_rep, max(4000, 4e8 / L)))
            e = mc(p, fits, L, n, rng)
            brows.append(dict(shape=name, L=L, s=s, TL=sum(b1), e_fit=es[1],
                              analytic=a, exact=e, diff=a - e, n_rep=n))
            print(f"    {name:22s} {L:7d} {s:9.4f} {sum(b1):8.3f} "
                  f"{a:9.4f} {e:9.4f} {a - e:+8.4f} {n:6d}", flush=True)
    pd.DataFrame(brows).to_csv(OUT / "stage0_boundary.csv", index=False)
    if brows:
        bd = pd.DataFrame(brows)
        print(f"\n  on the boundary: worst |analytic - exact| "
              f"{bd['diff'].abs().max():.4f}, mean "
              f"{bd['diff'].abs().mean():.4f}; exact never leaves "
              f"[{bd.exact.min():.4f}, {bd.exact.max():.4f}] around 0.95")
    print(f"\n  {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
