"""Stage 1: does profiling the branch lengths change the answer?

Stage 0 evaluated every topology's likelihood at a fixed parameter vector --
the truth's for the true topology, the KL optima for the two alternatives --
and got the right error rate.  That is an approximation, and it is the one
Stage 0 explicitly did not test.  Real inference maximises over the branch
lengths, so the quantity a program actually compares is

    max_theta log P_{T_a}(x ; theta),

not log P_{T_a}(x ; theta_a^*).  Each maximum exceeds its fixed-parameter
value by O_p(1), and the gap between the topologies is Theta(L) -- but the
boundary is not where that gap is large.  The boundary is exactly where the
gap has shrunk to the order of the sampling noise, Theta(sqrt(L)), so an
O_p(1) perturbation does NOT become negligible there and the argument that
it must is wrong.  This script measures how far the boundary actually moves.

Note also what the O_p(1) term is.  It is not a shift.  Profiling gains
about 2.5 nats for the true topology and about 2.3 for each alternative --
all three are maximised, so all three gain -- and the *difference* of two
such gains has a standard deviation of order 1 that adds to the sampling
variance rather than cancelling it.  At L = 10^4 on a configuration whose
fixed-parameter sd is 0.66, the profiled decision statistic has sd 1.27.

The profiler is checked against a tighter one in s28's own history, and
against the fixed-parameter likelihoods: the fitted branch lengths converge
to the population optima, the gains land near chi-square/2 with 5 degrees of
freedom, and re-running at ftol=1e-16 with four restarts reproduces the
numbers exactly.  That matters because a sloppy optimiser would manufacture
precisely this variance inflation out of nothing.

Warm starts matter.  Each replicate starts from the population optimum for
that topology, which is where the sample optimum is for large L.
"""
from __future__ import annotations

import importlib.util
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "results" / "mc_boundary"

_spec = importlib.util.spec_from_file_location(
    "s27", str(HERE / "s27_mc_boundary.py"))
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)

SPLITS = m.SPLITS
LABELS = list(SPLITS)
BOUNDS = [(1e-6, 30.0)] * 5
FLOOR = 1e-300


def profile(counts, labels, x0):
    """max_theta sum_x n_x log P_{T}(x; theta), and the argmax.

    Returns the LOG-likelihood, not its negative.  The optimiser minimises
    the negative one and an earlier version of this function passed r.fun
    straight through, which made every comparison below a comparison of
    negated likelihoods -- so the winner was the worst topology, and the
    script reported that profiling destroys the inference.  A log-likelihood
    is at most zero here, so a positive value anywhere is the signature.
    """
    logs = {}

    def negll(th):
        key = tuple(np.round(th, 12))
        if key not in logs:
            logs[key] = np.log(np.maximum(m.pattern_batch(th, labels)[0],
                                          FLOOR))
        return -float(counts @ logs[key])

    best = (np.inf, np.asarray(x0, float))
    for trial in (np.asarray(x0, float),
                  np.maximum(np.asarray(x0, float) * 0.5, 1e-5)):
        r = minimize(negll, trial, method="L-BFGS-B", bounds=BOUNDS,
                     options=dict(maxiter=400, maxfun=4000))
        if r.fun < best[0]:
            best = (float(r.fun), np.asarray(r.x, float))
    ll = -best[0]
    if ll > 1e-8:
        raise AssertionError(f"log-likelihood came out positive ({ll}); "
                             f"the sign convention has flipped again")
    return ll, best[1]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(90210)
    n_rep = 300
    LS = (3000, 10000, 30000)
    t0 = time.time()

    print("STAGE 1.  profiled maximum likelihood AT the stage-0 boundary.")
    print(f"{n_rep} replicates per point, three topologies, five branch "
          f"lengths each.\n")
    print("  Stage 0 compares topology 1 at its TRUE parameters against the")
    print("  two alternatives at their KL optima -- an oracle, since it is")
    print("  handed the answer.  Real inference profiles all three, and the")
    print("  O_p(1) it gains is not only a shift: it carries O_p(1) variance")
    print("  of its own, which does not vanish relative to sqrt(v2 L) at the")
    print("  boundary.  What matters is how far the boundary moves.\n")
    print(f"  {'shape':22s} {'L':>7} {'TL*':>7} {'stage0':>8} "
          f"{'profiled':>9} {'shift':>8} {'gap sd':>7} {'fixed sd':>9} "
          f"{'se':>6}")

    rows = []
    for name, shape in m.SHAPES.items():
        for L in LS:
            s = m.boundary(shape, L)
            if not np.isfinite(s):
                print(f"  {name:22s} {L:7d}       no crossing at any tree "
                      f"length")
                continue
            b1 = tuple(v * s for v in shape)
            p_true, fits, es = m.three_fits(b1)
            starts = [np.asarray(b1, float)]
            for lab in LABELS[1:]:
                _, bb = m.best_fit(p_true, SPLITS[lab], b1)
                starts.append(np.asarray(bb, float))

            counts = rng.multinomial(L, p_true, size=n_rep)
            wins = np.zeros(n_rep, dtype=int)
            gaps = np.zeros(n_rep)
            fgaps = np.zeros(n_rep)
            lgf = np.log(np.maximum(np.array(fits), 1e-300))
            for r_i in range(n_rep):
                vals = [profile(counts[r_i], SPLITS[lab], starts[a_i])[0]
                        for a_i, lab in enumerate(LABELS)]
                gaps[r_i] = vals[0] - vals[1]
                fg = counts[r_i] @ lgf.T
                fgaps[r_i] = fg[0] - fg[1]
                top = max(vals)
                tied = [i for i, v in enumerate(vals) if top - v <= 1e-6]
                wins[r_i] = int(rng.choice(tied))
            pc = float((wins == 0).mean())
            se = float(np.sqrt(max(pc * (1 - pc), 1e-12) / n_rep))
            rows.append(dict(shape=name, L=L, scale=s, TL=sum(b1),
                             stage0=0.95, profiled=pc, shift=pc - 0.95,
                             gap_sd=float(gaps.std()),
                             fixed_sd=float(fgaps.std()), se=se,
                             n_rep=n_rep))
            print(f"  {name:22s} {L:7d} {sum(b1):7.3f} {0.95:8.4f} "
                  f"{pc:9.4f} {pc - 0.95:+8.4f} {gaps.std():7.3f} "
                  f"{fgaps.std():9.3f} {se:6.4f}", flush=True)
            pd.DataFrame(rows).to_csv(OUT / "stage1.csv", index=False)
    print(f"\n  {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
