"""The shuffle as a null model for the decision statistic.

Under the fully shuffled distribution the population fits of all three
topologies coincide -- they are the same star -- so in the POPULATION the
statistic l1 - l2 is identically zero and its null variance is zero.  That is
degeneracy, not a result, and it is why a fully shuffled null looks useless
at first glance.

In a SAMPLE it is not degenerate.  The three MLEs all converge to the same
population point but land in different places at finite L, so l1 - l2 has mean
zero and a nonzero variance whose only source is estimation.  That variance is
the noise floor of the decision, and it is exactly the quantity Stage 1
measured but could not predict: profiling raised the sd of the decision
statistic from 0.66 to 1.27 on one configuration, an excess of sigma_p = 1.08
attributed to the profile with no independent handle on it.

So this script asks a sharp question.  Does the shuffle null reproduce that
number?

If it does, three things are gained at once:

  * a distribution-free calibration of the decision statistic, which does not
    rest on the normal approximation;
  * an independent measurement of the profile noise, closing the open item
    from Stage 1;
  * a way to express any signal as a multiple of the noise floor that is
    measured on the data rather than assumed.

If it does not, the profile excess comes from somewhere else and the
attribution in Stage 1 was wrong, which is worth knowing too.
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

OUT = PART2 / "results" / "shuffle_null"

_spec27 = importlib.util.spec_from_file_location(
    "s27", str(HERE / "s27_mc_boundary.py"))
m = importlib.util.module_from_spec(_spec27)
_spec27.loader.exec_module(m)

_spec34 = importlib.util.spec_from_file_location(
    "s34", str(HERE / "s34_iss_null_path.py"))
s34 = importlib.util.module_from_spec(_spec34)
_spec34.loader.exec_module(s34)

_spec28 = importlib.util.spec_from_file_location(
    "s28", str(HERE / "s28_stage1_profile.py"))
s28 = importlib.util.module_from_spec(_spec28)
_spec28.loader.exec_module(s28)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261015)

    b1 = (0.05, 0.05, 1.0, 1.0, 0.02)
    L = 30000
    n = 400

    p_true = m.pattern_batch(b1, m.SPLITS["12|34"])[0]
    q = s34.permuted_average(p_true)

    print("Is the shuffle null the profile noise floor?")
    print(f"  tree a=b=0.05, c=d=1.0, e=0.02 (TL={sum(b1):.2f}), L={L}, "
          f"{n} replicates\n")

    rows = []
    for label, dist in (("true data  P", p_true),
                        ("shuffled   Q", q),
                        ("half-half  P/2+Q/2", 0.5 * p_true + 0.5 * q)):
        counts = rng.multinomial(L, dist, size=n)
        wins = np.zeros(n, dtype=int)
        gaps = np.zeros((n, 2))

        # The starting points matter more than they should.  Started only
        # from the ORIGINAL tree's parameters, the alternatives stick at a
        # tree-like optimum that assigns near-zero probability to patterns a
        # shuffled alignment actually contains, and the log-likelihood
        # collapses by ~14000 nats -- which reads as a massive signal against
        # topology 1 rather than as an optimiser failure.  A saturated start
        # is the honest one for a distribution with no tree structure.
        starts = []
        for a, lab in enumerate(m.SPLITS):
            cand = [np.asarray(b1, float), np.full(5, 3.0), np.full(5, 1.0),
                    np.full(5, 0.4)]
            if a > 0:
                _, bb = m.best_fit(dist, m.SPLITS[lab], b1)
                cand.insert(1, np.asarray(bb, float))
            starts.append(cand)

        for i in range(n):
            vals = [max(s28.profile(counts[i], m.SPLITS[lab], st)[0]
                        for st in starts[a])
                    for a, lab in enumerate(m.SPLITS)]
            gaps[i] = (vals[0] - vals[1], vals[0] - vals[2])
            top = max(vals)
            tied = [j for j, v in enumerate(vals) if top - v <= 1e-6]
            wins[i] = int(rng.choice(tied))

        d = gaps[:, 0]
        rows.append(dict(dist=label, mean=float(d.mean()),
                         sd=float(d.std(ddof=1)),
                         recovery=float((wins == 0).mean()), n=n))
        print(f"  {label}")
        print(f"      l1-l2 : mean {d.mean():+9.4f}   sd {d.std(ddof=1):8.4f}")
        print(f"      recovery of topology 1   {np.mean(wins == 0):7.4f}   "
              f"(1/3 = {1/3:.4f})", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "shuffle_null.csv", index=False)

    print("\n  Stage 1 measured, on the true data, that profiling raised the")
    print("  sd of l1-l2 from 0.66 to 1.27, an excess of about 1.08 with no")
    print("  independent handle on it.  The shuffled row above is that")
    print("  handle: its sd is the noise floor with the signal removed.")
    print("\n  Caveat: at full shuffling the POPULATION fits coincide, so the")
    print("  null variance here is pure estimation.  That is what makes it a")
    print("  calibration of the estimator rather than of the data.")


if __name__ == "__main__":
    main()
