"""What the three-nucleotide weight should be, from our own theory.

CASTER-site scores a quartet topology by summing a fixed vote over sites: +1
to the topology a two-base site resolves, and a weight w3 to the topology
that groups the pair sharing a base at a three-base site, with one-base and
four-base sites abstaining.  Because the vote depends only on the site, and
sites are i.i.d., the three topology scores are sums of i.i.d. three-vectors.
They are therefore *exactly* trivariate normal at every L, not merely
asymptotically, and

    Pr[correct] = Pr[S1 - S2 > 0 and S1 - S3 > 0]

is an exact bivariate normal orthant probability with no simulation and no
fitted parameter.  That is the same orthant machinery section 3 uses for
maximum likelihood, applied to a linear statistic instead of a likelihood.

This makes the CASTER weight an object we can optimise rather than accept.
The paper fixes w3 = -1/2 and says why: a negative vote repairs
Felsenstein-zone trees but "chips away at the score gap for Farris trees".
Under our own theory that trade-off is quantitative and has a location:

  * the value of w3 that maximises Pr[correct] depends on the tree, so no
    single constant is optimal for both zones;
  * below some critical w3 the Farris-zone margin E[S1 - S2] changes sign and
    the estimator stops being consistent, which is a hard boundary rather
    than a gradual loss.

Both are computed here, in closed form, alongside the maximum-likelihood
power on the same configurations so the two are read off one scale.

For maximum likelihood the per-site contribution is log P_i(pattern) with P_i
the projection of the true distribution onto topology i -- the fitted
distribution each topology converges to.  That is the oracle version of the
statistic, so the ML column is an upper bound on ML's power, which is the
conservative direction for the comparison.

Outputs
    results/caster_weights/weights.csv
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

OUT = PART2 / "results" / "caster_weights"

_spec = importlib.util.spec_from_file_location("s36", str(HERE / "s36_factors.py"))
s36 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s36)
m = s36.m

LABELS = ("12|34", "13|24", "14|23")
# a cherry pair -> the topology it indicates.  The complement pair of a
# two-base site indicates the same topology, which is why (2,3) and (0,1)
# agree.
PAIR_TOPO = {(0, 1): 0, (2, 3): 0,
             (0, 2): 1, (1, 3): 1,
             (0, 3): 2, (1, 2): 2}


def digits(idx):
    return ((idx >> 6) & 3, (idx >> 4) & 3, (idx >> 2) & 3, idx & 3)


def vote_table(w3):
    """(256, 3) per-site votes.  Scheme: a two-base site votes +1 for the one
    topology it resolves, a three-base site votes w3 for the topology that
    groups the pair sharing a base, and one- and four-base sites abstain.

    A two-base site with a 3+1 division resolves nothing and abstains: all
    three topologies place some pair of the majority trio together, so the
    site is symmetric between them.  That is the u = 0 convention.
    """
    V = np.zeros((256, 3))
    for idx in range(256):
        d = digits(idx)
        groups = {}
        for t, b in enumerate(d):
            groups.setdefault(b, []).append(t)
        if len(groups) == 2:
            parts = sorted(groups.values(), key=len)
            if len(parts[0]) == 2:
                V[idx, PAIR_TOPO[tuple(parts[0])]] = 1.0
        elif len(groups) == 3:
            pair = [v for v in groups.values() if len(v) == 2][0]
            V[idx, PAIR_TOPO[tuple(sorted(pair))]] = w3
    return V


def power(V, p, L):
    """Exact Pr[correct] for a linear site score, by orthant probability."""
    mu = V.T @ p                                   # (3,)
    second = V.T @ (V * p[:, None])                # E[v_i v_j] = sum_k p_k v_ki v_kj
    cov = second - np.outer(mu, mu)
    d12 = np.array([1.0, -1.0, 0.0])               # S1 - S2
    d13 = np.array([1.0, 0.0, -1.0])               # S1 - S3
    mx, my = L * (d12 @ mu), L * (d13 @ mu)
    vx, vy = L * (d12 @ cov @ d12), L * (d13 @ cov @ d13)
    cxy = L * (d12 @ cov @ d13)
    if vx <= 0 or vy <= 0:
        return float("nan")
    sx, sy = np.sqrt(vx), np.sqrt(vy)
    rho = float(np.clip(cxy / (sx * sy), -1.0, 1.0))
    return m.orthant(mx / sx, my / sy, rho)


def ml_table(p, init):
    """(256, 3) log-likelihood contributions from the three projections."""
    projs = []
    for lab in LABELS:
        _, bb = m.best_fit(p, m.SPLITS[lab], init)
        projs.append(m.pattern_batch(bb, m.SPLITS[lab])[0])
    P = np.stack(projs, axis=1)
    return np.log(np.maximum(P, 1e-300))


FARRIS = (1.0, 1.0, 0.05, 0.05, 0.02)       # long branches are sisters
FELSEN = (1.0, 0.05, 1.0, 0.05, 0.02)       # long branches are not


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    W3 = np.round(np.arange(-1.5, 0.65, 0.05), 4)
    print("The three-nucleotide weight, by our own theory.")
    print("  Pr[correct] is exact: the score is a sum of i.i.d. site votes,")
    print("  so it is exactly normal and the answer is an orthant probability.\n")

    rows = []
    for zone, shape in (("Farris", FARRIS), ("Felsenstein", FELSEN)):
        p0 = m.pattern_batch(shape, m.SPLITS["12|34"])[0]
        Vml = ml_table(p0, shape)
        for L in (10_000, 100_000, 1_000_000):
            ml = power(Vml, p0, L)
            curve = [(w3, power(vote_table(w3), p0, L)) for w3 in W3]
            ok = [c for c in curve if np.isfinite(c[1])]
            best_w3, best_pr = max(ok, key=lambda c: c[1])

            # The vote table is linear in w3, so the margin E[S1 - S2] is too:
            # margin(w3) = m0 + w3 * (m1 - m0).  Its zero is the w3 below which
            # the wrong topology has the larger expected score and the
            # estimator is no longer consistent -- a hard boundary, not a
            # gradual loss of power.
            def margin(w):
                V = vote_table(w)
                return float((V[:, 0] - V[:, 1]) @ p0)

            m0, m1 = margin(0.0), margin(1.0)
            wc = -m0 / (m1 - m0) if m1 != m0 else float("nan")
            rows.append(dict(zone=zone, L=L, ml_power=ml,
                             caster_w3_half=power(vote_table(-0.5), p0, L),
                             best_w3=best_w3, best_power=best_pr,
                             margin_w3_0=m0,
                             margin_w3_minus_half=margin(-0.5),
                             consistency_w3=wc))
        print(f"  {zone} zone: ML power {ml:.4f} at L=1e6")

    d = pd.DataFrame(rows)
    d.to_csv(OUT / "weights.csv", index=False)

    print(f"\n  {'zone':>12} {'L':>9} {'ML':>7} {'CASTER(-1/2)':>13} "
          f"{'best w3':>8} {'best P':>7} {'w3 at margin 0':>15}")
    for _, r in d.iterrows():
        print(f"  {r.zone:>12} {int(r.L):9d} {r.ml_power:7.4f} "
              f"{r.caster_w3_half:13.4f} {r.best_w3:8.2f} {r.best_power:7.4f} "
              f"{r.consistency_w3:15.3f}")

    print("\n  margin E[S1 - S2] on the Farris-zone tree, per site:")
    pF = m.pattern_batch(FARRIS, m.SPLITS["12|34"])[0]
    for w3 in (-1.5, -1.0, -0.75, -0.5, -0.25, 0.0, 0.5):
        V = vote_table(w3)
        print(f"    w3={w3:+.2f}   E[S1-S2]={float((V[:, 0] - V[:, 1]) @ pF):+.6e}"
              f"   E[S1-S3]={float((V[:, 0] - V[:, 2]) @ pF):+.6e}")


if __name__ == "__main__":
    main()
