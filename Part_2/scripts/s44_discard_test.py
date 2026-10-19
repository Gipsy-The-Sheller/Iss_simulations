"""Does discarding the saturated sites help, and by what criterion?

Duchene et al. (2022) report that excluding loci with substitution saturation
improves inferences from phylogenomic data.  In the framework of this paper
that claim has a precise form.  The decision rests on

    z = Delta sqrt(L) / sqrt(v2),

so filtering a fraction f of the sites helps if and only if

    (Delta_f / sqrt(v2_f)) sqrt(f)  >  Delta / sqrt(v2),

where the subscripted quantities are recomputed on the kept sites.  Both
sides are measurable, exactly, without simulating anything.

The filter is the one Iss implies: a site's contribution to Iss is its own
entropy h(x), so "saturated sites" are the high-h sites and excluding them is
a threshold on h.  That is not an approximation to what a user of Iss does at
locus level; at quartet scale it is the same operation.

Three outcomes are possible and they mean different things:

  * removal helps -- the discarded sites were carrying more noise than
    signal, and Duchene et al.'s practice is sound in this framework;
  * removal hurts -- the discarded sites carried the signal, and the practice
    is backwards;
  * removal helps only above some threshold, or the answer changes with the
    tree -- the benefit is not a property of saturation but of which sites
    happen to be informative, which is what this paper argues.

The last case is the one to look for, and the way to find it is to sweep the
threshold rather than test one.
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

OUT = PART2 / "results" / "discard_test"

_spec = importlib.util.spec_from_file_location(
    "s36", str(HERE / "s36_factors.py"))
s36 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s36)

m = s36.m
from issbench import entropy                                        # noqa: E402

H_PAT = np.empty(256)
for _i in range(256):
    _d = [(_i >> 6) & 3, (_i >> 4) & 3, (_i >> 2) & 3, _i & 3]
    H_PAT[_i] = entropy.h_from_counts(np.bincount(_d, minlength=4))

HMAX = float(np.log(4))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    # a Farris-zone tree at its critical scale for L = 30000
    shape = (1.0, 1.0, 0.05, 0.05, 0.02)
    s_star = s36.boundary_fast(shape, 30000, target=0.95)
    b1 = tuple(v * s_star for v in shape)
    p0 = m.pattern_batch(b1, m.SPLITS["12|34"])[0]

    print("Discarding saturated sites, on the criterion of this paper.")
    print(f"  tree a=b=1.0, c=d=0.05, e=0.02 at its critical scale "
          f"({s_star:.4f}, TL={sum(b1):.3f}) for L=30000\n")
    print(f"  {'h threshold':>11} {'kept f':>7} {'Delta_f':>11} "
          f"{'sqrt(v2_f)':>11} {'signal':>9} {'z ratio':>9}")

    rows = []
    # recompute Delta and v2 on the filtered distribution
    for thr in (HMAX, 1.30, 1.10, 0.90, 0.70, 0.58, 0.0):
        keep = H_PAT <= thr + 1e-12
        f = float(p0[keep].sum())
        if f <= 1e-9:
            continue
        pf = p0 * keep
        pf = pf / pf.sum()
        # refit the alternatives on the filtered distribution: the true
        # topology's fit is the mixture itself only when the mixture is on the
        # tree, which after filtering it is not, so topology 1 is fitted too
        fits = []
        for lab in ("12|34", "13|24", "14|23"):
            _, bb = m.best_fit(pf, m.SPLITS[lab], b1)
            fits.append(m.pattern_batch(bb, m.SPLITS[lab])[0])
        mu, cov = m.moments(pf, m.logs_of(fits))
        m2, m3, v2, v3, c23 = m.params(mu, cov)
        sig = m2 / np.sqrt(v2) if v2 > 0 else float("nan")
        rows.append(dict(thr=thr, kept=f, delta=m2, sqrt_v2=float(np.sqrt(v2)),
                         signal=sig))
        print(f"  {thr:11.4f} {f:7.4f} {m2:11.4e} {np.sqrt(v2):11.4e} "
              f"{sig:9.4f} {'--':>9}", flush=True)

    d = pd.DataFrame(rows)
    d.to_csv(OUT / "discard_test.csv", index=False)

    if len(d) > 1:
        base = d.iloc[0]
        print(f"\n  reference (nothing discarded): signal = {base.signal:.4f}")
        print(f"  {'h threshold':>11} {'kept f':>7} {'z ratio vs none':>17} "
              f"{'verdict':>10}")
        for _, r in d.iterrows():
            ratio = (r.signal / base.signal) * np.sqrt(r.kept)
            verdict = "helps" if ratio > 1.001 else (
                "hurts" if ratio < 0.999 else "flat")
            print(f"  {r.thr:11.4f} {r.kept:7.4f} {ratio:17.4f} "
                  f"{verdict:>10}")


if __name__ == "__main__":
    main()
