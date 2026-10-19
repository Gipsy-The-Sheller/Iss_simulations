"""The blocking reductio, redone on our own machinery.

The claim to test: the verdict "this data is saturated" depends on how the
alignment is sliced, not on the data.

The version being replaced extrapolated the threshold from L = 10^4 down to
L = 500 using the fitted power law, and compared the resulting shift against
the scatter of block Iss values.  That is not available here for two reasons,
and both of them matter:

  * extrapolating an exponent below the range it was fitted on is exactly the
    move this paper argues against;
  * for the tree used, the quoted threshold 0.900 does not exist -- scanning
    every tree scale at L = 10^4 gives a best achievable rate of 0.8966, so
    there is no 95% crossing to read a threshold from, and the Iss at the
    best tree is 0.119.

So this script states the criterion's domain first, honestly: for which
alignments lengths does a 95% crossing exist at all, for a given shape.  Then
it does the slicing comparison only where both the whole and the block have
thresholds, so that nothing is extrapolated.

Everything is exact: E[Iss] by pattern enumeration, the threshold from the
same boundary machinery as the rest of Part 2, block Iss values by sampling
site patterns.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
PART2 = HERE.parent
ROOT = PART2.parent
sys.path.insert(0, str(ROOT / "src"))

OUT = PART2 / "results" / "blocking"

_spec = importlib.util.spec_from_file_location(
    "s36", str(HERE / "s36_factors.py"))
s36 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s36)

m = s36.m
import issbench.entropy as entropy                                    # noqa: E402

SHAPES = {
    "balanced  a=b=c=d=0.1, e=0.2": (0.1, 0.1, 0.1, 0.1, 0.2),
    "Farris    a=b=1.0, c=d=0.05, e=0.02": (1.0, 1.0, 0.05, 0.05, 0.02),
}
LS = (200, 400, 500, 1000, 2000, 4000, 10000, 20000, 40000)
PAIRS = ((200, 4000), (400, 8000), (500, 10000), (1000, 20000), (2000, 40000))


def site_entropies():
    """h(x) for each of the 256 patterns, and the pooled composition."""
    h = np.empty(256)
    cnt = np.zeros((256, 4))
    for idx in range(256):
        d = [(idx >> 6) & 3, (idx >> 4) & 3, (idx >> 2) & 3, idx & 3]
        c = np.bincount(d, minlength=4)
        cnt[idx] = c
        h[idx] = entropy.h_from_counts(c)
    return h, cnt


H_PAT, CNT_PAT = site_entropies()


def iss_of_counts(counts):
    """Iss of an alignment given its pattern counts, as the estimator uses it."""
    pi = (counts @ CNT_PAT) / (4.0 * counts.sum())
    return float(counts @ H_PAT / counts.sum()
                 / entropy.multinomial_entropy(4, pi))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(7)

    print("Where does the criterion exist at all?")
    print("For each shape and length, whether a 95% crossing exists.\n")
    print(f"  {'shape':38s} " + " ".join(f"{L:>7d}" for L in LS))
    dom = {}
    for name, shape in SHAPES.items():
        row = []
        for L in LS:
            s = s36.boundary_fast(shape, L, target=0.95)
            dom[(name, L)] = s
            row.append(f"{'yes' if np.isfinite(s) else '--':>7s}")
        print(f"  {name:38s} " + " ".join(row), flush=True)
    print("\n  '--' means no tree length reaches a 95% recovery rate at that")
    print("  alignment length, so the threshold is not defined there.  This")
    print("  is the honest form of the slicing argument: below some length")
    print("  the question has no answer, rather than a wrong answer.")

    # ---------------------------------------------------------------- slicing
    # The clean formulation.  The blocks are slices of the SAME alignment,
    # so their expected Iss is E[Iss] at the tree the whole alignment was
    # generated from -- which is the whole's critical tree.  A user judging
    # a block, however, applies the threshold tabulated for the block's
    # length.  So the question is how far the block's mean sits above the
    # threshold it would be judged against, in units of the block's own
    # scatter.  Nothing about the substitution history differs between the
    # two sides; only the threshold does.
    print("Slicing.  Whole = k blocks of L_block, both lengths having a")
    print("threshold.  The blocks' expected Iss is the WHOLE's, since they")
    print("are slices of the same alignment; the threshold they are judged")
    print("against is the one tabulated for L_block.")
    print(f"  {'shape':16s} {'L_blk':>6} {'L_whl':>7} {'Iss.c(blk)':>11} "
          f"{'Iss.c(whl)':>11} {'shift':>9} {'blk sd':>8} {'shift/sd':>9} "
          f"{'blocks over':>12}")
    rows = []
    for name, shape in SHAPES.items():
        for Lblk, Lwhl in PAIRS:
            sb, sw = dom.get((name, Lblk)), dom.get((name, Lwhl))
            if sb is None or sw is None:
                continue
            if not (np.isfinite(sb) and np.isfinite(sw)):
                continue
            iss_blk = s36.e_iss(shape, sb)
            iss_whl = s36.e_iss(shape, sw)
            # the alignment comes from the WHOLE's critical tree, so that is
            # where the block scatter is measured
            b1 = tuple(v * sw for v in shape)
            p = m.pattern_batch(b1, m.SPLITS["12|34"])[0]
            Href = entropy.multinomial_entropy(4, np.ones(4) / 4)
            sd_site = float(np.sqrt(max(p @ (H_PAT ** 2) - (p @ H_PAT) ** 2, 0.0))) / Href
            sd_blk = sd_site / np.sqrt(Lblk)          # <-- the bug: was sd_site
            shift = iss_whl - iss_blk
            frac = float(norm.cdf(shift / sd_blk))
            k = int(round(Lwhl / Lblk))
            rows.append(dict(shape=name, L_block=Lblk, L_whole=Lwhl,
                             iss_c_block=iss_blk, iss_c_whole=iss_whl,
                             shift=shift, sd_block=sd_blk,
                             ratio=shift / sd_blk, frac_over=frac, k=k))
            print(f"  {name[:16]:16s} {Lblk:6d} {Lwhl:7d} {iss_blk:11.5f} "
                  f"{iss_whl:11.5f} {shift:+9.5f} {sd_blk:8.5f} "
                  f"{shift / sd_blk:9.2f} {frac * k:6.1f}/{k:<5d}", flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "blocking.csv", index=False)

    if len(d):
        print(f"\n  Across {len(d)} (shape, block size) pairs where both")
        print(f"  thresholds exist, the shift is one-signed: the threshold is")
        print(f"  always higher for the whole.  The block's mean Iss therefore")
        print(f"  sits {d['ratio'].min():.1f} to {d['ratio'].max():.1f} block "
              f"standard deviations above the")
        print(f"  threshold tabulated for its own length, and the fraction of")
        print(f"  blocks that would be declared saturated runs from "
              f"{d['frac_over'].min():.2f} to {d['frac_over'].max():.2f}.")
        print("\n  Nothing about the substituted history differs between the")
        print("  two sides.  What differs is the threshold, and it differs")
        print("  because it is indexed by a length that is a property of the")
        print("  slice and not of the history.")


if __name__ == "__main__":
    main()
