"""An Iss-preserving path from full recoverability to none.

Take an alignment and permute the taxa WITHIN each column, independently per
column.  This is the construction that makes the point, because of what it
does and does not change:

  * each column keeps its multiset of bases, so the per-site entropy h is
    exactly unchanged;
  * the pooled base composition pi-hat is unchanged, since each column
    contributes the same counts;
  * therefore Iss = mean(h) / H_FSS(pi-hat) is EXACTLY unchanged, not
    approximately.

  * the correspondence between taxa across columns is destroyed, so the
    phylogenetic signal is gone.

And it is not only the fully shuffled alignment that preserves Iss.  Any
mixture of shuffled and unshuffled columns does, because Iss is linear in
the site distribution.  So for

    P_f = (1 - f) P + f Q,      Q = (1/N!) sum_sigma P o sigma^{-1},

Iss(P_f) = Iss(P) for every f in [0, 1] -- a one-parameter family along
which Iss is constant to the last bit.

At f = 1, Q is invariant under every relabelling, so the three topologies
are exchangeable and recovery is exactly 1/3.  At f = 0 it is whatever the
tree gives.  Along the path, Iss does not move at all.

That is a constructive refutation of reading Iss as a saturation threshold
that governs whether the topology can be recovered, and it is stronger than
any amount of curve fitting: it does not argue that Iss and Iss.c disagree,
it exhibits two data sets with the same Iss and completely different
recoverability, with a continuous path between them.

Everything here is exact.  P and Q are pattern vectors, P_f is their
mixture, Iss is the exact pattern-enumerated expectation, and Pr[correct]
comes from the same Delta, v2 machinery as the rest of Part 2 -- with a
Monte Carlo cross-check on P_f itself.
"""
from __future__ import annotations

import importlib.util
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PART2 = HERE.parent
ROOT = PART2.parent
sys.path.insert(0, str(ROOT / "src"))

OUT = PART2 / "results" / "iss_null_path"

_spec = importlib.util.spec_from_file_location(
    "s27", str(HERE / "s27_mc_boundary.py"))
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)

from issbench import entropy                                        # noqa: E402

TRUE = m.SPLITS["12|34"]
PI = m.PI


def perm_index_arrays():
    """For each of the 24 leaf relabellings, the index map it induces.

    A pattern is packed as A,B,C,D.  Relabelling by sigma sends x to y with
    y_j = x_{sigma^-1(j)}, which is a permutation of the 256 indices.
    """
    out = []
    for sig in itertools.permutations(range(4)):
        inv = [0] * 4
        for j, s in enumerate(sig):
            inv[s] = j
        perm = np.empty(256, dtype=int)
        for idx in range(256):
            d = [(idx >> 6) & 3, (idx >> 4) & 3, (idx >> 2) & 3, idx & 3]
            e = [d[inv[j]] for j in range(4)]
            perm[idx] = (e[0] << 6) | (e[1] << 4) | (e[2] << 2) | e[3]
        out.append(perm)
    return out


PERMS = perm_index_arrays()


def permuted_average(p):
    """Q = mean over relabellings of P o sigma^{-1}."""
    q = np.zeros(256)
    for perm in PERMS:
        q[perm] += p
    return q / len(PERMS)


def e_iss_vector(p):
    """Iss from a pattern vector: mean site entropy over H_FSS(pi-hat).

    Uses issbench's own h_from_counts and multinomial entropy so the
    convention matches the rest of Part 2 exactly.
    """
    pi = np.zeros(4)
    eh = 0.0
    for idx in range(256):
        if p[idx] == 0.0:
            continue
        d = [(idx >> 6) & 3, (idx >> 4) & 3, (idx >> 2) & 3, idx & 3]
        cnt = np.bincount(d, minlength=4)
        pi += p[idx] * cnt
        eh += p[idx] * entropy.h_from_counts(cnt)
    pi /= 4.0
    return eh / entropy.multinomial_entropy(4, pi)


def recovery(p, b1, L, rng=None):
    """Pr[correct] for the mixture distribution p, analytically and by MC.

    ALL THREE topologies are fitted, including the first.  The rest of Part 2
    can take topology 1 at the truth because the truth lies on it; a mixture
    P_f lies on no topology's manifold once f > 0, and handing it to topology
    1 unfitted gives it a KL of zero while the alternatives still pay theirs
    -- which drives the rate to 1 exactly as the shuffle destroys the signal.
    At f = 0 the two conventions agree, since the KL projection of a tree
    onto its own topology is itself.
    """
    fits = []
    for lab in m.SPLITS:
        _, bb = m.best_fit(p, m.SPLITS[lab], b1)
        fits.append(m.pattern_batch(bb, m.SPLITS[lab])[0])
    a = m.analytic(p, fits, L)
    e = float("nan")
    if rng is not None:
        e = m.mc(p, fits, L, 20000, rng)
    return a, e


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(4242)

    # A tree in the Farris zone, at a length where recovery is good.
    b1 = (0.05, 0.05, 1.0, 1.0, 0.02)          # TL = 2.12
    L = 30000
    p = m.pattern_batch(b1, TRUE)[0]
    q = permuted_average(p)

    print("An Iss-preserving path from full recovery to none.")
    print("Within-column taxon permutation, f = the fraction of columns")
    print("shuffled.  Iss is meant to be exactly constant in f.\n")
    print(f"  tree a=b=0.05, c=d=1.0, e=0.02 (TL={sum(b1):.2f}), L={L}")
    print(f"  Iss(P) = {e_iss_vector(p):.8f}")
    print(f"  Iss(Q) = {e_iss_vector(q):.8f}")
    print(f"  difference = {e_iss_vector(p) - e_iss_vector(q):+.3e}")
    print("  (Q is checked, not assumed, to preserve Iss)\n")

    print(f"  {'f':>6} {'Iss(P_f)':>12} {'Pr[correct]':>12} {'MC check':>10}")
    rows = []
    for f in (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0):
        pf = (1 - f) * p + f * q
        iss = e_iss_vector(pf)
        a, e = recovery(pf, b1, L, rng)
        rows.append(dict(f=f, iss=iss, p_correct=a, mc=e))
        print(f"  {f:6.2f} {iss:12.8f} {a:12.5f} {e:10.5f}", flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "iss_null_path.csv", index=False)

    span = d.iss.max() - d.iss.min()
    print(f"\n  Iss moves by {span:.2e} across the whole path.")
    print(f"  Pr[correct] moves from {d.p_correct.iloc[0]:.4f} to "
          f"{d.p_correct.iloc[-1]:.4f}.")
    print("\n  At f = 1 the distribution is invariant under every relabelling")
    print("  of the four leaves, so the three topologies are exchangeable and")
    print("  the rate must be 1/3 -- which is a check on the machinery, not a")
    print("  result.  The result is that Iss does not move while it happens.")


if __name__ == "__main__":
    main()
