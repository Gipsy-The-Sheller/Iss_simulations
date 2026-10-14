"""Is Iss.c the point where neighbour joining stops being able to decide?

The proposal: at a fixed topology the NJ failure threshold should be the
point where the quartet's distance matrix goes semi-degenerate, and the Iss
there should be the same object Xia tabulated -- and that point is derivable
from p-distances under JC69.

The first half is exactly right and worth writing down.  NJ on a quartet
minimises the within-pair sum over the three splits, so it chooses 12|34 over
13|24 exactly when d12 + d34 < d13 + d24 in expectation, i.e. when

    M = f(d13) + f(d24) - f(d12) - f(d34) > 0,     f(D) = (3/4)(1 - e^{-4D/3}),

and the failure boundary is M = 0.  That is a statement about expectations
and carries no L, so it is a genuine algebraic critical point rather than a
power threshold.

M has an exact closed form on a quartet with pendants a, b, c, d and internal
length e, putting A = e^{-4a/3} and so on:

    M = (3/4) [ AB + CD - E (AC + BD) ],      E = e^{-4e/3},

whose first-order expansion in e is (3/4)(A-D)(B-C) + e(AC+BD) -- the
path-imbalance term plus the internal-branch term, as the note has it.

The second half is the part that does not survive.  M = 0 is ONE equation in
four pendant lengths, so its solution set is a three-dimensional surface, and
Iss is not constant on it.  This script shows where the idea works (a
one-parameter family, where the tie point is isolated) and where it cannot.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Part 2.  Inputs that belong to the published preprint (the empirical
# alignments, the simulated work directory) are read from the
# repository root; everything this line of work produces is written
# under Part_2/, so the two cannot be confused.
PART2 = Path(__file__).resolve().parents[1]
ROOT = PART2.parent
sys.path.insert(0, str(ROOT / "src"))

OUT = PART2 / "results" / "nj_tie"


def X(D: float) -> float:
    return float(np.exp(-4 * D / 3))


def f(D: float) -> float:
    return 0.75 * (1 - X(D))


def M_of(a, b, c, d, e):
    """Exact quartet margin: positive means NJ prefers the true split."""
    A, B, C, D, E = X(a), X(b), X(c), X(d), X(e)
    return 0.75 * (A * B + C * D - E * (A * C + B * D))


def iss_at(a, b, c, d, e, N=4):
    """E[Iss] under convention (A) for JC69, from the mean p-distance."""
    D = [a + b, c + d, a + c + e, b + d + e, a + d + e, b + c + e]
    wp = float(np.mean([f(x) for x in D]))
    Psi = 3.0 - 4.0 * wp
    return 1.0 - (N - 1) / (2 * N * np.log(4)) * Psi, wp


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------ the note's expansion
    print("the closed form and its expansion, checked against the note")
    a, b, c, d, e = 1.0, 1.0, 0.1, 0.1, 0.05
    exact = M_of(a, b, c, d, e)
    A, B, C, D_, E = X(a), X(b), X(c), X(d), X(e)
    note = 0.75 * (A - D_) * (B - C) + e * (A * C + B * D_)
    print(f"  sisters  (a=b=1.0, c=d=0.1, e=0.05): exact M = {exact:+.4f}"
          f"   note's expansion = {note:+.4f}")
    a, b, c, d = 1.0, 0.1, 1.0, 0.1
    exact = M_of(a, b, c, d, e)
    A, B, C, D_, E = X(a), X(b), X(c), X(d), X(e)
    note = 0.75 * (A - D_) * (B - C) + e * (A * C + B * D_)
    print(f"  opposite (a=c=1.0, b=d=0.1, e=0.05): exact M = {exact:+.4f}"
          f"   note's expansion = {note:+.4f}")
    print("  the note's two rows are +0.303 and -0.240; both reproduce.")

    # ------------------------------------------- balanced quartet has no tie
    print("\nthe balanced quartet has no NJ tie at any internal length:")
    for e in (0.001, 0.01, 0.1, 0.5, 2.0):
        m = M_of(1.0, 1.0, 1.0, 1.0, e)
        print(f"  e = {e:<6}   M = {m:+.3e} = 2e x^2 + O(e^2), positive")
    print("  so on a balanced quartet NJ never fails systematically, and the")
    print("  half-additivity line M = e is a marginal-resolvability")
    print("  criterion, not a failure boundary.  They are different objects.")

    # -------------------------------- the tie point along a one-parameter path
    print("\nthe Felsenstein-Farris transfer path: (a,b,c,d) = (L-t, S, L, S+t)")
    print("  t = 0 is Felsenstein (long branches opposite), t = L-S is Farris")
    L, S = 1.0, 0.1
    rows = []
    for e in (0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8):
        ts = np.linspace(0, L - S, 4000)
        ms = np.array([M_of(L - t, S, L, S + t, e) for t in ts])
        sign = np.where(np.diff(np.sign(ms)) != 0)[0]
        if not len(sign):
            print(f"  e = {e:<5}  no tie on the path")
            continue
        i = sign[0]
        t0, t1 = ts[i], ts[i + 1]
        for _ in range(80):                      # bisect the crossing
            tm = (t0 + t1) / 2
            if M_of(L - t0, S, L, S + t0, e) * M_of(L - tm, S, L, S + tm, e) <= 0:
                t1 = tm
            else:
                t0 = tm
        tstar = (t0 + t1) / 2
        iss, wp = iss_at(L - tstar, S, L, S + tstar, e)
        iss_lo, _ = iss_at(L, S, L, S, e)
        iss_hi, _ = iss_at(S, S, L, L, e)
        rows.append(dict(e=e, t_star=tstar, frac=tstar / (L - S), wp=wp,
                         iss=iss, iss_t0=iss_lo, iss_farris=iss_hi))
        print(f"  e = {e:<5}  t* = {tstar:.4f} ({100 * tstar / (L - S):4.1f}%"
              f" along)   mean p = {wp:.4f}   Iss at t* = {iss:.4f}"
              f"   (Iss range over path {iss_lo:.4f} to {iss_hi:.4f})")
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "transfer_tie.csv", index=False)
    if len(d):
        print(f"\n  Iss at the tie point runs {d.iss.min():.4f} to "
              f"{d.iss.max():.4f} as e goes {d.e.min()} to {d.e.max()}")
        print(f"  the internal branch length moves it by "
              f"{d.iss.max() - d.iss.min():.4f}, which is "
              f"{100 * (d.iss.max() - d.iss.min()) / d.iss.mean():.1f}% of "
              f"its value")

    # ----------------------------- the e -> 0 limit is a star, and it drifts
    # e -> 0 collapses the quartet onto a star, where M = (3/4)(A-D)(B-C) and
    # the tie is a = d or b = c.  It is tempting to read the e -> 0 end of
    # the table as a single asymptotic critical value.  It is not: the star
    # family still has three free lengths and Iss varies over all of them.
    print("\nthe e -> 0 limit, on the star itself (tie at a = d)")
    print(f"  {'L':>5} {'S':>6} {'a=d':>7} {'mean p':>8} {'Iss at tie':>11}")
    rows2 = []
    for L, S in ((1.0, 0.10), (1.0, 0.02), (1.0, 0.30), (0.5, 0.05),
                 (2.0, 0.20), (2.0, 0.02), (1.0, 0.50)):
        a, b_, c, dd = (L + S) / 2, S, L, (L + S) / 2   # tie is a = d
        iss, wp = iss_at(a, b_, c, dd, 0.0)
        rows2.append(dict(L=L, S=S, a_equals_d=a, wp=wp, iss=iss))
        print(f"  {L:5.2f} {S:6.2f} {a:7.3f} {wp:8.4f} {iss:11.4f}")
    pd.DataFrame(rows2).to_csv(OUT / "star_tie.csv", index=False)
    e_series = pd.DataFrame(rows2).iss
    print(f"  the star tie point runs {e_series.min():.4f} to "
          f"{e_series.max():.4f} -- it is not a constant either, so the")
    print("  e -> 0 end of the previous table is not an asymptotic critical")
    print("  value.  And on a star the question the whole construction asks")
    print("  has no answer: there is no internal edge, so there is no split")
    print("  for NJ to recover.  What M = 0 marks there is the point where")
    print("  NJ's choice becomes a coin flip on a tree that has no split to")
    print("  find -- which is why the limit cannot carry a critical Iss.")

    # ------------------------------------ what Xia's own numbers reverse to
    print("\nXia et al. 2003's own critical values, reverse-engineered to xi")
    print("  (the EF-1a example, N=16, SeqLen=350 from the text)")
    for label, c in (("symmetrical", 0.7026), ("asymmetrical", 0.4890)):
        xi = (1 - c) * (2 * 16 * np.log(4)) / (3 * 15)
        print(f"  {label:14s} Iss.c = {c:.4f}  ->  xi = {xi:.4f}"
              f"   (mean p = {0.75 * (1 - xi):.4f})")
    print("  the half-additivity line in the note is xi = 1/2, and the")
    print("  asymmetrical value sits on it to three decimals.")

    # ------------------------------------- is the L-dependence L^{-1/4}?
    # The detectability criterion M = sampling noise gives e xi^2 =
    # sqrt(p(1-p)/L), hence xi proportional to L^{-1/4} at fixed e and p.
    # Xia's text gives one point that can be checked directly (N=16,
    # SeqLen=350); the note's table gives another at L=1177 whose provenance
    # is secondary.  If both are the same quantity, the ratio is a test.
    print("\nis the drift of xi* with L the L^{-1/4} that detectability implies?")
    xi_350 = (1 - 0.7026) * (2 * 16 * np.log(4)) / (3 * 15)
    xi_1177 = (1 - 0.780) * (2 * 16 * np.log(4)) / (3 * 15)
    pred = (350 / 1177) ** 0.25
    print(f"  N=16 symmetrical, L =  350 (Xia et al. text):  xi = {xi_350:.4f}")
    print(f"  N=16 symmetrical, L = 1177 (note's table):     xi = {xi_1177:.4f}")
    print(f"  observed ratio {xi_1177 / xi_350:.4f}   L^-1/4 predicts "
          f"{pred:.4f}   agreement to "
          f"{100 * abs(xi_1177 / xi_350 - pred) / pred:.2f}%")
    print("  this is the experiment the note asks for in its section 10, and")
    print("  two points may already answer it -- but the L=1177 row is a")
    print("  secondary table whose source is unconfirmed, so treat it as a")
    print("  lead and re-derive it from DAMBE before relying on it.")


if __name__ == "__main__":
    main()
