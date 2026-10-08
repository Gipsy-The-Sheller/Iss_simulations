"""Verify the gpt6 review's counterexamples with the project's exact machinery.

R-3.2  two 4x4 alignments with equal mean p-distance but different Iss
R-3.5  exact E[Iss] on a JC69 star, N=4, at d = 0.1 / 0.4 / 1.0
R-7    Iss > 1 with no gaps: every site 1:1:1:1
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench import exact, models, trees                              # noqa: E402
from issbench.entropy import h, multinomial_entropy                    # noqa: E402
from issbench.stats import iss, mean_p_distance                        # noqa: E402

CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def align(rows) -> np.ndarray:
    """rows: list of site strings, four taxa left to right."""
    return np.array([[CODE[c] for c in r] for r in rows]).T


def show(tag, M):
    print(f"  {tag:<12} pbar={mean_p_distance(M):.6f}  hbar={-np.mean([np.mean([0])]):.0f}"
          if False else
          f"  {tag:<12} pbar={mean_p_distance(M):.6f}  "
          f"Iss_A={iss(M,'A'):.6f}  Iss_B={iss(M,'B'):.6f}")


HR = "=" * 74
print(HR + "\nR-3.2   equal mean p-distance, different Iss\n" + HR)
A = align(["AAAC", "CCCG", "GGGT", "TTTA"])
B = align(["CCGG", "GGTT", "TTCC", "AAAA"])
show("alignment A", A)
show("alignment B", B)
print(f"  review claims  Iss_A=0.612753  Iss_B=0.566470  (diff 0.0463)")

print("\n" + HR + "\nR-3.5   exact E[Iss] vs mean p-distance, JC69 star, N=4, conv A\n" + HR)
Q = models.build_Q("JC69")
PI = np.ones(4) / 4
print(f"  {'d':>6} {'E[pbar]':>10} {'E[Iss] exact':>14}")
prev = None
pts = []
for d in (0.1, 0.4, 1.0):
    e_iss = exact.E_Iss(4, d, Q, PI, convention="A")
    e_pbar = 0.75 * (1 - np.exp(-8 * d / 3))
    pts.append((e_pbar, e_iss))
    print(f"  {d:>6.1f} {e_pbar:>10.5f} {e_iss:>14.5f}")
s1 = (pts[1][1] - pts[0][1]) / (pts[1][0] - pts[0][0])
s2 = (pts[2][1] - pts[1][1]) / (pts[2][0] - pts[1][0])
print(f"  secant slopes: {s1:.4f} and {s2:.4f}   (review: 1.3005 / 1.4111)")
print(f"  review's Iss values: 0.21788 / 0.62927 / 0.91995")

print("\n" + HR + "\nR-7   Iss > 1 with no gaps: every site 1:1:1:1\n" + HR)
M = align(["ACGT"] * 20)
show("all 1:1:1:1", M)
print(f"  review claims Iss = 1.510586")
print(f"  H_FSS^(A) at N=4 uniform = {multinomial_entropy(4, PI):.6f}")
print(f"  exact saturation ceiling  = {multinomial_entropy(4, PI)/h(PI):.6f} (conv B)")
print(f"  E[h] under Mult(4,unif)   = {multinomial_entropy(4, PI):.6f}")

print("\n  how often does the 1:1:1:1 pattern occur at true saturation?")
rng = np.random.default_rng(0)
n = rng.multinomial(4, PI, size=200000)
allfour = (n > 0).sum(axis=1) == 4
print(f"    P(all four states present) = {allfour.mean():.4f}")
print(f"    E[Hbar] at saturation      = "
      f"{(-np.sum(np.where(n>0,(n/4)*np.log(np.where(n>0,n/4,1)),0),axis=1)).mean():.6f}")
