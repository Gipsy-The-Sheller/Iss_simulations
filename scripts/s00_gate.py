"""G0 -- pipeline gate.  Nothing downstream runs until every check here passes.

Checks
  1  AliSim executes at all, and --seed makes it bit-reproducible
  2  --branch-scale is multiplicative (else we write scaled Newicks ourselves)
  3  measured mean p-distance matches the JC69 closed form
  4  measured Iss matches Part I Theorem 2 exactly, under BOTH conventions
  5  which Gamma discretisation AliSim's +G actually uses (mean vs median)
  6  exact.py: star closed form == small-tree pattern enumeration
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench import entropy, exact, gamma, trees                       # noqa: E402
from issbench.alisim import run_alisim, iqtree_path, alisim_help        # noqa: E402
from issbench.models import build_Q, spectrum                           # noqa: E402
from issbench.stats import iss, mean_p_distance, read_fasta             # noqa: E402

WORK = ROOT / "results" / "gate"
WORK.mkdir(parents=True, exist_ok=True)
RESULT: list[tuple[str, bool, str]] = []


def record(name: str, ok: bool, detail: str = ""):
    RESULT.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))


def md5(p: Path) -> str:
    return hashlib.md5(p.read_bytes()).hexdigest()


# --------------------------------------------------------------------- check 1
def check_alisim_runs():
    tree = trees.star(16, 0.2)
    trees.write_newick(tree, WORK / "star16.nwk")
    outs = []
    for trial in (1, 2):
        f = run_alisim(WORK / f"seed{trial}", WORK / "star16.nwk", "JC69",
                       length=2000, n_alignments=1, seed=42, threads=1,
                       log=WORK / f"seed{trial}.log")
        outs.append(f)
    if not all(outs):
        record("1 AliSim runs", False, "no output files")
        return None
    same = md5(outs[0][0]) == md5(outs[1][0])
    record("1a AliSim runs", True, f"{outs[0]}")
    record("1b --seed reproducible", same,
           "identical md5" if same else "DIFFERENT output for the same seed")
    return outs[0][0]


# --------------------------------------------------------------------- check 2
def check_branch_scale():
    t1 = trees.star(16, 0.2)
    trees.write_newick(t1, WORK / "t1.nwk")
    t2 = trees.star(16, 0.4)
    trees.write_newick(t2, WORK / "t2.nwk")
    try:
        a = run_alisim(WORK / "scale_a", WORK / "t1.nwk", "JC69", length=20000,
                       seed=7, extra=["--branch-scale", "2"], log=WORK / "scale_a.log")
        b = run_alisim(WORK / "scale_b", WORK / "t2.nwk", "JC69", length=20000,
                       seed=7, log=WORK / "scale_b.log")
    except Exception as e:                                       # noqa: BLE001
        record("2 --branch-scale multiplicative", False, f"error: {e}")
        return
    da = mean_p_distance(read_fasta(a[0])[1])
    db = mean_p_distance(read_fasta(b[0])[1])
    ok = abs(da - db) < 0.01
    record("2 --branch-scale multiplicative", ok, f"p̄(scale 2)={da:.4f}  p̄(2x tree)={db:.4f}")


# ------------------------------------------------------------- checks 3 and 4
def check_closed_forms(d=0.2, L=200_000, rep_file=None):
    Q = build_Q("JC69")
    pi = np.ones(4) / 4
    N = 16
    f = rep_file or run_alisim(WORK / "cf", WORK / "star16.nwk", "JC69",
                               length=L, seed=11, log=WORK / "cf.log")[0]
    _, M = read_fasta(f)
    obs_p = mean_p_distance(M)
    pred_p = 0.75 * (1.0 - np.exp(-8.0 * d / 3.0))
    record("3 p-distance matches closed form", abs(obs_p - pred_p) < 5e-3,
           f"obs={obs_p:.5f}  pred={pred_p:.5f}")

    for conv in ("A", "B"):
        obs = iss(M, convention=conv)
        pred = exact.E_Iss(N, d, Q, pi, convention=conv)
        record(f"4 Iss matches Theorem 2 (convention {conv})",
               abs(obs - pred) < 5e-3, f"obs={obs:.5f}  exact={pred:.5f}  d={obs-pred:+.5f}")


# --------------------------------------------------------------------- check 5
GAMMA_D = 3.2


def check_gamma_discretisation():
    Q = build_Q("JC69")
    pi = np.ones(4) / 4
    N, alpha = 16, 1.0
    trees.write_newick(trees.star(N, GAMMA_D), WORK / "star16_d3.nwk")
    try:
        f = run_alisim(WORK / "gam", WORK / "star16_d3.nwk", f"JC69+G4{{{alpha}}}",
                       length=200_000, seed=13, log=WORK / "gam.log")[0]
    except Exception as e:                                       # noqa: BLE001
        record("5 Gamma discretisation", False, f"AliSim model string failed: {e}")
        return
    obs = iss(read_fasta(f)[1], convention="A")
    preds = {}
    for meth in ("mean", "median"):
        r = gamma.discrete_gamma_rates(alpha, 4, meth)
        preds[meth] = exact.E_Iss(N, GAMMA_D, Q, pi, rates=r,
                                  weights=gamma.gamma_weights(4), convention="A")
    best = min(preds, key=lambda m: abs(preds[m] - obs))
    record(f"5 AliSim +G4 uses the {best} discretisation",
           abs(preds[best] - obs) < 5e-3,
           f"obs={obs:.5f}  meanG4={preds['mean']:.5f}  medianG4={preds['median']:.5f}")


# ------------------------------------------------------------------- check 7/8
def check_freqs():
    """+F must actually install the requested stationary composition."""
    try:
        f = run_alisim(WORK / "freq", WORK / "star16.nwk",
                       "HKY+F{0.4,0.3,0.2,0.1}+G4{1.0}", length=200_000,
                       seed=17, log=WORK / "freq.log")[0]
    except Exception as e:                                       # noqa: BLE001
        record("7 +F frequencies honoured", False, f"AliSim failed: {e}")
        return
    from issbench.stats import pi_hat
    obs = pi_hat(read_fasta(f)[1])
    ok = np.allclose(obs, [0.4, 0.3, 0.2, 0.1], atol=5e-3)
    record("7 +F frequencies honoured", ok, f"observed {np.round(obs,4)}")


def check_gtr_rates():
    """Custom exchangeabilities (GTR) must be accepted."""
    try:
        run_alisim(WORK / "gtr", WORK / "star16.nwk",
                   "GTR{1.2,3.1,0.7,2.4,1.0}+F{0.3,0.2,0.2,0.3}+G4{0.5}",
                   length=5_000, seed=19, log=WORK / "gtr.log")
        record("8 custom GTR Q accepted", True, "R={1.2,3.1,0.7,2.4,1.0}")
    except Exception as e:                                       # noqa: BLE001
        record("8 custom GTR Q accepted", False, f"AliSim failed: {e}")


# --------------------------------------------------------------------- check 6
def check_exact_agreement():
    """Star closed form (any N) vs pattern enumeration (N=4) -- same numbers."""
    Q = build_Q("JC69")
    pi = np.ones(4) / 4
    d, r = 0.3, 0.05
    star4 = trees.star(4, d)
    q = trees.quartet(d, d, d, d, 0.0)
    # a star of 4 with leaf length d  ==  quartet with zero internal branch
    a = exact.E_Iss(4, d, Q, pi, convention="A")
    b = exact.E_Iss(4, q, Q, pi, convention="A")
    record("6 closed form == pattern enumeration (N=4)",
           abs(a - b) < 5e-6, f"star={a:.6f}  quartet(r=0)={b:.6f}")
    # and the quartet with r>0 must differ
    c = exact.E_Iss(4, trees.quartet(d, d, d, d, r), Q, pi, convention="A")
    record("6b exact.py sensitive to topology", abs(c - b) > 1e-6,
           f"r=0:{b:.6f}  r={r}:{c:.6f}")


if __name__ == "__main__":
    print(f"IQ-TREE: {iqtree_path()}\n")
    f = check_alisim_runs()
    check_branch_scale()
    check_closed_forms(rep_file=f)
    check_gamma_discretisation()
    check_exact_agreement()
    check_freqs()
    check_gtr_rates()

    print("\n" + "=" * 68)
    npass = sum(ok for _, ok, _ in RESULT)
    print(f"G0: {npass}/{len(RESULT)} checks passed")
    for name, ok, det in RESULT:
        if not ok:
            print(f"  FAILED: {name}  {det}")
    sys.exit(0 if npass == len(RESULT) else 1)
