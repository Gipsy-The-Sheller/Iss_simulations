"""Item 4 -- one parameter at a time.

No composite power index.  Each knob is swept on its own, everything else held
fixed, and Iss / mean p-distance / C value / lambda are reported side by side so
the trends can be read directly.

Knobs: alpha (Gamma shape), kappa (ts/tv), composition skew, internal-branch
fraction, sequence length.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench import trees                                             # noqa: E402
from issbench.runner import dump_meta, run_cells                       # noqa: E402
from issbench.trees import star, write_newick                          # noqa: E402

EXP = "exp4_perturb"
TREES = ROOT / "results" / EXP / "trees"
N_REPS = 200
L = 1000
# Total tree length must sit in the RESPONSIVE range.  For a 4-taxon tree
# TL = 8 gives mean pairwise distance 3(8-r)/6 + 4r/6 ~= 4.5, i.e. pbar ~ 0.75:
# every arm is saturated and Iss is pinned near 1 whichever knob moves.
# TL ~ 2 lands at pbar ~ 0.6, where the knobs actually separate.
TL = 2.0

ALPHAS = [np.inf, 4.0, 2.0, 1.0, 0.5, 0.25]
KAPPAS = [1.0, 2.0, 5.0, 10.0, 20.0]
PIS = [("uniform", [0.25, 0.25, 0.25, 0.25]),
       ("mild", [0.30, 0.20, 0.20, 0.30]),
       ("atrich", [0.35, 0.15, 0.15, 0.35]),
       ("extreme", [0.40, 0.10, 0.10, 0.40])]
FIS = [0.05, 0.10, 0.20, 0.35, 0.60]
TLS = [0.5, 1.0, 1.5, 2.0, 3.0, 4.5]   # TL sweep, same fixed shape otherwise
LS = [100, 300, 1000, 3000, 10000]


def _q4(fi, TL):
    r = fi * TL
    l = (TL - r) / 4.0
    return trees.quartet(l, l, l, l, r)


def cells():
    TREES.mkdir(parents=True, exist_ok=True)

    def emit(tag, node, model, L_):
        f = TREES / f"{tag}.nwk"
        write_newick(node, f)
        return dict(tag=f"{tag}_L{L_}", tree_file=f, tree_nwk=trees.to_newick(node),
                    model=model, L=L_, n_reps=N_REPS)

    base = _q4(0.35, TL)
    for a in ALPHAS:
        lab = "inf" if np.isinf(a) else str(a)
        # the homogeneous arm must NOT carry a +G4 modifier: AliSim rejects
        # "-m JC69+G4{inf}" outright.
        model = "JC69" if np.isinf(a) else f"JC69+G4{{{a}}}"
        yield emit(f"alpha{lab}", base, model, L)
    for k in KAPPAS:
        yield emit(f"kappa{k}", base, f"K80{{{k}}}", L)
    for lab, pi in PIS:
        s = ",".join(f"{x:.2f}" for x in pi)
        yield emit(f"pi_{lab}", base, f"HKY+F{{{s}}}+G4{{1.0}}", L)
    for fi in FIS:
        yield emit(f"fi{fi}", _q4(fi, TL), "JC69", L)
    for L_ in LS:
        yield emit("Lsweep", base, "JC69", L_)
    for tl in TLS:
        yield emit(f"TL{tl}", _q4(0.35, tl), "JC69", L)


if __name__ == "__main__":
    dump_meta(EXP, dict(item=4, TL=TL, n_reps=N_REPS, L=L,
                        TLs=TLS,
                        alphas=[str(a) for a in ALPHAS], kappas=KAPPAS,
                        pis=[p[0] for p in PIS], fis=FIS, Ls=LS,
                        note="one knob at a time; no composite power index"))
    n = run_cells(EXP, cells(), ROOT)
    print(f"\n{EXP}: {n} new replicate rows")
