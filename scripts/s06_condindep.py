"""Experiment 6 -- does Iss carry information about Y beyond mean p-distance?

    Pr(Y | pbar, Iss)   vs   Pr(Y | pbar)

Two modes:
    python scripts/s06_condindep.py            # simulate the crossed design
    python scripts/s06_condindep.py --analyze  # fit and compare the models

WHY THE DESIGN IS CROSSED.  If the scale s is the only knob, Iss and pbar are
deterministically linked (to second order they are an affine map), so within a
pbar bin there is nothing to condition on and the test is vacuous.  Power comes
from *different configurations* that land on the same pbar but on different Iss:
tree shape (internal-branch fraction, pendant arrangement), composition, and
model.  All three are varied here.

TARGETS.  Y1 = parsimony-informative site fraction (exact, no inference).
Y2 = recovery under the uncorrected four-point rule (in-process, cheap).
ML recovery can be added by pointing Y2 at issbench.infer.ptrue(..., "ML").
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench import trees                                             # noqa: E402
from issbench.models import build_Q                                    # noqa: E402
from issbench.runner import append_rows, dump_meta, load_done, out_path  # noqa: E402
from issbench.sim import simulate, cell_key                            # noqa: E402
from issbench.stats import parsimony_informative_fraction, read_fasta  # noqa: E402
from issbench.trees import write_newick                                # noqa: E402

EXP = "exp6_condindep"
TREES = ROOT / "results" / EXP / "trees"
N_REPS = 150
L = 1000

FIS = [0.08, 0.15, 0.35, 0.80]
LADDER = [1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 3.0, 4.0]      # pendant ratios per shape
CONDS = [
    ("JC69",        "JC69",                    [0.25, 0.25, 0.25, 0.25]),
    ("HKY_uni",     "HKY+F{0.25,0.25,0.25,0.25}+G4{1.0}", [0.25] * 4),
    ("HKY_at",      "HKY+F{0.35,0.15,0.15,0.35}+G4{1.0}", [0.35, 0.15, 0.15, 0.35]),
    ("HKY_ex",      "HKY+F{0.40,0.10,0.10,0.40}+G4{2.0}", [0.40, 0.10, 0.10, 0.40]),
    ("JC_a05",      "JC69+G4{0.5}",            [0.25] * 4),
]
S_GRID = np.exp(np.linspace(np.log(0.25), np.log(25.0), 18))


def shape(fi: float, ratio: float, TL: float):
    """Quartet with internal = fi*TL and pendants weighted by `ratio`."""
    r = fi * TL
    rest = TL - r
    if ratio == 1.0:
        p = [rest / 4.0] * 4
    else:
        w = np.array([1.0, ratio, 1.0, ratio])
        p = list(rest * w / w.sum())
    return trees.quartet(*p, r)


def cells():
    TREES.mkdir(parents=True, exist_ok=True)
    for cname, model, pi in CONDS:
        for fi in FIS:
            for ratio, rlab in ((1.0, "eq"), (3.0, "lad")):
                for s in S_GRID:
                    tag = f"{cname}_{fi:g}_{rlab}_{s:.2f}"
                    f = TREES / f"{tag}.nwk"
                    node = shape(fi, ratio, float(s))
                    write_newick(node, f)
                    yield dict(tag=tag, tree_file=f,
                               tree_nwk=trees.to_newick(node),
                               model=model, L=L, n_reps=N_REPS,
                               cond=cname, fi=fi, ratio=ratio, s=float(s))


def run_sim():
    import pandas as pd
    from issbench.infer import TRUE_SPLIT, four_point_split
    from issbench.sim import score_cell

    dump_meta(EXP, dict(item=6, n_reps=N_REPS, L=L, fis=FIS,
                        conds=[c[0] for c in CONDS],
                        s_grid=[float(x) for x in S_GRID],
                        targets=["pi_site_fraction", "Ptrue_4P"],
                        note="crossed design: shape x composition x model x scale"))
    d = out_path(EXP).parent
    path = out_path(EXP)
    done = load_done(path)
    total = 0
    for c in cells():
        todo = [r for r in range(N_REPS) if (c["tag"], r) not in done]
        if not todo:
            continue
        seed = 1 + (int(cell_key(**c), 16) % 900000)
        wd = ROOT / "work" / EXP / c["tag"]
        try:
            files = simulate(Path(c["tree_file"]), c["model"], L, len(todo),
                             seed, wd, tag=c["tag"])
        except Exception as e:                                        # noqa: BLE001
            print(f"  [{c['tag']}] FAILED: {e}", flush=True)
            continue

        rows = score_cell(files)
        targets = []
        for r, f, rep in zip(rows, files, todo):
            r.update(exp=EXP, cell=c["tag"], rep=rep, seed=seed,
                     tree=c["tree_nwk"], model=c["model"], L=L)
            _n, M = read_fasta(f)
            targets.append(dict(exp=EXP, cell=c["tag"], rep=rep,
                                cond=c["cond"], fi=c["fi"], ratio=c["ratio"],
                                s=c["s"],
                                pi_site=parsimony_informative_fraction(M),
                                Ptrue_4P=int(four_point_split(M) == TRUE_SPLIT)))
        append_rows(path, rows)
        ex = d / "targets.csv"
        pd.DataFrame(targets).to_csv(ex, mode="a", header=not ex.exists(),
                                     index=False)
        total += len(rows)
        print(f"  [{c['tag']}] {len(rows)} reps", flush=True)
    print(f"\n{EXP}: {total} new replicate rows")


def analyze():
    """Fit logit(Y) ~ f(pbar) and ~ f(pbar) + f(Iss) on held-out halves."""
    import pandas as pd
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_score
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import PolynomialFeatures, StandardScaler

    d = out_path(EXP).parent
    cells_df = pd.read_csv(d / "cells.csv")
    targ = pd.read_csv(d / "targets.csv")
    df = targ.merge(cells_df[["cell", "rep", "pbar", "iss_A", "iss_B"]],
                    on=["cell", "rep"], how="inner")
    print(f"merged {len(df)} rows")

    def fit(cols, y):
        X = df[cols].to_numpy(float)
        m = make_pipeline(PolynomialFeatures(2, include_bias=False),
                          StandardScaler(),
                          LogisticRegression(max_iter=5000))
        return cross_val_score(m, X, df[y].to_numpy(int), cv=5,
                               scoring="roc_auc").mean()

    out = []
    for y in ("Ptrue_4P",):
        auc_p = fit(["pbar"], y)
        auc_pi = fit(["pbar", "iss_A"], y)
        out.append(dict(target=y, auc_pbar=auc_p, auc_pbar_plus_iss=auc_pi,
                        delta=auc_pi - auc_p))
        print(f"  {y}: AUC(pbar)={auc_p:.4f}  AUC(pbar,Iss)={auc_pi:.4f}"
              f"  delta={auc_pi-auc_p:+.4f}")
    # Continuous target, treated the SAME way as the binary one: compare
    # out-of-sample R^2 of  f(pbar)  against  f(pbar) + Iss, with polynomial
    # features so that pbar enters flexibly.  Two ad-hoc conditioning
    # statistics (binned residual vs partial correlation) disagreed in sign on
    # this data, which is why the fit comparison is the one that is reported.
    from sklearn.linear_model import RidgeCV
    from sklearn.model_selection import KFold, cross_val_score
    from scipy.stats import spearmanr

    y = df["pi_site"].to_numpy(float)
    r_p_spear = spearmanr(df["pbar"], y).correlation

    def r2(cols):
        X = df[cols].to_numpy(float)
        m = make_pipeline(PolynomialFeatures(2, include_bias=False),
                          StandardScaler(), RidgeCV(alphas=np.logspace(-3, 3, 13)))
        kf = KFold(5, shuffle=True, random_state=0)
        return float(cross_val_score(m, X, y, cv=kf, scoring="r2").mean())

    r2_p = r2(["pbar"])
    r2_pi = r2(["pbar", "iss_A"])
    print(f"  pi_site: rho(pbar)={r_p_spear:+.4f}  CV R2(pbar)={r2_p:.4f}"
          f"  CV R2(pbar,Iss)={r2_pi:.4f}  delta={r2_pi-r2_p:+.4f}")
    out.append(dict(target="pi_site_fraction", auc_pbar=r2_p,
                    auc_pbar_plus_iss=r2_pi, delta=r2_pi - r2_p))
    pd.DataFrame(out).to_csv(d / "condindep.csv", index=False)
    print(f"wrote {d/'condindep.csv'}")


if __name__ == "__main__":
    if "--analyze" in sys.argv:
        analyze()
    else:
        run_sim()
