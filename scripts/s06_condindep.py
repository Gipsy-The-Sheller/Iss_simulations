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
Y2 = recovery under maximum parsimony, via MPBoot (see --recovery).  A third
target, ML recovery by IQ-TREE, can be added by pointing Y2 at
issbench.infer.ptrue(..., "ML"); it costs one IQ-TREE run per alignment and
was not run at the full 108,000.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench import trees                                             # noqa: E402
from issbench.results import res                                       # noqa: E402
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
    from issbench.sim import score_cell

    dump_meta(EXP, dict(item=6, n_reps=N_REPS, L=L, fis=FIS,
                        conds=[c[0] for c in CONDS],
                        s_grid=[float(x) for x in S_GRID],
                        targets=["pi_site_fraction", "MPB"],
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
        if len(rows) != len(todo):
            raise RuntimeError(
                f"[{c['tag']}] scored {len(rows)} alignments for {len(todo)} "
                f"pending replicates; a zip() here would mislabel the reps")
        targets = []
        for r, f, rep in zip(rows, files, todo):
            r.update(exp=EXP, cell=c["tag"], rep=rep, seed=seed,
                     tree=c["tree_nwk"], model=c["model"], L=L)
            _n, M = read_fasta(f)
            targets.append(dict(exp=EXP, cell=c["tag"], rep=rep,
                                cond=c["cond"], fi=c["fi"], ratio=c["ratio"],
                                s=c["s"],
                                pi_site=parsimony_informative_fraction(M)))
        append_rows(path, rows)
        ex = d / "targets.csv"
        pd.DataFrame(targets).to_csv(ex, mode="a", header=not ex.exists(),
                                     index=False)
        total += len(rows)
        print(f"  [{c['tag']}] {len(rows)} reps", flush=True)
    print(f"\n{EXP}: {total} new replicate rows")


def recovery_targets(methods=("MPB",)):
    """Second pass over the alignments already on disk: recovery under
    maximum parsimony as implemented by MPBoot.

    Incremental -- one row per (cell, rep, method); re-running costs only the
    missing cells.  MPBoot is a subprocess per alignment, so on Windows the
    process-spawn overhead, not the search, dominates.
    """
    import pandas as pd
    from issbench.infer import (TRUE_SPLIT, mpboot_mp, split_of_newick)

    out = out_path(EXP).parent / "recovery_targets.csv"
    cols = ["cell", "rep", *methods]
    # The file is WIDE -- one column per method -- and done-ness is per
    # (cell, rep, method), because a row may carry a value for one method and
    # be missing the other.  It used to be read as long format
    # (prev["method"]), which crashed on the first restart: the job could only
    # ever run start-to-finish in one pass.
    done = set()
    if out.exists():
        prev = pd.read_csv(out)
        for m in methods:
            if m not in prev.columns:
                continue
            ok = prev[m].notna()
            done |= {(c, int(r), m) for c, r in
                     zip(prev.loc[ok, "cell"], prev.loc[ok, "rep"])}
    wd_root = ROOT / "work" / EXP
    n_new = 0
    files = sorted(wd_root.glob("*/*_r*.fa"))
    print(f"found {len(files)} alignments under {wd_root}; "
          f"{len(done)} (cell, rep, method) targets already done", flush=True)
    for i, f in enumerate(files):
        cell = f.parent.name
        rep = int(f.stem.rsplit("_r", 1)[1])
        need = [m for m in methods if (cell, rep, m) not in done]
        if not need:
            continue
        row = {"cell": cell, "rep": rep}
        for m in need:
            s = split_of_newick(mpboot_mp(f, seed=1 + rep, workdir=f.parent))
            row[m] = int(s == TRUE_SPLIT)
        # explicit column order: a row for only one method would otherwise be
        # written under the other method's header
        pd.DataFrame([row], columns=cols).to_csv(
            out, mode="a", header=not out.exists(), index=False)
        n_new += 1
        if n_new % 500 == 0:
            print(f"  {n_new} new  ({i+1}/{len(files)})", flush=True)
    print(f"wrote {out}  ({n_new} new rows)", flush=True)


def analyze():
    """Does Iss add anything beyond pbar?  Compare out-of-sample fits.

    THE BASELINE MUST BE FLEXIBLE.  Comparing a spline-or-quadratic fit in pbar
    against the same plus Iss is only meaningful if the pbar-only fit is a fair
    summary of what pbar alone can do.  A degree-2 polynomial is not: in the
    transition band pbar in [0.63, 0.70] its cross-validated AUC is 0.37, i.e.
    *worse than chance*, because the quadratic is non-monotone and inverts the
    ordering there.  Against that baseline Iss looks like it hurts overall
    (delta = -0.010) while in fact it is repairing a mis-specified competitor.
    We therefore use a cubic spline in pbar throughout, and report the
    transition band separately from the full range, because the two answer
    different questions: over the full range pbar alone already separates
    almost perfectly and there is nothing left to add, whereas in the band
    where the recovery decision is actually made pbar alone is nearly useless.

    A third model adds the internal-branch fraction to identify the mechanism:
    it is a property of the true tree that no user observes, so it is not a
    competitor to Iss but a probe of what Iss is proxying for.
    """
    import pandas as pd
    from sklearn.linear_model import LogisticRegression, RidgeCV
    from sklearn.model_selection import KFold, cross_val_score
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import SplineTransformer, StandardScaler

    d = out_path(EXP).parent
    cells_df = pd.read_csv(res(d, "cells.csv"),
                           usecols=["cell", "rep", "pbar", "iss_A"])
    targ = pd.read_csv(res(d, "targets.csv"))
    rt = res(d, "recovery_targets.csv")
    if rt.exists():
        targ = targ.merge(pd.read_csv(rt), on=["cell", "rep"], how="left")
    ml = res(d, "recovery_ml.csv")
    if ml.exists():
        targ = targ.merge(pd.read_csv(ml), on=["cell", "rep"], how="left")
    df = targ.merge(cells_df, on=["cell", "rep"], how="inner")
    # cell tag = <cond>_<fi>_<ratio>_<s>; fi is the internal fraction
    parts = df["cell"].str.rsplit("_", n=3, expand=True)
    df["fi"] = parts[1].astype(float)
    have = [c for c in ("MPB", "ML") if c in df.columns]
    print(f"merged {len(df)} rows; recovery columns: {have}")

    KNOTS = 16
    BAND = (0.63, 0.70)

    def score(kind, cols, y, sub):
        X = sub[cols].to_numpy(float)
        yy = sub[y].to_numpy(int if kind == "auc" else float)
        if kind == "auc":
            m = make_pipeline(SplineTransformer(n_knots=KNOTS, degree=3),
                              StandardScaler(), LogisticRegression(max_iter=5000))
            return float(cross_val_score(m, X, yy, cv=5,
                                         scoring="roc_auc").mean())
        m = make_pipeline(SplineTransformer(n_knots=KNOTS, degree=3),
                          StandardScaler(),
                          RidgeCV(alphas=np.logspace(-3, 3, 13)))
        return float(cross_val_score(m, X, yy,
                                     cv=KFold(5, shuffle=True, random_state=0),
                                     scoring="r2").mean())

    band = df[(df.pbar >= BAND[0]) & (df.pbar <= BAND[1])]
    print(f"transition band pbar in {BAND}: n={len(band)}")
    out = []
    for y, kind in ([(c, "auc") for c in have]
                    + [("pi_site", "r2")]):
        n_na = int(df[y].isna().sum())
        if n_na:
            print(f"  {y}: skipped -- {n_na}/{len(df)} missing "
                  f"(recovery pass still running?)")
            continue
        full = df.dropna(subset=[y])
        bnd = band.dropna(subset=[y])
        r = dict(target=y)
        for tag, sub in (("full", full), ("band", bnd)):
            r[f"{tag}_pbar"] = score(kind, ["pbar"], y, sub)
            r[f"{tag}_pbar_iss"] = score(kind, ["pbar", "iss_A"], y, sub)
            r[f"{tag}_delta"] = r[f"{tag}_pbar_iss"] - r[f"{tag}_pbar"]
        # mechanism probe, band only
        r["band_pbar_fi"] = score(kind, ["pbar", "fi"], y, bnd)
        r["band_pbar_fi_iss"] = score(kind, ["pbar", "fi", "iss_A"], y, bnd)
        out.append(r)
        print(f"  {y}:")
        print(f"    full  pbar={r['full_pbar']:.4f}  +Iss={r['full_pbar_iss']:.4f}"
              f"  delta={r['full_delta']:+.4f}")
        print(f"    band  pbar={r['band_pbar']:.4f}  +Iss={r['band_pbar_iss']:.4f}"
              f"  delta={r['band_delta']:+.4f}")
        print(f"    band  +fi ={r['band_pbar_fi']:.4f}  +fi+Iss="
              f"{r['band_pbar_fi_iss']:.4f}")

    pd.DataFrame(out).to_csv(d / "condindep.csv", index=False)
    print(f"wrote {d/'condindep.csv'}")


if __name__ == "__main__":
    if "--analyze" in sys.argv:
        analyze()
    elif "--recovery" in sys.argv:
        recovery_targets()
    else:
        run_sim()
