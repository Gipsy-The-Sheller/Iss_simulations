"""Iss.c as a function of alignment length, measured on the fragment pool.

This is the quantity the whole paper is about, so it is measured the way the
paper insists on: by resampling, on our own tree, with our own recovery rule.
Nothing is read off Xia's tables.

For each tree scale s in the pool and each length L, replicate alignments of
length L are built by concatenating consecutive 500 bp fragments -- which is
valid because the fragments are i.i.d. draws from the same distribution, so
concatenation is just more sampling, not a new simulation.  The true 16-taxon
tree is then recovered by maximum likelihood and the recovery rate is the
fraction of replicates returning it exactly.

Iss.c(L) is the Iss at which that rate crosses 95%.  Because the crossing is
found by interpolation between grid scales and the grid is the same at every
L, the curve is comparable across L by construction.

Two things are reported besides the crossing:

  * the mean RF distance, which keeps discriminating after the exact-recovery
    rate has flattened at zero, and
  * a fit of Iss.c(L) to a power law in L, so the curve can be extrapolated
    past the lengths that were measured -- which is what lets the CASTER
    experiment place a threshold at a length it never simulated at.

Outputs
    results/fragment_pool/issc_curve.csv
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

OUT = PART2 / "results" / "fragment_pool"
WORK = PART2 / "_caster" / "pool"
IQ = Path("E:/Codes/YRTools/plugins/YR-MPE/YR_MPE/bin/iqtree3/iqtree3.exe")

_spec = importlib.util.spec_from_file_location("s46", str(HERE / "s46_fragment_pool.py"))
s46 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s46)

FRAG = s46.FRAG
N_TAXA = s46.N_TAXA
PER_SCALE = s46.N_FRAG_PER_SCALE
L_GRID = (500, 1000, 2000, 5000, 10000)
MAX_REPS = 50
BASES = "ACGT"


def main() -> None:
    pool = np.load(OUT / "fragments.npy", mmap_mode="r")
    meta = pd.read_csv(OUT / "meta.csv")
    newick = (OUT / "tree.nwk").read_text().strip()
    scales = sorted(meta.scale.unique())
    print(f"pool {pool.shape}, {len(scales)} scales")
    print(f"Iss.c(L) by resampling: 16 taxa, JC, IQ-TREE ML, "
          f"up to {MAX_REPS} replicates\n")

    rows = []
    for si, s in enumerate(scales):
        base = si * PER_SCALE
        iss_s = float(meta.loc[meta.scale == s, "iss"].mean())
        for L in L_GRID:
            k = L // FRAG
            n_rep = min(MAX_REPS, PER_SCALE // k)
            ok, rfs = 0, []
            for r in range(n_rep):
                idx = np.arange(base + r * k, base + r * k + k)
                a = np.asarray(pool[idx]).transpose(1, 0, 2).reshape(N_TAXA, -1)
                fa = WORK / f"c_{si}_{L}_{r}.fa"
                with open(fa, "w") as f:
                    for t in range(N_TAXA):
                        f.write(f">T{t + 1}\n"
                                f"{''.join(BASES[c] for c in a[t])}\n")
                pre = WORK / f"c_{si}_{L}_{r}"
                for ext in (".treefile", ".iqtree", ".log", ".ckp.gz"):
                    p = Path(str(pre) + ext)
                    if p.exists():
                        p.unlink()
                rr = s46.sh([IQ, "-s", fa, "-m", "JC", "-seed", 1, "-quiet",
                             "-fast", "-redo", "-pre", pre, "-nt", 2])
                tf = Path(str(pre) + ".treefile")
                if not tf.exists():
                    continue
                d = s46.rf_distance(tf.read_text().strip(), newick)
                rfs.append(d)
                ok += int(d == 0)
            rec = ok / max(len(rfs), 1)
            rows.append(dict(scale=s, L=L, n_rep=len(rfs), iss=iss_s,
                             recovery=rec, mean_rf=float(np.mean(rfs))))
            print(f"  s={s:6.3f} L={L:6d}  Iss={iss_s:.4f}  "
                  f"rec={rec:5.3f}  meanRF={np.mean(rfs):6.2f}", flush=True)

    d = pd.DataFrame(rows)
    d.to_csv(OUT / "issc_curve.csv", index=False)

    # ------------------------------------------------------------- Iss.c(L)
    print("\nIss.c(L): the Iss at which recovery first crosses 95%")
    print(f"  {'L':>7} {'Iss.c':>8} {'rec at crossing':>16}")
    curve = []
    for L in L_GRID:
        sub = d[d.L == L].sort_values("scale")
        cross = None
        for a, b in zip(sub.itertuples(), sub.iloc[1:].itertuples()):
            if a.recovery < 0.95 <= b.recovery:
                w = (0.95 - a.recovery) / max(b.recovery - a.recovery, 1e-12)
                cross = a.iss + w * (b.iss - a.iss)
                break
        curve.append(dict(L=L, iss_c=cross,
                          max_rec=float(sub.recovery.max())))
        print(f"  {L:7d} "
              f"{'--' if cross is None else f'{cross:8.4f}'}"
              f"{'':>8} {sub.recovery.max():16.3f}")
    c = pd.DataFrame(curve)
    c.to_csv(OUT / "issc_vs_L.csv", index=False)

    got = c.dropna(subset=["iss_c"])
    if len(got) >= 2:
        slope, intercept = np.polyfit(np.log(got.L), np.log(got.iss_c), 1)
        print(f"\n  power-law fit over the measured range: "
              f"Iss.c ~ L^({slope:+.4f})  (r = "
              f"{np.corrcoef(np.log(got.L), np.log(got.iss_c))[0, 1]:.4f})")
        for L in (10000, 100000, 1000000):
            if L > got.L.max():
                print(f"  extrapolated Iss.c({L}) = "
                      f"{np.exp(intercept) * L ** slope:.4f}")
    else:
        print("\n  fewer than two crossings -- the grid does not straddle 95%")
        print("  at enough lengths; widen SCALES in s46.")


if __name__ == "__main__":
    main()
