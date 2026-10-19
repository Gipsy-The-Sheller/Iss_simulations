"""The same band matrices, given to IQ-TREE as a real analysis.

The ML column of s48 is a forced-JC, -fast search.  That is a fine way to
read a likelihood off a known simulation, and it is what the CASTER
comparison needs, but it is not what a phylogenomicist would run, and on
these matrices it is misspecified in a way that matters.

The reason is the construction of the pool.  A band matrix is a concatenation
of 500 bp fragments simulated at scales spanning a factor of a hundred or so,
so the alignment is strongly *rate heterogeneous between fragments* while
each fragment is internally homogeneous.  A single-rate JC fit has one
substitution rate for the whole thing and must therefore compromise.  The
HIGH band makes this sharpest: its fragments are individually saturated, but
the band also reaches down to Iss 0.67, and a model that can put most of its
rate on the informative fragments is in a very different position from one
that cannot.

So this runs IQ-TREE properly on the identical matrices:

    MFP        model selection by ModelFinder, including rate heterogeneity
    full       no -fast; the heuristic search is what might be understating
               maximum likelihood
    -B 1000    ultrafast bootstrap, so the result comes with support values

The matrices are rebuilt from the same permutations as s48, so the comparison
is against exactly the same alignments, and the -fast JC numbers are read
back from s48's cache rather than recomputed.

Outputs
    results/caster_bands/iqtree_bands.csv
"""
from __future__ import annotations

import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PART2 = HERE.parent
ROOT = PART2.parent
sys.path.insert(0, str(ROOT / "src"))

_spec = importlib.util.spec_from_file_location("s48", str(HERE / "s48_caster_bands.py"))
s48 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s48)

OUT = s48.OUT
WORK = s48.WORK
IQ = s48.IQ
L_GRID = [1_000, 2_000, 10_000, 100_000, 1_000_000]


def best_model(iqtree_file: Path):
    if not iqtree_file.exists():
        return None
    m = re.search(r"Best-fit model(?: according to BIC)?: (\S+)",
                  iqtree_file.read_text(errors="replace"))
    return m.group(1) if m else None


def main() -> None:
    s48.TRUTH = (s48.POOL / "tree.nwk").read_text().strip()
    pool = np.load(s48.POOL / "fragments.npy", mmap_mode="r")
    meta = pd.read_csv(s48.POOL / "meta.csv")
    c500, c10k = s48.thresholds()
    iss = meta.iss.values

    raw = {"CLEAN": iss < c500,
           "MID": (iss >= c500) & (iss < c10k),
           "HIGH": iss >= c10k}
    idx_of = {k: np.flatnonzero(v) for k, v in raw.items()}
    conds = {                                   # identical to s48
        "CLEAN": s48.perm(idx_of["CLEAN"], 101),
        "MID": s48.perm(idx_of["MID"], 102),
        "HIGH": s48.perm(idx_of["HIGH"], 103),
        "ALL": s48.perm(np.arange(len(iss)), 7),
    }
    reach = min(len(v) for v in idx_of.values()) * s48.FRAG

    print("IQ-TREE as a real analysis on the same band matrices")
    print("  ModelFinder + rate heterogeneity, full search, 1000 UFBoot\n")
    print(f"  {'condition':>8} {'L':>9} {'ML -fast':>9} {'ML proper':>10} "
          f"{'model':>12} {'mean support':>13}")

    rows = []
    for name, draw in conds.items():
        for L in L_GRID:
            if L > reach:
                continue
            k = L // s48.FRAG
            arr = s48.matrix(pool, draw[:k])
            tag = f"IQF_{name}_{L}"
            fa = WORK / f"{tag}.fa"
            with open(fa, "w") as f:
                for t in range(s48.N_TAXA):
                    f.write(f">T{t + 1}\n"
                            f"{''.join(s48.BASES[c] for c in arr[t])}\n")
            pre = WORK / tag
            for ext in (".treefile", ".iqtree", ".log", ".ckp.gz", ".contree"):
                p = Path(str(pre) + ext)
                if p.exists():
                    p.unlink()
            subprocess.run([str(IQ), "-s", str(fa), "-m", "MFP", "-B", "1000",
                            "-seed", "1", "-redo", "-pre", str(pre),
                            "-nt", "4", "-quiet"],
                           capture_output=True, text=True, cwd=WORK)
            tf = Path(str(pre) + ".treefile")
            rf = (s48.s46.rf_distance(tf.read_text().strip(), s48.TRUTH)
                  if tf.exists() else None)
            model = best_model(Path(str(pre) + ".iqtree"))

            # mean support over the internal branches of the recovered tree
            sup = None
            ct = Path(str(pre) + ".contree")
            if ct.exists():
                vals = [float(x) for x in re.findall(r"\)([0-9.]+):", ct.read_text())]
                sup = float(np.mean(vals)) if vals else None

            fast = None
            c = OUT / "cache" / f"{name}_ML_{L}.rf"
            if c.exists():
                t = c.read_text().strip()
                fast = None if t == "NA" else int(t)
            rows.append(dict(condition=name, L=L, ml_fast=fast, ml_proper=rf,
                             model=model, support=sup, iss_agg=s48.matrix_iss(arr)))
            pd.DataFrame(rows).to_csv(OUT / "iqtree_bands.csv", index=False)
            print(f"  {name:>8} {L:9d} {str(fast):>9} {str(rf):>10} "
                  f"{str(model):>12} "
                  f"{'--' if sup is None else f'{sup:13.1f}'}", flush=True)


if __name__ == "__main__":
    main()
