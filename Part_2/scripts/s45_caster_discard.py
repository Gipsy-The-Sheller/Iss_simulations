"""Does discarding saturated sites help a method that is immune to LBR?

The discriminating experiment.  Duchene et al. (2022) report that excluding
loci with substitution saturation improves inferences from phylogenomic data.
This paper argues that what such a filter removes is not "saturation" but
sites carrying no topological signal, and that the two readings make different
predictions about a method that is not susceptible to long-branch
repulsion.

CASTER (Zhang, Nielsen & Mirarab 2025) is such a method: it provably corrects
long-branch attraction by assigning negative weights to sites with exactly
three distinct bases.  So:

  * if the benefit of filtering is a saturation effect, CASTER should benefit
    too -- removing saturated sites should help any method;
  * if the benefit is LBR repair, CASTER should not benefit, because there is
    nothing for it to repair.

Neither outcome refutes the other paper; they distinguish the two readings,
which is what is currently missing.

Simulation is exact: sites are drawn i.i.d. from the quartet's pattern
distribution, which is computed in closed form, so the alignment is a genuine
sample from the tree and not an approximation to one.  CASTER is then run
blind -- it sees only a FASTA file.
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

OUT = PART2 / "results" / "caster_discard"
WORK = PART2 / "_caster"
CASTER = Path("E:/Codes/ASTER/exe/caster-site.exe")

_spec = importlib.util.spec_from_file_location(
    "s36", str(HERE / "s36_factors.py"))
s36 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s36)

m = s36.m
from issbench import entropy                                        # noqa: E402

H_PAT = np.array([
    entropy.h_from_counts(np.bincount([(i >> 6) & 3, (i >> 4) & 3,
                                       (i >> 2) & 3, i & 3], minlength=4))
    for i in range(256)])
BASES = np.array(list("ACGT"))
TRUE_SPLIT = "{A,B}|{C,D}"


def write_fasta(seqs, path):
    with open(path, "w") as f:
        for i, n in enumerate("ABCD"):
            f.write(f">{n}\n{''.join(seqs[i])}\n")


SPLITS = {"{A,B}|{C,D}": ("A", "B", "C", "D"),
          "{A,C}|{B,D}": ("A", "C", "B", "D"),
          "{A,D}|{B,C}": ("A", "D", "B", "C")}


def recovered_split(newick):
    """Which of the three unrooted splits a rooted display represents.

    A rooted display of a four-taxon tree shows the split {i,j}|{k,l} by
    making either {i,j} or {k,l} a clade, depending on where the root landed.
    So the split is present iff one of the two pairs appears as a literal
    two-leaf clade.  An earlier version tested for the two pairs appearing
    *together*, which fails whenever the root falls so that only one of them
    is a clade -- the common case, and it made CASTER look as though it never
    recovered the tree at all.
    """
    s = re.sub(r":[0-9.eE+-]+", "", newick)      # branch lengths
    s = re.sub(r"\)[0-9.]+", ")", s)             # support values
    s = s.replace(" ", "")
    for name, (i, j, k, l) in SPLITS.items():
        if (f"({i},{j})" in s or f"({j},{i})" in s
                or f"({k},{l})" in s or f"({l},{k})" in s):
            return name
    return None


def run_caster(fa_path, tag):
    """The unrooted split CASTER reports, or None if it produced no tree."""
    lst = WORK / f"{tag}_list.txt"
    lst.write_text(str(fa_path.resolve()) + "\n")
    tre = WORK / f"{tag}.nw"
    if tre.exists():
        tre.unlink()
    subprocess.run([str(CASTER), "-i", str(lst), "-o", str(tre)],
                   capture_output=True, text=True, cwd=str(WORK))
    if not tre.exists():
        return None
    return recovered_split(tre.read_text().strip())


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    if not CASTER.exists():
        sys.exit(f"CASTER not found at {CASTER}")

    shape = (1.0, 1.0, 0.05, 0.05, 0.02)
    L = 200_000
    n_rep = 8
    rng = np.random.default_rng(20261016)

    print("Does discarding high-entropy sites help CASTER?")
    print(f"  tree a=b=1.0, c=d=0.05, e=0.02 (Farris), L={L}, "
          f"{n_rep} replicates per point")
    print(f"  CASTER: {CASTER}\n")
    print(f"  {'scale':>7} {'TL':>7} {'kept f':>7} "
          f"{'full':>8} {'filtered':>9}")

    rows = []
    for s in (0.6, 1.0, 1.6, 2.5):
        b1 = tuple(v * s for v in shape)
        p = m.pattern_batch(b1, m.SPLITS["12|34"])[0]
        keep = H_PAT <= 0.90 + 1e-12
        f_kept = float(p[keep].sum())
        pf = p * keep
        pf = pf / pf.sum()

        res = {}
        for cond, dist in (("full", p), ("filtered", pf)):
            ok = 0
            for rep in range(n_rep):
                idx = rng.choice(256, size=L, p=dist)
                digits = np.stack([(idx >> 6) & 3, (idx >> 4) & 3,
                                   (idx >> 2) & 3, idx & 3])
                seqs = BASES[digits]
                fa = WORK / f"t_{cond}_{s}_{rep}.fa"
                write_fasta(seqs, fa)
                got = run_caster(fa, f"t_{cond}_{s}_{rep}")
                ok += int(got == TRUE_SPLIT)
            res[cond] = ok / n_rep
        rows.append(dict(scale=s, TL=sum(b1), f_kept=f_kept,
                         full=res["full"], filtered=res["filtered"]))
        print(f"  {s:7.2f} {sum(b1):7.3f} {f_kept:7.4f} "
              f"{res['full']:8.3f} {res['filtered']:9.3f}", flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "caster_discard.csv", index=False)

    if len(d):
        delta = (d.filtered - d.full)
        print(f"\n  mean change in CASTER's recovery rate from filtering: "
              f"{delta.mean():+.3f}")
        print("  A positive change means CASTER benefits from removing the")
        print("  sites the index calls saturated, which would make the")
        print("  benefit generic rather than LBR-specific.")


if __name__ == "__main__":
    main()
