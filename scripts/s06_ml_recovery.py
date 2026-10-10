"""Experiment 6, third pass -- ML recovery of the true quartet by IQ-TREE.

The MPBoot pass (06b) covered maximum parsimony; this one covers maximum
likelihood, the criterion the critical tree length is calibrated against.  One
IQ-TREE call per alignment at ~0.09 s makes the full 108,000 cost ~2.6 h on one
core, so the work is sharded across processes instead:

    python scripts/s06_ml_recovery.py --shard i --of n     # one shard
    python scripts/s06_ml_recovery.py --merge              # combine shards

Shard i handles every n-th alignment of the sorted work list and appends to
`recovery_ml.part{i}.csv`.  Each part carries its own done-set, so an
interrupted shard resumes where it stopped; shards never write to the same
file.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench.infer import TRUE_SPLIT, iqtree_ml, split_of_newick  # noqa: E402

EXP = "exp6_condindep"
MODEL = "JC69"


def part_path(shard: int) -> Path:
    return ROOT / "results" / EXP / f"recovery_ml.part{shard}.csv"


def alignments() -> list[Path]:
    return sorted((ROOT / "work" / EXP).glob("*/*_r*.fa"))


def run_shard(shard: int, of: int) -> None:
    files = alignments()[shard::of]
    out = part_path(shard)
    done: set[tuple[str, int]] = set()
    if out.exists():
        prev = pd.read_csv(out)
        done = set(zip(prev["cell"].astype(str), prev["rep"].astype(int)))
    print(f"[shard {shard}/{of}] {len(files)} alignments, "
          f"{len(done)} already done", flush=True)
    rows = []
    for i, f in enumerate(files):
        cell = f.parent.name
        rep = int(f.stem.rsplit("_r", 1)[1])
        if (cell, rep) in done:
            continue
        # the seed is keyed to the replicate, not to the shard, so the same
        # alignment gets the same IQ-TREE seed whichever shard picks it up
        tree = iqtree_ml(f, MODEL, seed=1 + rep, workdir=f.parent)
        rows.append({"cell": cell, "rep": rep,
                     "ML": int(split_of_newick(tree) == TRUE_SPLIT)})
        if len(rows) >= 200:
            pd.DataFrame(rows).to_csv(out, mode="a",
                                      header=not out.exists(), index=False)
            rows = []
        if (i + 1) % 2000 == 0:
            print(f"[shard {shard}/{of}] {i + 1}/{len(files)}", flush=True)
    if rows:
        pd.DataFrame(rows).to_csv(out, mode="a",
                                  header=not out.exists(), index=False)
    print(f"[shard {shard}/{of}] wrote {out}", flush=True)


def merge(of: int) -> None:
    parts = [part_path(i) for i in range(of)]
    have = [p for p in parts if p.exists()]
    if not have:
        print("no shard files found")
        return
    df = pd.concat([pd.read_csv(p) for p in have], ignore_index=True)
    df = df.dropna(subset=["cell", "rep", "ML"])
    df["rep"] = df["rep"].astype(int)
    df = df.drop_duplicates(["cell", "rep"])
    out = ROOT / "results" / EXP / "recovery_ml.csv"
    df.to_csv(out, index=False)
    n = df.groupby("cell").size()
    print(f"merged {len(have)} parts -> {out}: {len(df)} rows, "
          f"{df.cell.nunique()} cells, reps/cell {n.value_counts().to_dict()}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int)
    ap.add_argument("--of", type=int, default=1)
    ap.add_argument("--merge", action="store_true")
    a = ap.parse_args()
    if a.merge:
        merge(a.of)
    elif a.shard is not None:
        run_shard(a.shard, a.of)
    else:
        ap.error("pass --shard i --of n, or --merge")
