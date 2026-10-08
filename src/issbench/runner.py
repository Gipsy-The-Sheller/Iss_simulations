"""Shared cell loop: simulate, score, infer, append to CSV, resume.

Every experiment writes one row per replicate to results/<name>/cells.csv with a
stable column set; a re-run skips (cell, replicate) pairs already present, so a
killed job restarts cleanly.
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Iterable

import pandas as pd

from .sim import cell_key, score_cell, simulate

COLS = ["exp", "cell", "rep", "tree", "model", "L", "n_seq", "seed",
        "pbar", "hbar", "iss_A", "iss_B", "pi_A", "pi_C", "pi_G", "pi_T",
        "file"]


def out_path(exp: str, root: Path | None = None) -> Path:
    root = Path(root or Path(__file__).resolve().parents[2] / "results")
    d = root / exp
    d.mkdir(parents=True, exist_ok=True)
    return d / "cells.csv"


def load_done(path: Path) -> set:
    if not Path(path).exists():
        return set()
    try:
        df = pd.read_csv(path, usecols=["cell", "rep"])
    except Exception:
        return set()
    return set(zip(df["cell"].astype(str), df["rep"].astype(int)))


def append_rows(path: Path, rows: list[dict]):
    df = pd.DataFrame(rows)
    for c in COLS:
        if c not in df.columns:
            df[c] = None
    df = df[COLS]
    header = not Path(path).exists()
    df.to_csv(path, mode="a", header=header, index=False)


def run_cells(exp: str, cells: Iterable[dict], workroot: Path, seed0: int = 1,
              threads: int = 1, quiet: bool = False) -> int:
    """cells: iterable of dicts with keys tree_file, model, L, n_reps, tag,
    plus anything else to be recorded verbatim (must be JSON-serialisable)."""
    path = out_path(exp)
    done = load_done(path)
    total = 0
    for c in cells:
        tag = c["tag"]
        n = int(c["n_reps"])
        todo = [r for r in range(n) if (tag, r) not in done]
        if not todo:
            continue
        seed = seed0 + (int(cell_key(**{k: v for k, v in c.items()
                                        if k not in ("n_reps",)}), 16) % 100000)
        wd = Path(workroot) / "work" / exp / tag
        t0 = time.time()
        try:
            files = simulate(Path(c["tree_file"]), c["model"], int(c["L"]),
                             len(todo), seed, wd, tag=tag, threads=threads)
        except Exception as e:                                    # noqa: BLE001
            if not quiet:
                print(f"  [{exp}/{tag}] SIMULATION FAILED: {e}", flush=True)
            continue
        rows = score_cell(files)
        for r, rep in zip(rows, todo):
            r["rep"] = rep
            r["cell"] = tag
            r["exp"] = exp
            r["seed"] = seed
            r["tree"] = c.get("tree_nwk", c["tree_file"])
            r["model"] = c["model"]
            r["L"] = int(c["L"])
        append_rows(path, rows)
        total += len(rows)
        if not quiet:
            print(f"  [{exp}/{tag}] {len(rows)} reps in {time.time()-t0:.1f}s"
                  f"  iss_A~{sum(x['iss_A'] for x in rows)/len(rows):.4f}",
                  flush=True)
    return total


def dump_meta(exp: str, meta: dict, root: Path | None = None):
    d = out_path(exp, root).parent
    (d / "meta.json").write_text(json.dumps(meta, indent=1, default=str),
                                 encoding="utf-8")
