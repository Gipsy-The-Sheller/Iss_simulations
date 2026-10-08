"""AliSim (IQ-TREE >= 2.x, HPC mode included) driver.

AliSim lives inside the iqtree binary -- there is no separate executable:

    iqtree3 --alisim <prefix> -t <tree.nwk> -m <MODEL> --length L \
            --num-alignments n --seed S -af fasta

Seeds are supplied explicitly for every call so that any run is bit-reproducible.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Sequence

DEFAULT_IQTREE = os.environ.get(
    "IQTREE",
    r"E:\Codes\YRTools\plugins\YR-MPE\YR_MPE\bin\iqtree3\iqtree3.exe",
)


class AliSimError(RuntimeError):
    pass


def iqtree_path() -> str:
    p = os.environ.get("IQTREE", DEFAULT_IQTREE)
    if not Path(p).exists():
        found = shutil.which("iqtree3") or shutil.which("iqtree")
        if not found:
            raise AliSimError(f"IQ-TREE not found; set $IQTREE (tried {p!r})")
        return found
    return p


def run_alisim(prefix: str | Path, tree: str | Path, model: str, length: int,
               n_alignments: int = 1, seed: int = 1, seqtype: str = "DNA",
               fmt: str = "fasta", single_output: bool = True,
               threads: int = 1, log: Optional[str | Path] = None,
               extra: Sequence[str] = ()) -> list[Path]:
    """Simulate `n_alignments` alignments of `length` sites and return their paths."""
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        iqtree_path(), "--alisim", str(prefix),
        "-t", str(tree), "-m", model,
        "--length", str(length),
        "--num-alignments", str(n_alignments),
        "--seed", str(seed),
        "--seqtype", seqtype,
        "-af", fmt,
        "-T", str(threads),
    ]
    if single_output:
        cmd.append("--single-output")
    cmd.extend(extra)
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=prefix.parent,
                       encoding="utf-8", errors="replace")
    if log is not None:
        Path(log).write_text((r.stdout or "") + "\n" + (r.stderr or ""),
                             encoding="utf-8")
    if r.returncode != 0:
        raise AliSimError(f"AliSim failed (rc={r.returncode}):\n{r.stdout}\n{r.stderr}")
    ext = {"fasta": ".fa", "phy": ".phy"}.get(fmt, "." + fmt)
    hits = sorted(prefix.parent.glob(prefix.name + "*" + ext))
    if not hits:
        hits = sorted(prefix.parent.glob(prefix.name + "*"))
    return hits


def alisim_help() -> str:
    r = subprocess.run([iqtree_path(), "--alisim", "-h"],
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    return r.stdout + r.stderr
