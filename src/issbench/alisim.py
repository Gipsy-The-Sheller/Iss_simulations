"""AliSim (IQ-TREE >= 2.x, HPC mode included) driver.

AliSim lives inside the iqtree binary -- there is no separate executable:

    iqtree3 --alisim <prefix> -t <tree.nwk> -m <MODEL> --length L \
            --num-alignments n --seed S -af fasta

Seeds are supplied explicitly for every call so that any run is bit-reproducible.
"""
from __future__ import annotations

import os
import re
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


def _own(prefix: Path, ext: str):
    """Regex matching only this cell's outputs: {prefix}{ext}, {prefix}_r{i}{ext}."""
    return re.compile(re.escape(prefix.name) + r"(?:_r\d+)?" + re.escape(ext) + "$")


def clear_stale(prefix: str | Path, ext: str) -> None:
    """Delete this cell's own outputs from a previous attempt.

    The glob at the end of run_alisim is anchored on the prefix, so without
    this a re-run of a cell returns files written by the earlier attempt --
    including the concatenated multi-alignment output, which if scored as a
    single alignment gives the fully saturated statistics of a 600-sequence
    "alignment".  Only exact `{prefix}{ext}` and `{prefix}_r{i}{ext}` are
    touched: a bare `prefix*` glob would also match a different cell whose tag
    merely starts with this one.
    """
    prefix = Path(prefix)
    pat = _own(prefix, ext)
    for f in prefix.parent.iterdir():
        if f.is_file() and pat.match(f.name):
            f.unlink()


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
    ext = {"fasta": ".fa", "phy": ".phy"}.get(fmt, "." + fmt)
    clear_stale(prefix, ext)
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=prefix.parent,
                       encoding="utf-8", errors="replace")
    if log is not None:
        Path(log).write_text((r.stdout or "") + "\n" + (r.stderr or ""),
                             encoding="utf-8")
    if r.returncode != 0:
        raise AliSimError(f"AliSim failed (rc={r.returncode}):\n{r.stdout}\n{r.stderr}")
    pat = _own(prefix, ext)
    return sorted(f for f in prefix.parent.iterdir() if pat.match(f.name))


def alisim_help() -> str:
    r = subprocess.run([iqtree_path(), "--alisim", "-h"],
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    return r.stdout + r.stderr
