"""Locate a results table, tolerating the archived (gzipped) form.

The dense per-replicate tables are what make the repository large -- the two
`cells.csv` files alone are 60 MB and compress by 7x -- so the archived copy
of a table is kept as `<name>.csv.gz`.  pandas reads and writes gzip
transparently, so nothing but the path changes.

Both forms may be on disk at once, because every stage writes the plain file
and the gzip copy is made afterwards.  Whichever is newer wins, so re-running
a stage picks up its own fresh output instead of the stale archive.
"""
from __future__ import annotations

from pathlib import Path


def res(d, name: str) -> Path:
    """Path to table `name` in directory `d`, gzipped copy preferred if newer."""
    d = Path(d)
    plain, gz = d / name, d / (name + ".gz")
    have = [p for p in (plain, gz) if p.exists()]
    if not have:
        return plain                       # let the caller fail on a real path
    return max(have, key=lambda p: p.stat().st_mtime)
