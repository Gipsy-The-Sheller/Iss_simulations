"""Move the dense result tables between their plain and archived forms.

Everything under results/ is written as plain CSV by the stage that produces
it, which is what makes a run inspectable while it is happening.  Committing
those files verbatim would put ~100 MB of text in the repository; gzipped they
are ~15 MB, and pandas reads either form, so the archive copy is the one that
gets tracked.

    python scripts/s20_pack_results.py            # report what is there
    python scripts/s20_pack_results.py --gzip     # plain -> .csv.gz
    python scripts/s20_pack_results.py --unpack   # .csv.gz -> plain

Only files above --min-mb are touched: the small summary tables the manuscript
cites by value stay readable in the repository without a decompression step.
The plain file is removed after a successful compress, and vice versa; the two
forms are never both left in place by this script.
"""
from __future__ import annotations

import argparse
import gzip
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def candidates(suffix: str, min_bytes: int) -> list[Path]:
    out = []
    for p in sorted(RESULTS.rglob(f"*{suffix}")):
        if p.is_file() and p.stat().st_size >= min_bytes:
            out.append(p)
    return out


def compress(p: Path) -> None:
    tmp = p.with_suffix(p.suffix + ".gz.part")
    with open(p, "rb") as src, gzip.open(tmp, "wb", compresslevel=6) as dst:
        shutil.copyfileobj(src, dst, 1 << 20)
    tmp.replace(p.with_suffix(p.suffix + ".gz"))
    p.unlink()


def decompress(p: Path) -> None:
    out = p.with_suffix("")                 # drops the .gz
    tmp = out.with_suffix(out.suffix + ".part")
    with gzip.open(p, "rb") as src, open(tmp, "wb") as dst:
        shutil.copyfileobj(src, dst, 1 << 20)
    tmp.replace(out)
    p.unlink()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gzip", action="store_true")
    ap.add_argument("--unpack", action="store_true")
    ap.add_argument("--min-mb", type=float, default=0.5)
    a = ap.parse_args()
    if a.gzip and a.unpack:
        sys.exit("--gzip and --unpack are mutually exclusive")
    floor = int(a.min_mb * 1024 * 1024)

    done = 0
    for p in (candidates(".csv", floor) if a.gzip else []):
        n = p.stat().st_size
        compress(p)
        print(f"  gzip   {n / 1048576:7.2f} -> "
              f"{p.with_suffix(p.suffix + '.gz').stat().st_size / 1048576:6.2f} MB"
              f"   {p.relative_to(ROOT)}")
        done += 1
    for p in (candidates(".csv.gz", 0) if a.unpack else []):
        decompress(p)
        print(f"  unpack {p.relative_to(ROOT)}")
        done += 1

    if not done:
        tot = sum(p.stat().st_size for p in candidates(".csv", 0))
        totg = sum(p.stat().st_size for p in candidates(".csv.gz", 0))
        print(f"plain {tot / 1048576:.1f} MB, archived {totg / 1048576:.1f} MB "
              f"-- pass --gzip or --unpack to move the ones over "
              f"{a.min_mb:g} MB")


if __name__ == "__main__":
    main()
