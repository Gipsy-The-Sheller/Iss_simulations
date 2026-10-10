"""Fetch the GigaDB concatenated alignment, in resumable parallel chunks.

The file is 4.6 GB and the link to the S3 endpoint runs at ~47 KB/s per
connection, so a single stream would take overnight.  Range requests are
honoured, so the file is split into fixed-size chunks fetched concurrently and
concatenated at the end; each chunk is verified against its expected length
before being kept, so an interrupted run resumes rather than restarting.

    python scripts/s09_fetch_gigadb.py            # download (resumes)
    python scripts/s09_fetch_gigadb.py --verify    # check the assembled file
"""
from __future__ import annotations

import argparse
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "work" / "empirical" / "ALL.aln.gz"
CHUNKS = ROOT / "work" / "empirical" / "chunks"
URL = ("https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/"
       "10.5524/101001_102000/101041/Concatenated_alignments/"
       "examl-concatenations/WGT.ExaML/ALL.aln.gz")
SIZE = 4_598_391_425
CHUNK = 4 << 20           # 4 MiB: 64 MiB chunks took ~38 min each at
                          # the observed ~900 KB/s aggregate, which is
                          # far too coarse to report or resume
WORKERS = 12         # 32 concurrent connections got the endpoint to
                     # refuse new ones after ~800 chunks
RETRIES = 8


def head_size() -> int:
    req = urllib.request.Request(URL, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers["Content-Length"])


def fetch(args) -> tuple[int, int]:
    idx, start, end = args
    dest = CHUNKS / f"{idx:05d}"
    want = end - start + 1
    if dest.exists() and dest.stat().st_size == want:
        return idx, want                       # already have it
    last = None
    for attempt in range(RETRIES):
        try:
            req = urllib.request.Request(
                URL, headers={"Range": f"bytes={start}-{end}"})
            with urllib.request.urlopen(req, timeout=180) as r:
                data = r.read()
            if len(data) != want:
                raise IOError(f"short read {len(data)} != {want}")
            tmp = dest.with_suffix(".part")
            tmp.write_bytes(data)
            tmp.replace(dest)
            return idx, want
        # catch broadly: http.client.IncompleteRead is an HTTPException, not
        # an OSError, so a narrow catch let one truncated read kill the run
        # at 1074/1097 chunks
        except Exception as e:                                        # noqa: BLE001
            last = e
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"chunk {idx} failed after {RETRIES} attempts: {last}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--workers", type=int, default=WORKERS)
    a = ap.parse_args()

    if a.verify:
        if not OUT.exists():
            sys.exit(f"{OUT} missing")
        got = OUT.stat().st_size
        print(f"{OUT}: {got:,} bytes (expected {SIZE:,}) "
              f"-> {'OK' if got == SIZE else 'MISMATCH'}")
        return

    CHUNKS.mkdir(parents=True, exist_ok=True)
    n = (SIZE + CHUNK - 1) // CHUNK
    jobs = [(i, i * CHUNK, min((i + 1) * CHUNK - 1, SIZE - 1))
            for i in range(n)]
    have = sum(1 for i, _, _ in jobs if (CHUNKS / f"{i:05d}").exists())
    print(f"{SIZE:,} bytes in {n} chunks of {CHUNK >> 20} MiB, "
          f"{have} already present, {a.workers} workers", flush=True)

    t0 = time.time()
    done = 0
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(fetch, j): j for j in jobs}
        for f in as_completed(futs):
            f.result()
            done += 1
            if done % 20 == 0 or done == n:
                el = time.time() - t0
                rate = done * CHUNK / el / 1024 if el else 0
                print(f"  {done}/{n} chunks  {rate:,.0f} KB/s", flush=True)

    print("assembling ...")
    with OUT.open("wb") as fh:
        for i in range(n):
            fh.write((CHUNKS / f"{i:05d}").read_bytes())
    got = OUT.stat().st_size
    print(f"{OUT}: {got:,} bytes "
          f"({'OK' if got == SIZE else 'MISMATCH!'}) in "
          f"{(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
