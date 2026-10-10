"""Stream the GigaDB concatenated alignment and index it sitewise.

The file is a 4.6 GB gzip of a 51 x 322,150,876 PHYLIP-style alignment, one
unwrapped line per taxon.  Nothing here holds the alignment: reading it once
and keeping a per-site count vector over the four states costs
5 * 322M bytes = 1.6 GB, which is the only large object, and every quantity
the index needs is a function of those counts:

    Hbar        = mean over sites of the entropy of the count vector
    p-distance  = mean over sites of (1 - sum_a (c_a/N)^2)
    informative = sites with two states each occurring at least twice

The count form also makes the mean p-distance exact without ever forming a
pairwise matrix: at a site with counts c, the number of ordered mismatched
pairs is N^2 - sum_a c_a^2.

Sites carrying a gap or an ambiguity code are excluded, which is what the
index is applied to in practice.

Windows are runs of `--window` CONSECUTIVE VALID sites, so each window is
built from the same number of usable sites however sparse the region is.  The
normalising HFSS stays the global one, so a window's index varies only through
its own Hbar and windows are comparable with each other.
"""
from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GZ = ROOT / "work" / "empirical" / "ALL.aln.gz"
OUT = ROOT / "results" / "empirical"

# ASCII -> state; 0..3 = A,C,G,T, 4 = gap or ambiguity, 5 = whitespace.
# The default MUST be 4, not 5: only categories 0..4 are counted, so anything
# unmapped that defaults to 5 is silently dropped instead of excluding its
# site -- which is how a first version reported "gap/ambiguity: 0" on an
# alignment that visibly contains dashes.
_LUT = np.full(256, 4, dtype=np.uint8)
for _ws in (9, 10, 13, 32):
    _LUT[_ws] = 5
for ch, k in (("A", 0), ("C", 1), ("G", 2), ("T", 3)):
    _LUT[ord(ch)] = k
    _LUT[ord(ch.lower())] = k


def lines_of(fh):
    """Yield the header and then one (name, sequence-bytes) pair per taxon."""
    head = fh.readline().decode("ascii", "replace").strip()
    n_tax, n_sites = (int(x) for x in head.split()[:2])
    yield n_tax, n_sites, None, None
    for i in range(n_tax):
        line = fh.readline()
        if not line:
            raise EOFError(f"taxon {i}: stream ended early")
        parts = line.strip().split(None, 1)
        if len(parts) != 2:
            raise ValueError(f"taxon {i}: expected 'name sequence', got "
                             f"{len(parts)} field(s)")
        name, seq = parts
        yield n_tax, n_sites, name.decode("ascii", "replace"), seq


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", type=int, default=1000,
                    help="valid sites per window")
    ap.add_argument("--limit", type=int, default=0,
                    help="stop after this many taxa (debugging only)")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    counts = None
    n_tax = n_sites = 0
    with gzip.open(GZ, "rb") as fh:
        for i, (nt, ns, name, seq) in enumerate(lines_of(fh)):
            if name is None:
                n_tax, n_sites = nt, ns
                counts = np.zeros((5, n_sites), dtype=np.uint8)
                print(f"header: {n_tax} taxa x {n_sites:,} sites", flush=True)
                continue
            if len(seq) != n_sites:
                raise ValueError(f"{name}: {len(seq)} sites, expected {n_sites}")
            cat = _LUT[np.frombuffer(seq, dtype=np.uint8)]
            for k in range(5):
                counts[k] += (cat == k)
            if (i + 1) % 10 == 0 or i + 1 == n_tax:
                print(f"  {i + 1}/{n_tax} taxa  {time.time() - t0:.0f}s",
                      flush=True)
            if a.limit and i + 1 >= a.limit:
                break

    N = i + 1 if a.limit else n_tax
    print(f"read {N} taxa in {time.time() - t0:.0f}s; computing", flush=True)

    valid = counts[4] == 0
    n_valid = int(valid.sum())
    n_other = int((~valid).sum())
    print(f"usable sites: {n_valid:,} ({100 * n_valid / n_sites:.1f}%); "
          f"gap/ambiguity: {n_other:,}")

    # Everything is accumulated in blocks over the valid sites.  Materialising
    # counts[:4][:, valid] as float64 would be 4 x 32M x 8 = 1 GB, and the same
    # again for P; blocking keeps the peak to one block plus the counts array.
    W = a.window
    BLK = W * 4000                      # multiple of W, so blocks align to windows
    at = np.arange(BLK)

    def windows_over(sel: np.ndarray, label: str) -> tuple[pd.DataFrame, dict]:
        """Per-window statistics over a subset of the valid sites."""
        n = len(sel)
        n_win = n // W
        n_use = n_win * W
        wh = np.empty(n_win); wp = np.empty(n_win)
        wpi = np.empty(n_win); wt = np.empty(n_win)
        sum_h = sum_p = sum_pi = 0.0
        q_acc = np.zeros(4, dtype=np.int64)
        for s in range(0, n_use, BLK):
            e = min(s + BLK, n_use)
            c = counts[:4][:, sel[s:e]].astype(np.float32)
            tot = c.sum(axis=0)                # taxa present per valid site
            P = c / tot                        # valid sites have tot >= 1
            h = -(np.where(P > 0, P * np.log(P), 0.0)).sum(axis=0)
            p = 1.0 - (P ** 2).sum(axis=0)
            pi = ((c >= 2).sum(axis=0) >= 2).astype(np.float32)
            sum_h += float(h.sum()); sum_p += float(p.sum())
            sum_pi += float(pi.sum())
            q_acc += c.sum(axis=1).astype(np.int64)
            ws, we = s // W, e // W
            if we > ws:
                k = we - ws
                wh[ws:we] = h.reshape(k, W).mean(axis=1)
                wp[ws:we] = p.reshape(k, W).mean(axis=1)
                wpi[ws:we] = pi.reshape(k, W).mean(axis=1)
                wt[ws:we] = tot.reshape(k, W).mean(axis=1)
            del c, tot, P, h, p, pi

        sys.path.insert(0, str(ROOT / "src"))
        from issbench.entropy import (h as entropy,                # noqa: E402
                                      multinomial_entropy)
        q = q_acc / q_acc.sum()
        hbar, pbar = sum_h / n_use, sum_p / n_use
        hfss_a, hfss_b = multinomial_entropy(N, q), entropy(q)
        df = pd.DataFrame(dict(window=np.arange(n_win), hbar=wh, pbar=wp,
                               pi_frac=wpi, mean_taxa=wt))
        df["iss_A"] = df["hbar"] / hfss_a
        df["iss_B"] = df["hbar"] / hfss_b
        info = dict(subset=label, n_sites=n, window=W, hbar=hbar, pbar=pbar,
                    hfss_A=hfss_a, hfss_B=hfss_b, iss_A=hbar / hfss_a,
                    iss_B=hbar / hfss_b, frac_pi=sum_pi / n_use)
        print(f"\n{label}: {n:,} sites, {n_win:,} windows")
        print(f"  Hbar={hbar:.6f} pbar={pbar:.4f}  Iss(A)={info['iss_A']:.4f}"
              f"  Iss(B)={info['iss_B']:.4f}  composition {np.round(q, 4)}")
        return df, info

    valid_idx = np.flatnonzero(valid).astype(np.int32)

    # The parsimony-informative subset of the complete columns, which is the
    # restriction Duchene et al. recommend.  Same counts array, so no extra
    # pass over the alignment is needed.
    pi_sel = np.zeros(n_valid, dtype=bool)
    for s in range(0, n_valid, BLK):
        e = min(s + BLK, n_valid)
        c = counts[:4][:, valid_idx[s:e]]
        pi_sel[s:e] = (c >= 2).sum(axis=0) >= 2
    print(f"parsimony-informative: {int(pi_sel.sum()):,} of {n_valid:,} "
          f"complete sites ({100 * pi_sel.mean():.1f}%)")

    out_rows = []
    for label, sel in (("all complete columns", valid_idx),
                       ("parsimony-informative", valid_idx[pi_sel])):
        d, info = windows_over(sel, label)
        d.to_csv(OUT / (f"gigadb_windows{'' if label.startswith('all') else '_pi'}"
                        ".csv"), index=False)
        out_rows.append(info)
    for r in out_rows:
        r.update(n_tax=N, n_sites=n_sites, n_valid=n_valid)
    pd.DataFrame(out_rows).to_csv(OUT / "gigadb_global.csv", index=False)
    print(f"\ntotal {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
