"""Three bands of Iss, three methods, one length axis.

The fragment pool is partitioned by Iss at the threshold that the paper's own
calibration actually produced:

    CLEAN   Iss <  Iss.c(10k)                clean by the 10 kbp standard
    MID     Iss.c(10k) <= Iss < Iss.c(500)   saturated by the 500 bp standard,
                                             not by the 10 kbp one
    HIGH    Iss >= Iss.c(500)                saturated by both

The two cuts come from different places on purpose.  The lower one is ours,
measured by resampling at 10 kbp on the tree the fragments were simulated
from.  The upper one is Xia's, read off Fig. 2 at t = 16 and L = 500: 0.72
for the symmetric tree, 0.51 for the asymmetric, the lower of the two used
here.  Our own 500 bp calibration has no crossing at all on this tree,
because recovery there peaks below 95% at every branch length.  That is a
fact about the tree rather than about 500 bp -- the same t and the same L
give Xia a threshold and give us none -- and it is one more instance of the
paper's point, so it is used rather than worked around.

MID is therefore exactly the band the two calibrations disagree about, which
is the band this experiment exists to price.

Every condition is truncated to exactly the same number of sites, so the only
thing that differs between the columns is *which* fragments are in the
alignment.  That is the whole design: matched on length, differing only in
composition.

Three methods are run on every alignment:

    ML           IQ-TREE under JC, the method Iss.c is defined against
    CASTER-site  site weights, provable under F84 and hence under JC
    CASTER-pair  weights on pairs of sites, provable for the lumpable class

Predictions, stated before the runs:

  * MID is not dead.  Fragments the index calls saturated still reconstruct
    the tree once enough of them are concatenated, because Iss is a per-site
    statistic and does not know how many sites there are.  If so, the filter
    is discarding data that is good at the scale it is actually used at.
  * HIGH is dead, so the index is not useless -- it locates a real tail.
  * And the discriminating one: whatever MID costs ML, it costs CASTER less,
    because CASTER's guarantee does not depend on the thing the filter is
    trying to repair.

The ALL condition is a uniform random draw from the pool.  It must be a
draw and not a prefix, because the pool is stored in scale order and a prefix
would be a low-scale subset, not a sample of the whole.

Results are cached per (condition, method, length, replicate), so a killed
run resumes where it stopped.
"""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PART2 = HERE.parent
ROOT = PART2.parent
sys.path.insert(0, str(ROOT / "src"))

OUT = PART2 / "results" / "caster_bands"
WORK = PART2 / "_caster" / "bands"
POOL = PART2 / "results" / "fragment_pool"
IQ = Path("E:/Codes/YRTools/plugins/YR-MPE/YR_MPE/bin/iqtree3/iqtree3.exe")
CASTER_SITE = Path("E:/Codes/ASTER/exe/caster-site.exe")
CASTER_PAIR = Path("E:/Codes/ASTER/exe/caster-pair.exe")

_spec = importlib.util.spec_from_file_location("s46", str(HERE / "s46_fragment_pool.py"))
s46 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s46)

from issbench import entropy                                        # noqa: E402

FRAG = s46.FRAG
N_TAXA = s46.N_TAXA
BASES = "ACGT"

METHODS = ["ML", "CASTER-site", "CASTER-pair"]
L_GRID = [500, 1_000, 2_000, 5_000, 10_000, 30_000, 100_000, 300_000, 1_000_000]
JACK_LENGTHS = [500, 10_000, 100_000, 1_000_000]
JACK_K = 10
# Iss.c(500) for t = 16, read off Xia et al. (2003) Fig. 2: 0.72 for the
# symmetric tree and 0.51 for the asymmetric one.  The lower of the two is
# used as the primary cut, so "what a 500 bp filter would discard" is the
# larger and more generous band.
XIA_C500_ASYM = 0.51
XIA_C500_SYM = 0.72

TRUTH = ""


def perm(idx, seed):
    """A random order of a condition's fragments.

    Any prefix of it is then a uniform random subset of the condition, so the
    length axis is a nested family of proper samples and the curves at
    different lengths are not confounded with which fragments were drawn.
    Taking the pool in its stored order instead would have made each prefix a
    low-scale subset, because the pool is stored in scale order.
    """
    return np.random.default_rng(seed).permutation(idx)


def matrix(pool, idx):
    return np.asarray(pool[idx]).transpose(1, 0, 2).reshape(N_TAXA, -1)


def matrix_iss(arr):
    """Iss of a whole alignment, convention A.

    This is the number the argument turns on.  Iss is a per-site average, so
    bundling fragments does not change it -- which is exactly why the same
    sequences can be called saturated at 500 bp and unsaturated at 10 kbp.
    Recording it per condition lets that be checked rather than asserted.
    """
    n_tax = arr.shape[0]
    onehot = arr[:, :, None] == np.arange(4)[None, None, :]
    c = onehot.sum(axis=0).astype(float)              # (n_sites, 4)
    p = c / n_tax
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(p > 0, p * np.log(np.where(p > 0, p, 1.0)), 0.0)
    hbar = float(-terms.sum(axis=1).mean())
    pi = c.sum(axis=0) / c.sum()
    return hbar / entropy.multinomial_entropy(4, pi)


def run(tag, method, arr):
    """Cache-and-run one alignment; returns the RF distance or None."""
    cache = OUT / "cache" / f"{tag}.rf"
    if cache.exists():
        txt = cache.read_text().strip()
        return None if txt == "NA" else int(txt)
    cache.parent.mkdir(parents=True, exist_ok=True)

    fa = WORK / f"{tag}.fa"
    with open(fa, "w") as f:
        for t in range(N_TAXA):
            f.write(f">T{t + 1}\n{''.join(BASES[c] for c in arr[t])}\n")

    got = None
    if method == "ML":
        pre = WORK / tag
        for ext in (".treefile", ".iqtree", ".log", ".ckp.gz"):
            p = Path(str(pre) + ext)
            if p.exists():
                p.unlink()
        subprocess.run([str(IQ), "-s", str(fa), "-m", "JC", "-seed", "1",
                        "-quiet", "-fast", "-redo", "-pre", str(pre),
                        "-nt", "2"], capture_output=True, text=True, cwd=WORK)
        tf = Path(str(pre) + ".treefile")
        if tf.exists():
            got = tf.read_text().strip()
    else:
        exe = CASTER_SITE if method == "CASTER-site" else CASTER_PAIR
        lst = WORK / f"{tag}_list.txt"
        lst.write_text(str(fa.resolve()) + "\n")
        tre = WORK / f"{tag}.nw"
        if tre.exists():
            tre.unlink()
        subprocess.run([str(exe), "-i", str(lst), "-o", str(tre)],
                       capture_output=True, text=True, cwd=WORK)
        if tre.exists() and tre.read_text().strip():
            # CASTER's -o file holds one tree: the final one
            got = tre.read_text().strip().splitlines()[-1].strip()

    val = None if got is None else s46.rf_distance(got, TRUTH)
    cache.write_text("NA" if val is None else str(val))
    return val


def thresholds():
    """Both cuts from our own resampling, under Xia's falling-crossing rule.

    s47b produces these; the alternative of reading Iss.c(500) off Xia's
    Fig. 2 is kept in the constants below for reference but is not used,
    because those values are for a different tree and the whole point is that
    the tree matters.  Ours is measured on the tree the fragments came from.

    The endpoints of the curve are well determined; the intermediate lengths
    are not, at 50 replicates.  Only the endpoints are needed here.
    """
    up = pd.read_csv(POOL / "issc_upper.csv")

    def t(L, col="iss_c_uniform"):
        row = up[up.L == L]
        if not len(row):
            return None
        v = row[col].iloc[0]
        return None if v is None or not np.isfinite(v) else float(v)

    return t(500), t(10_000)


def main() -> None:
    global TRUTH
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    (OUT / "cache").mkdir(exist_ok=True)
    TRUTH = (POOL / "tree.nwk").read_text().strip()

    pool = np.load(POOL / "fragments.npy", mmap_mode="r")
    meta = pd.read_csv(POOL / "meta.csv")
    c500, c10k = thresholds()
    if c500 is None or c10k is None:
        sys.exit("cannot place the thresholds -- run s47 first")

    iss = meta.iss.values
    # "saturated at length L" means Iss > Iss.c(L), and Iss.c RISES with L, so
    # the band the two calibrations disagree about lies between the cuts: too
    # saturated for a 500 bp reading, acceptable to a 10 kbp one.
    raw = {"CLEAN": iss < c500,
           "MID": (iss >= c500) & (iss < c10k),
           "HIGH": iss >= c10k}
    idx_of = {k: np.flatnonzero(v) for k, v in raw.items()}
    print(f"Iss.c(500)  = {c500:.4f}   Iss.c(10k) = {c10k:.4f}"
          "   (both measured here)")
    print("  Iss.c rises with SeqLen, so MID -- Iss between the two cuts --"
          "\n  is exactly the data a 500 bp reading rejects and a 10 kbp"
          "\n  reading accepts.")
    print(f"\n  {'band':>8} {'fragments':>10} {'sites':>10} {'mean Iss':>9} "
          f"{'min':>7} {'max':>7}")
    for k, v in raw.items():
        print(f"  {k:>8} {v.sum():10d} {v.sum() * FRAG:10d} "
              f"{iss[v].mean():9.4f} {iss[v].min():7.4f} {iss[v].max():7.4f}")
    print(f"  {'ALL':>8} {len(iss):10d} {len(iss) * FRAG:10d} "
          f"{iss.mean():9.4f} {iss.min():7.4f} {iss.max():7.4f}")

    conds = {
        "CLEAN": perm(idx_of["CLEAN"], 101),
        "MID": perm(idx_of["MID"], 102),
        "HIGH": perm(idx_of["HIGH"], 103),
        "ALL": perm(np.arange(len(iss)), 7),
    }
    # the grid is capped by the smallest band, since every condition must
    # reach every plotted length
    reach = min(len(v) for v in idx_of.values()) * FRAG
    grid = [L for L in L_GRID if L <= reach]
    print(f"\n  matched length axis: {grid}  (smallest band: {reach} sites)")
    print(f"  methods: {', '.join(METHODS)}")

    rows = []
    for name, draw in conds.items():
        for L in grid:
            k = L // FRAG
            arr = matrix(pool, draw[:k])
            rec = {meth: run(f"{name}_{meth}_{L}", meth, arr) for meth in METHODS}
            rows.append(dict(condition=name, L=L, n_frag=k,
                             iss_agg=matrix_iss(arr),
                             **{m.replace("-", "_"): rec[m] for m in METHODS}))
            pd.DataFrame(rows).to_csv(OUT / "bands_curve.csv", index=False)
            print(f"  {name:>8} {L:9d} Iss={rows[-1]['iss_agg']:.4f}  "
                  + "  ".join(
                      f"{m}={rec[m] if rec[m] is not None else -1}"
                      for m in METHODS), flush=True)

    if "--curve-only" in sys.argv:
        return

    # ------------------------------------------------------ block jackknife
    print(f"\nblock jackknife over the fragment order, K={JACK_K}")
    print(f"  {'condition':>8} {'L':>9} {'method':>12} {'full':>5} "
          f"{'lobo mean':>10} {'jack SE':>8}")
    jrows = []
    for name, draw in conds.items():
        for L in [x for x in JACK_LENGTHS if x in grid]:
            k = L // FRAG
            prefix = draw[:k]
            bks = np.array_split(prefix, JACK_K)
            for meth in METHODS:
                full = run(f"{name}_{meth}_{L}", meth, matrix(pool, prefix))
                vals = []
                for b, blk in enumerate(bks):
                    keep = np.setdiff1d(prefix, blk, assume_unique=True)
                    if len(keep) < 2:
                        continue
                    v = run(f"{name}_{meth}_{L}_j{b}", meth, matrix(pool, keep))
                    if v is not None:
                        vals.append(v)
                if len(vals) < 2:
                    continue
                v = np.array(vals, float)
                K = len(v)
                # delete-a-block jackknife: SE of the full-sample estimate
                se = float(np.sqrt((K - 1) / K * ((v - v.mean()) ** 2).sum()))
                jrows.append(dict(condition=name, L=L, method=meth,
                                  full=full, lobo_mean=float(v.mean()),
                                  jack_se=se, lobo_sd=float(v.std(ddof=1)),
                                  n=K, values=list(v)))
                pd.DataFrame(jrows).to_csv(OUT / "bands_jackknife.csv",
                                           index=False)
                print(f"  {name:>8} {L:9d} {meth:>12} "
                      f"{full if full is not None else -1:5d} "
                      f"{v.mean():10.3f} {se:8.3f}", flush=True)


if __name__ == "__main__":
    main()
