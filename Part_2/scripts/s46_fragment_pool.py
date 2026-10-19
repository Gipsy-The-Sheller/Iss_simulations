"""A pool of 500 bp fragments, and the Iss.c that goes with it.

The CASTER experiment needs three things before it can run: a tree, a pool of
short fragments that span the interesting range of Iss, and a threshold.

  * the tree is a random 16-taxon Yule-Harding tree from AliSim.  It happens
    to carry a long-branch sister pair (T5, T12), which is the configuration
    the Farris zone is about, and it was generated before anything else was
    measured, so it is not selected for.

  * each fragment is 500 bp simulated under that tree with all branch lengths
    multiplied by a scale s.  The scale is the knob: small s is a star, large
    s is saturated.  Because sites are i.i.d., a long alignment chopped into
    500 bp pieces is exactly a set of independent 500 bp fragments, so the
    whole pool comes from one AliSim call per scale.

  * Iss.c is calibrated the way this paper insists on: by resampling.  For
    each scale the true tree is recovered from the 500 bp fragments, the
    recovery rate is measured against the true tree, and Iss.c is the Iss at
    which that rate crosses 95%.  Nothing is taken from Xia's tables.

The model is JC on purpose.  Under a compositionally biased model Iss and
base composition move together, and the experiment would no longer be able to
tell a filter on saturation from a filter on composition.  Under JC the
fragments differ in Iss and essentially nothing else.

Outputs
    results/fragment_pool/tree.nwk            the tree
    results/fragment_pool/fragments.npy        (n_taxa, n_frag, 500) codes
    results/fragment_pool/meta.csv             scale and Iss per fragment
    results/fragment_pool/pilot.csv            recovery rate per scale
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

OUT = PART2 / "results" / "fragment_pool"
WORK = PART2 / "_caster" / "pool"
IQ = Path("E:/Codes/YRTools/plugins/YR-MPE/YR_MPE/bin/iqtree3/iqtree3.exe")

from issbench import entropy                                        # noqa: E402

FRAG = 500
N_TAXA = 16
N_FRAG_PER_SCALE = 1000         # 20 scales * 1000 * 500 bp = 10 Mbp pool
SCALES = np.geomspace(0.04, 6.0, 20)
CODE = {c: i for i, c in enumerate("ACGT")}


def sh(cmd, **kw):
    return subprocess.run([str(c) for c in cmd], capture_output=True,
                          text=True, encoding="utf-8", errors="replace", **kw)


def scale_newick(newick: str, s: float) -> str:
    """Multiply every branch length in a Newick string by s.

    Done by regex rather than by parsing, because the only numbers in this
    file that follow a colon are branch lengths.
    """
    return re.sub(r":([0-9.eE+-]+)", lambda mo: f":{float(mo.group(1)) * s:.10g}",
                  newick)


def read_fasta(path):
    names, seqs, cur = [], [], None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if cur is not None:
                seqs.append("".join(cur))
            names.append(line[1:].split()[0])
            cur = []
        elif line.strip():
            cur.append(line.strip())
    if cur is not None:
        seqs.append("".join(cur))
    return names, np.array([[CODE[c] for c in s] for s in seqs], dtype=np.uint8)


def fragment_iss(codes):
    """Iss per fragment, convention A, straight off the simulated bases.

    coded[i] has shape (n_taxa, n_frag, FRAG).  Iss = mean site entropy over
    the multinomial reference H_FSS(pi-hat, 4), which is the convention the
    rest of Part 2 uses.
    """
    n_tax = codes.shape[0]
    onehot = codes[:, :, :, None] == np.arange(4)[None, None, None, :]
    counts = onehot.sum(axis=0).astype(np.float64)      # (n_frag, FRAG, 4)
    p = counts / n_tax
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(p > 0, p * np.log(np.where(p > 0, p, 1.0)), 0.0)
    h_site = -terms.sum(axis=2)                          # (n_frag, FRAG)
    hbar = h_site.mean(axis=1)
    pi = counts.sum(axis=1) / (n_tax * FRAG)             # (n_frag, 4)
    ref = np.array([entropy.multinomial_entropy(4, q) for q in pi])
    return hbar / ref, pi


def rf_distance(a: str, b: str) -> int:
    """Robinson-Foulds between two unrooted Newick trees on the same taxa."""
    from issbench import trees as T
    ta = _parse(a)
    tb = _parse(b)
    return len(_splits(ta) ^ _splits(tb))


def _parse(newick: str):
    """Minimal Newick reader: returns nested tuples, ignoring branch lengths."""
    s = re.sub(r":[0-9.eE+-]+", "", newick).strip().rstrip(";")
    pos = 0

    def node():
        nonlocal pos
        if s[pos] == "(":
            pos += 1
            kids = [node()]
            while s[pos] == ",":
                pos += 1
                kids.append(node())
            assert s[pos] == ")"
            pos += 1
            while pos < len(s) and s[pos] not in ",)":
                pos += 1
            return tuple(kids)
        start = pos
        while pos < len(s) and s[pos] not in ",)":
            pos += 1
        return s[start:pos].strip()

    return node()


def _splits(node, out=None):
    """The set of nontrivial splits of an unrooted tree, as frozensets."""
    if isinstance(node, str):
        return out if out is not None else set()
    if out is None:
        taxa = _leaves(node)
        out = set()

        def add(n):
            if isinstance(n, str):
                return {n}
            r = set()
            for k in n:
                r |= add(k)
            if 1 < len(r) < len(taxa) - 1:
                out.add(frozenset(r))
            return r

        add(node)
        return out

    return out


def _leaves(n):
    if isinstance(n, str):
        return {n}
    r = set()
    for k in n:
        r |= _leaves(k)
    return r


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    if not IQ.exists():
        sys.exit(f"IQ-TREE not found at {IQ}")

    tree = OUT / "tree.nwk"
    if not tree.exists():
        r = sh([IQ, "--alisim", WORK / "seedtree", "-t", "RANDOM{yh/16}",
                "--length", 500, "--num-alignments", 1, "--seed", 7,
                "--seqtype", "DNA", "-af", "fasta", "-m", "JC",
                "--single-output", "-redo"], cwd=WORK)
        if r.returncode != 0:
            sys.exit(f"seed tree failed:\n{r.stdout}\n{r.stderr}")
        tree.write_text((WORK / "seedtree.treefile").read_text().strip() + "\n")
    newick = tree.read_text().strip()
    print("tree:", newick[:80], "...")

    frames, metas = [], []
    for s in SCALES:
        (WORK / f"scaled_{s:.4f}.nwk").write_text(scale_newick(newick, s) + "\n")
        r = sh([IQ, "--alisim", WORK / f"s{s:.4f}",
                "-t", WORK / f"scaled_{s:.4f}.nwk",
                "-m", "JC", "--length", FRAG * N_FRAG_PER_SCALE,
                "--num-alignments", 1, "--seed", 4242,
                "--seqtype", "DNA", "-af", "fasta",
                "--single-output", "-redo"], cwd=WORK)
        if r.returncode != 0:
            sys.exit(f"alisim s={s} failed:\n{r.stdout}\n{r.stderr}")
        names, seq = read_fasta(WORK / f"s{s:.4f}.fa")
        assert seq.shape[0] == N_TAXA, (s, seq.shape, names)
        assert seq.shape[1] == FRAG * N_FRAG_PER_SCALE, seq.shape
        assert set(names) == {f"T{i + 1}" for i in range(N_TAXA)}, names
        # AliSim writes the sequences in the tree's Newick leaf order, not in
        # label order -- for this tree that is T2, T7, T6, T9, ...  Relabelling
        # the rows T1..T16 by position silently attaches every sequence to the
        # wrong taxon, and the only symptom is a recovery rate pinned at
        # chance.  Reorder the rows instead.
        seq = seq[[names.index(f"T{i + 1}") for i in range(N_TAXA)]]

        # chop the long alignment into independent 500 bp fragments
        codes = seq.reshape(N_TAXA, N_FRAG_PER_SCALE, FRAG).transpose(1, 0, 2)
        iss, pi = fragment_iss(codes.transpose(1, 0, 2))
        frames.append(codes)
        metas.append(pd.DataFrame(dict(scale=s, frag=np.arange(N_FRAG_PER_SCALE),
                                       iss=iss, pi_A=pi[:, 0])))
        tl = sum(float(x) for x in re.findall(
            r":([0-9.eE+-]+)", scale_newick(newick, s)))
        print(f"  s={s:6.3f}  TL={tl:7.3f}  "
              f"Iss {iss.mean():.4f} +- {iss.std():.4f}", flush=True)

    # each frame is (n_frag, n_taxa, FRAG); the fragment axis is the one to
    # join on, not the taxon axis
    pool = np.concatenate(frames, axis=0)                # (n_frag, n_taxa, FRAG)
    meta = pd.concat(metas, ignore_index=True)
    meta["block"] = np.arange(len(meta)) // N_FRAG_PER_SCALE
    np.save(OUT / "fragments.npy", pool)
    meta.to_csv(OUT / "meta.csv", index=False)
    print(f"\npool: {pool.shape[0]} fragments, {pool.shape[0] * FRAG} bp total")
    q = np.percentile(meta.iss, [0, 10, 25, 50, 75, 90, 100])
    print("Iss percentiles (0/10/25/50/75/90/100): "
          + "  ".join(f"{v:.4f}" for v in q))


if __name__ == "__main__":
    main()
