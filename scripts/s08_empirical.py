"""Empirical check: where does real data sit on the index curve?

The open question left by the theory is whether the rate-heterogeneity
miscalibration and the resolution collapse bite in practice, and that depends
on where real alignments place the critical point on the index curve.  The
data set is the one Xia & Lemey (2009) use to introduce the index:
eight vertebrate mitochondrial COI sequences, 1512 nt, no gaps.

The design is the one the index is actually applied in.  Third codon
positions of a protein-coding mitochondrial gene are the textbook saturated
partition and first+second are the textbook conserved one, and the reference
topology ships with the data, so each partition can be scored on both the
index and whether it recovers the tree.

Per partition we record the index under both conventions, the mean pairwise
p-distance, the pooled composition (so the two conventions can be related by
the pointwise factor c of Eq. 3), the phylogenetic term Psi, and the
resolution-collapse sensitivity |d ln Iss / d ln Psi| of Eq. 15.  Recovery is
delegated to IQ-TREE, and the inferred topology is compared with the
reference by the symmetric difference of their splits.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from issbench.entropy import h, multinomial_entropy                # noqa: E402
from issbench.infer import iqtree_ml                               # noqa: E402
from issbench.stats import (h_bar, mean_p_distance, pi_hat,        # noqa: E402
                            parsimony_informative_fraction, read_fasta,
                            site_counts)

DATA = ROOT / "data" / "VertCOI.FAS"
REF = ROOT / "data" / "VertCOI.reference.nwk"
OUT = ROOT / "results" / "empirical"

# Xia & Lemey apply the index to the whole alignment and to codon partitions;
# we add the two combinations that separate the codon-position effect from the
# length effect (1st+2nd is the same 1008 sites as 1st,2nd pooled).
PARTITIONS = {
    "all":        (0, 1, 2),
    "pos1":       (0,),
    "pos2":       (1,),
    "pos3":       (2,),
    "pos1+2":     (0, 1),
}


def read_reference() -> str:
    return REF.read_text(encoding="utf-8").strip().rstrip(";")


def splits_of(nwk: str) -> set[frozenset[str]]:
    """Unrooted splits of a Newick string, canonicalised by complement.

    Recursive descent rather than a token loop: an earlier version collapsed
    each clade to a single representative leaf, so every nested clade came out
    incomplete and the shared-split counts were meaningless.
    """
    t = nwk.strip().rstrip(";")
    clades: list[frozenset[str]] = []

    def parse(i: int) -> tuple[frozenset[str], int]:
        if t[i] == "(":
            i += 1
            leaves: set[str] = set()
            while True:
                sub, i = parse(i)
                leaves |= sub
                if i < len(t) and t[i] == ",":
                    i += 1
                    continue
                break
            assert t[i] == ")", (t[i - 20:i + 20], i)
            i += 1
            while i < len(t) and t[i] not in ",()":
                i += 1                       # skip label / support / length
            fs = frozenset(leaves)
            clades.append(fs)
            return fs, i
        j = i
        while j < len(t) and t[j] not in ",()":
            j += 1
        return frozenset([t[i:j].split(":")[0].strip()]), j

    all_leaves, _ = parse(0)
    n = len(all_leaves)
    out = set()
    for c in clades:
        if 1 < len(c) < n - 1:               # non-trivial, and its own complement
            other = all_leaves - c
            out.add(c if min(c) < min(other) else frozenset(other))
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    names, M = read_fasta(DATA)
    N, L = M.shape
    print(f"{DATA.name}: {N} sequences, {L} sites ({L // 3} codons)")
    print("taxa:", ", ".join(names))

    ref = read_reference()
    ref_splits = splits_of(ref)
    print(f"reference tree has {len(ref_splits)} non-trivial splits")

    rows = []
    for tag, positions in PARTITIONS.items():
        sub = np.concatenate([M[:, p::3] for p in positions], axis=1)
        q = pi_hat(sub)
        hb = h_bar(sub)
        c = h(q) / multinomial_entropy(N, q)
        # Psi = <tr P> - 1 is not observable without a tree; report instead the
        # observable quantity the sensitivity formula needs, the mean p-distance
        # and the index, and flag the sensitivity separately (see below).
        pbar = mean_p_distance(sub)
        counts = site_counts(sub)
        p_inf = parsimony_informative_fraction(sub)
        rows.append(dict(partition=tag, n_sites=sub.shape[1], N=N,
                         pbar=pbar, hbar=hb, pi_A=q[0], pi_C=q[1], pi_G=q[2],
                         pi_T=q[3], c=c,
                         iss_A=hb / multinomial_entropy(N, q),
                         iss_B=hb / h(q),
                         n_informative=int(round(p_inf * sub.shape[1])),
                         n_constant=int((counts.max(axis=1) == N).sum()),
                         p_informative=p_inf))
        # write the partition for IQ-TREE
        fa = OUT / f"{tag}.fa"
        with fa.open("w", encoding="utf-8") as fh:
            for nm, seq in zip(names, sub):
                fh.write(f">{nm}\n" + "".join("ACGT"[x] for x in seq) + "\n")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "partitions.csv", index=False)
    print()
    print(df[["partition", "n_sites", "pbar", "iss_A", "iss_B", "c"]]
          .to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    # recovery, by IQ-TREE, against the reference topology
    rec = []
    for tag in PARTITIONS:
        tree = iqtree_ml(OUT / f"{tag}.fa", "GTR+G", seed=1, workdir=OUT)
        if tree is None:
            print(f"  {tag}: IQ-TREE failed")
            continue
        got = splits_of(tree)
        common = ref_splits & got
        rec.append(dict(partition=tag,
                        n_ref=len(ref_splits), n_got=len(got),
                        n_shared=len(common),
                        rf=len(ref_splits - got) + len(got - ref_splits)))
        (OUT / f"{tag}.ml.nwk").write_text(tree.strip() + "\n", encoding="utf-8")
    rdf = pd.DataFrame(rec)
    rdf.to_csv(OUT / "recovery.csv", index=False)
    print()
    print(rdf.to_string(index=False))


if __name__ == "__main__":
    main()
