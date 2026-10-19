"""Tree construction, Newick I/O, and the Felsenstein/Farris transfer path."""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional, Sequence


@dataclass
class Node:
    name: str = ""
    length: float = 0.0
    children: List["Node"] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        return not self.children

    def leaves(self) -> List["Node"]:
        if self.is_leaf:
            return [self]
        out = []
        for c in self.children:
            out.extend(c.leaves())
        return out


def to_newick(node: Node) -> str:
    """Write a rooted Newick string (root gets no branch length). No trailing ';'."""
    if node.is_leaf:
        return f"{node.name}:{node.length:.10g}"
    inner = ",".join(to_newick(c) for c in node.children)
    return f"({inner}){':' + format(node.length, '.10g') if node.length else ''}"


def write_newick(node: Node, path) -> "Path":
    """Write a Newick file ready for AliSim / IQ-TREE (trailing ';' required)."""
    from pathlib import Path
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(to_newick(node) + ";\n", encoding="utf-8")
    return path


def quartet(a: float, b: float, c: float, d: float, r: float,
            labels: Sequence[str] = ("A", "B", "C", "D")) -> Node:
    """Unrooted quartet ((A,B),(C,D)) with pendants a,b,c,d and internal r.

    Rooted at the midpoint of the internal branch so that to_newick() gives a
    valid rooted representation of the unrooted tree.
    """
    A = Node(labels[0], a)
    B = Node(labels[1], b)
    C = Node(labels[2], c)
    D = Node(labels[3], d)
    n1 = Node("", r / 2.0, [A, B])
    n2 = Node("", r / 2.0, [C, D])
    return Node("", 0.0, [n1, n2])


def star(n: int, d: float, labels: Optional[Sequence[str]] = None) -> Node:
    """Star tree: n leaves, every root-to-leaf branch of length d."""
    labels = labels or [f"T{i+1}" for i in range(n)]
    return Node("", 0.0, [Node(labels[i], d) for i in range(n)])


def balanced(n: int, a: float, r: float,
             labels: Optional[Sequence[str]] = None) -> Node:
    """Balanced binary tree of n tips: pendants a, every internal edge r/2.

    The same construction as scripts/s19_scope.py, which keeps its own copy
    so that its published numbers do not move; this is the shared home for
    the n-taxon experiments.
    """
    labels = labels or [f"T{i+1}" for i in range(n)]
    if n == 1:
        return Node(labels[0], a)
    half = n // 2
    return Node("", r / 2.0,
                [balanced(half, a, r, labels[:half]),
                 balanced(n - half, a, r, labels[half:])])


def scale(node: Node, s: float) -> Node:
    """Every edge length multiplied by s -- the scaling family of a tree."""
    return Node(node.name, node.length * s,
                [scale(c, s) for c in node.children])


def total_length(node: Node) -> float:
    """Sum of all edge lengths below `node` (the root's own length is 0)."""
    return sum(c.length + total_length(c) for c in node.children)


def transfer_path(L: float, S: float, r: float, t: float) -> Node:
    """Felsenstein -> Farris transfer, exactly as specified for experiment 3.

        (a, b, c, d) = (L - t, S, L, S + t),   t in [0, L - S]

    t = 0        -> (L, S, L, S): long branches A, C separated by r  (Felsenstein)
    t = L - S    -> (S, S, L, L): long branches C, D adjacent        (Farris)
    """
    if not 0.0 <= t <= L - S + 1e-12:
        raise ValueError("t must lie in [0, L - S]")
    return quartet(L - t, S, L, S + t, r)


def pairwise_distances(node: Node) -> np.ndarray:
    """True patristic distance matrix between leaves.

    Uses the most recent common ancestor, so leaves in the same clade do NOT
    get their distance inflated by the path through the root.  (An earlier
    version summed root-to-tip depths, which is wrong within a clade.)
    """
    leaves = node.leaves()
    idx = {id(l): i for i, l in enumerate(leaves)}
    n = len(leaves)
    D = np.zeros((n, n))
    path: dict[int, list] = {}          # leaf index -> [(node id, depth from root)]

    def walk(v: Node, anc: list, depth: float):
        here = anc + [(id(v), depth)]
        if v.is_leaf:
            path[idx[id(v)]] = here
            return
        for c in v.children:
            walk(c, here, depth + c.length)

    walk(node, [], 0.0)

    for i in range(n):
        for j in range(i + 1, n):
            pi_, pj = path[i], path[j]
            common = 0
            for (ai, _), (aj, _) in zip(pi_, pj):
                if ai != aj:
                    break
                common += 1
            lca = pi_[common - 1][1] if common else 0.0
            D[i, j] = D[j, i] = (pi_[-1][1] - lca) + (pj[-1][1] - lca)
    return D
