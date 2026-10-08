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
    """Patristic distance matrix between leaves."""
    leaves = node.leaves()
    idx = {id(l): i for i, l in enumerate(leaves)}
    n = len(leaves)
    D = np.zeros((n, n))

    def walk(v: Node, depth: float):
        if v.is_leaf:
            return [(idx[id(v)], depth)]
        acc = []
        for c in v.children:
            acc.extend(walk(c, depth + c.length))
        return acc

    tips = walk(node, node.length)
    for i, di in tips:
        for j, dj in tips:
            D[i, j] = di + dj
    return D
