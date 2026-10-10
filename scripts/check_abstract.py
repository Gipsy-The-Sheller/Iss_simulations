"""Word count of the abstract, as the journal's limit counts it.

LaTeX commands and their arguments are stripped before counting, so
\\citep{xia2003iss} does not count as a word; inline maths counts as one word
per $...$ group, which is how a copy-editor would read it.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

TEX = Path(__file__).resolve().parents[1] / "manuscript.tex"

BODY = re.compile(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", re.S)
MATH = re.compile(r"\$[^$]*\$")
CMD = re.compile(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?")
BRACES = re.compile(r"[{}\\$]")


def main() -> None:
    m = BODY.search(TEX.read_text(encoding="utf-8"))
    if not m:
        sys.exit("no abstract found")
    t = m.group(1)
    t = MATH.sub(" X ", t)          # each formula reads as one word
    t = CMD.sub(" ", t)
    t = BRACES.sub(" ", t)
    words = [w for w in t.split() if any(c.isalnum() for c in w)]
    print(f"abstract: {len(words)} words")
    if len(words) > 300:
        sys.exit(f"{len(words)} > 300")


if __name__ == "__main__":
    main()
