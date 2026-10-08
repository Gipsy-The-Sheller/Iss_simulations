"""Cheap sanity checks on manuscript.tex -- a substitute for not having LaTeX here."""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
tex = (ROOT / "manuscript.tex").read_text(encoding="utf-8")
bib = (ROOT / "refs" / "refs_iss_clean.bib").read_text(encoding="utf-8")

BS = chr(92)
BSR = chr(92) * 2   # regex-escaped backslash
keys = set(re.findall(r"^@\w+\{([^,]+),", bib, re.M))

used: set[str] = set()
for m in re.finditer(BSR + "cite[pt]?\{([^}]*)\}", tex):
    used |= {k.strip() for k in m.group(1).split(",")}

print("cited but missing from bib :", sorted(used - keys) or "none")
print("bib entries / cited        :", len(keys), "/", len(used))
print("unused bib entries         :", len(keys - used))

print("brace balance              :", tex.count("{") - tex.count("}"))

body = tex.split(BS + "begin{document}")[1]
body = re.sub(r"\$[^$]*\$", "", body)          # drop inline math
body = re.sub(r"%[^\n]*", "", body)            # drop comments
for ch in ("_", "&", "#"):
    hits = [ln for ln in body.split("\n") if ch in ln]
    print(f"unescaped {ch!r}             : {len(hits)} line(s)")
    for h in hits[:5]:
        print("      ", h.strip()[:110])

labels = re.findall(BSR + "label\{([^}]*)\}", tex)
refs = re.findall(BSR + "ref\{([^}]*)\}", tex)
print("labels                     :", labels)
print("refs                       :", refs)
print("dangling refs              :", sorted(set(refs) - set(labels)) or "none")

envs = re.findall(BSR + "begin\{(\w+\*?)\}", tex)
print("environments used          :", sorted(set(envs)))
