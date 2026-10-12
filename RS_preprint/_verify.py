# -*- coding: utf-8 -*-
"""Check the generated preprint against ../manuscript.tex.

Two checks:

  1. fresh -- manuscript.tex on disk is what _build_rs.py produces from the
     source right now, so nothing has gone stale or been hand-edited.
  2. body  -- undoing the additions the README documents leaves the source
     body character for character, whitespace and comments aside.  This is
     the evidence for the README's "body text is unchanged" claim.

Unlike MPE_sumb/_verify.py this does not re-run a build by importing it: the
additions are imported as data, and build() is a pure function that is never
called for its side effects.
"""
import io
import re
import sys

import _build_rs as b


def norm(t):
    """Comments out, whitespace collapsed."""
    t = re.sub(r'(?m)^[ \t]*%[^\n]*\n?', '', t)
    return re.sub(r'\s+', ' ', t).strip()


def undo(gen):
    """Remove everything the preprint adds, leaving the source text."""
    out = gen
    for addition in (b.HEADER, b.GRAPHICSPATH, b.POSITIONING_TEX,
                     b.AI_DISCLOSURE, b.DECLARATIONS):
        if addition not in out:
            sys.exit('documented addition missing from %s:\n%s'
                     % (b.OUT, addition[:200]))
        out = out.replace(addition, '', 1)
    out = out.replace(b.EMAIL_CLAUSE, '', 1)
    out = out.replace('\\bibliography{../refs/refs_iss_clean}',
                      '\\bibliography{refs/refs_iss_clean}', 1)
    return out


def context(a, s):
    for i, (x, y) in enumerate(zip(a, s)):
        if x != y:
            return i, a[max(0, i - 60):i + 60], s[max(0, i - 60):i + 60]
    i = min(len(a), len(s))
    return i, a[i:i + 60], s[i:i + 60]


def main():
    src = io.open(b.SRC, encoding='utf-8').read()
    gen = io.open(b.OUT, encoding='utf-8').read()

    if norm(b.build(src)) != norm(gen):
        sys.exit('%s is stale or hand-edited -- rerun _build_rs.py' % b.OUT)
    print('fresh: %s is what _build_rs.py produces from %s' % (b.OUT, b.SRC))

    mine = [l for l in (b.HEADER + b.GRAPHICSPATH).split('\n')
            if l.lstrip().startswith('%')]
    stray = [l for l in gen.split('\n')
             if l.lstrip().startswith('%') and l not in mine]
    if stray:
        sys.exit('drafting comment left in %s:\n%s' % (b.OUT, stray[0]))

    a, s = norm(undo(gen)), norm(src)
    if a != s:
        i, ca, cs = context(a, s)
        sys.exit('body differs from %s at normalized character %d:\n'
                 '  generated: ...%s...\n  source   : ...%s...' % (b.SRC, i, ca, cs))
    print('body:  identical to %s once the documented additions are removed'
          % b.SRC)
    print('ok')


if __name__ == '__main__':
    main()
