# -*- coding: utf-8 -*-
"""Check that the MPE version differs from the original only where intended.

Compares the body of the two .tex files (Introduction ... before the
bibliography) line by line, after applying the British->American map to the
original, so that anything left over is a structural change and must be listed
in the README.
"""
import io, re, sys

sys.path.insert(0, '.')
from _build_mpe import american  # reuse the exact same spelling map

def body(path):
    s = io.open(path, encoding='utf-8').read()
    i = s.index('\\section{Introduction}')
    j = s.index('\\bibliographystyle')
    s = s[i:j]
    s = re.sub(r'(?m)^\s*%.*\n?', '', s)          # drop drafting comments
    s = re.sub(r'\\addcontentsline\{[^}]*\}\{[^}]*\}\{[^}]*\}\n', '', s)
    return s

def lines(t):
    return [re.sub(r'\s+', ' ', l).strip() for l in t.split('\n') if l.strip()]

orig = lines(american(body('../manuscript.tex')))
new = lines(body('manuscript.tex'))

from collections import Counter
co, cn = Counter(orig), Counter(new)
print('original body lines: %d   MPE body lines: %d' % (len(orig), len(new)))
print('\n--- only in the original ---')
for l in sorted((co - cn).elements()):
    print('  ' + l[:110])
print('\n--- only in the MPE version ---')
for l in sorted((cn - co).elements()):
    print('  ' + l[:110])
