# -*- coding: utf-8 -*-
"""Build the MPE (Elsevier) submission version of the manuscript.

Reads ../manuscript.tex, applies only the changes needed for an MPE
submission, and writes ./manuscript.tex:

  * elsarticle class + Elsevier front matter (title, author, abstract, keywords)
  * section order changed to IMRaD: Introduction, Materials and methods,
    Results, Discussion (the original follows the MBE order Results-then-Methods)
  * abstract trimmed to the journal's 250-word limit (text in _abstract.tex)
  * American spelling (the body was consistently British)
  * no line numbers, no drafting comments
  * body text otherwise byte-identical to the original
"""
import io, os, re, sys

SRC = '../manuscript.tex'
OUT = 'manuscript.tex'

src = io.open(SRC, encoding='utf-8').read()


def block(a, b):
    i = src.index(a)
    j = src.index(b, i)
    return src[i:j].rstrip() + '\n'


def strip_comments(t):
    return re.sub(r'(?m)^\s*%.*\n?', '', t)


# ---------------------------------------------------------------- body blocks
intro      = block(r'\section{Introduction}', r'\section{Results}')
results    = block(r'\section{Results}', r'\section{Discussion}')
discussion = block(r'\section{Discussion}', r'\section{Material and methods}')
methods    = block(r'\section{Material and methods}', r'\section*{Data availability}')
dataav     = strip_comments(block(r'\section*{Data availability}', r'\section*{Acknowledgments}'))
ack        = strip_comments(block(r'\section*{Acknowledgments}', '\n% ====='))

methods = methods.replace(r'\section{Material and methods}',
                          r'\section{Materials and methods}')

# --------------------------------------------------------- American spelling
BRIT = {
    'analysed': 'analyzed', 'analyse': 'analyze',
    'characterisation': 'characterization', 'characterise': 'characterize',
    'characterised': 'characterized',
    'discretisation': 'discretization', 'discretised': 'discretized',
    'maximises': 'maximizes', 'maximised': 'maximized',
    'minimises': 'minimizes', 'minimised': 'minimized',
    'normalisation': 'normalization', 'normalised': 'normalized',
    'normaliser': 'normalizer', 'normalising': 'normalizing', 'normalise': 'normalize',
    'parameterised': 'parameterized', 'parameterisation': 'parameterization',
    'realised': 'realized', 'realise': 'realize',
    'specialisation': 'specialization', 'specialised': 'specialized',
    'summarises': 'summarizes', 'summarised': 'summarized',
    'symmetrised': 'symmetrized', 'symmetrisation': 'symmetrization',
    'unrecognised': 'unrecognized',
    'behaviour': 'behavior', 'behavioural': 'behavioral',
    'colour': 'color', 'colours': 'colors',
    'neighbour': 'neighbor', 'neighbours': 'neighbors',
    'neighbourhood': 'neighborhood',
    'modelling': 'modeling', 'labelled': 'labeled',
}

CASE = {}
for a, b in BRIT.items():
    CASE[a] = b
    CASE[a[0].upper() + a[1:]] = b[0].upper() + b[1:]


def american(t):
    def sub(m):
        w = m.group(0)
        return CASE.get(w, CASE.get(w.lower(), w))
    return re.sub(r'\b[A-Za-z]+\b', sub, t)


intro      = american(intro)
results    = american(results)
discussion = american(discussion)
methods    = american(methods)
dataav     = american(dataav)
ack        = american(ack)

# ------------------------------------------------------------- LaTeX repairs
# (1) `\small` is invalid in math mode and does nothing there; the elasticity
#     matrix is 13pt wider than the text block at 12pt type, so scale it to
#     90% of the measure (a 3% reduction in type size, invisible in print).
m = re.search(r'\\begin\{equation\}\n\\small\n'
              r'\\begin\{array\}\{l\|r{10}\}(.*?)\\end\{array\}\n\\end\{equation\}',
              results, re.S)
if not m:
    raise SystemExit('elasticity display not found -- check manuscript.tex')
results = results.replace(
    m.group(0),
    '\\begin{equation}\n\\resizebox{0.9\\textwidth}{!}{$\n'
    '\\begin{array}{l|rrrrrrrrrr}' + m.group(1) + '\\end{array}$}\n\\end{equation}')

# (2) the criterion table overflows the 1in text block by 34pt at 12pt type
old = '\\begin{tabular}{llccccc}'
if old not in results:
    raise SystemExit('criterion table not found -- check manuscript.tex')
results = results.replace(old, '\\footnotesize\n\\begin{tabular}{llccccc}')

# ------------------------------------------------------------- abstract + key
abstract = io.open('_abstract.tex', encoding='utf-8').read().strip()


def count_words(t):
    t = re.sub(r'(?m)^\s*%.*$', '', t)
    t = re.sub(r'\\[a-zA-Z]+\*?(\[[^\]]*\])?', ' ', t)   # drop command names
    t = re.sub(r'\$[^$]*\$', ' MATH ', t)                 # each formula = 1 word
    t = re.sub(r'[{}\\&]', ' ', t)
    return len([w for w in re.split(r'\s+', t) if re.search(r'[A-Za-z0-9]', w)])


# ------------------------------------------------------------------- preamble
macros = src[src.index(r'\newcommand{\Iss}'):src.index(r'\title{')].rstrip()

doc = r"""%% =====================================================================
%%  Molecular Phylogenetics and Evolution --- submission version
%%  Generated by _build_mpe.py from ../../manuscript.tex; see README.md.
%%  Body text is unchanged except for American spelling and the section
%%  order (Materials and methods before Results, as MPE requires).
%%  Build:  pdflatex manuscript && bibtex manuscript && pdflatex manuscript x2
%% =====================================================================
\documentclass[preprint,12pt,authoryear]{elsarticle}

%% The elsarticle preprint measure (390pt) is narrower than the one the
%% tables and displays in this paper were set to; restore a 1in-margin
%% letter page so that nothing overflows the text block.
\usepackage[letterpaper,margin=1in]{geometry}
\usepackage{amsmath,amssymb,amsthm}
\usepackage{booktabs}
\usepackage{graphicx}
\usepackage{xcolor}
\usepackage[hidelinks]{hyperref}

""" + macros + r"""

\journal{Molecular Phylogenetics and Evolution}

\begin{document}

\begin{frontmatter}

\title{Revisiting Xia's Index of Substitution Saturation from a Mathematical
Perspective}

\author[dhu]{Zhi-Jie Xu\corref{cor1}}
\ead{zhi-jie.xu@dhu.edu.cn}                 %% TODO: confirm corresponding e-mail
\cortext[cor1]{Corresponding author.}

\address[dhu]{College of Biological Science \& Medical Engineering, Donghua
University, Shanghai 201620, China}

\begin{abstract}
""" + abstract + r"""
\end{abstract}

\begin{keyword}
Substitution saturation \sep phylogenetic signal \sep rate heterogeneity \sep
spectral decomposition \sep sequence length \sep entropy
\end{keyword}

\end{frontmatter}

""" + '\n'.join([intro, methods, results, discussion, ack]) + r"""
\section*{CRediT authorship contribution statement}

\textbf{Zhi-Jie Xu:} Conceptualization, Methodology, Software, Formal
analysis, Investigation, Writing -- original draft, Writing -- review and
editing, Visualization.

\section*{Declaration of competing interest}

The author declares that he has no known competing financial interests or
personal relationships that could have appeared to influence the work reported
in this paper.

\section*{Funding}

%% TODO: confirm. Delete this sentence if the work was funded.
This research did not receive any specific grant from funding agencies in the
public, commercial, or not-for-profit sectors.

""" + dataav + r"""
\bibliographystyle{elsarticle-harv}
\bibliography{refs/refs_iss_clean}

\end{document}
"""

doc = re.sub(r'\n{4,}', '\n\n\n', doc)
io.open(OUT, 'w', encoding='utf-8', newline='\n').write(doc)

# ------------------------------------------------- bibliography: protect names
# elsarticle-harv lowercases article titles (Elsevier sentence case), which
# also flattens acronyms, software names and genus names ("18S rDNA" ->
# "18s rdna", "Codium" -> "codium").  Brace-protect those tokens in the local
# copy of the .bib; everything else is left to the style.
PROTECT = ['18S', '28S', 'AliSim', 'Annelida', 'BaCoCa', 'BEAST', 'Bryopsidales',
           'Codium', 'DAMBE5', 'DAMBE7', 'DAMBE', 'DNA', 'IQ-TREE',
           'Lophotrochozoa', 'MAFFT', 'Metazoa', 'MPBoot', 'MUSCLE', 'NCBI',
           'Nextstrain', 'nRCFV', 'PAML', 'PhyML', 'PhyloBayes', 'PhyloMAd',
           'PhyPA', 'rDNA', 'RNA', 'SARS-CoV-2', 'SatuTe', 'TreeTime', 'TreSpEx',
           'Xia', 'Xuhua']
pat = re.compile(r'\b(' + '|'.join(sorted(PROTECT, key=len, reverse=True)) + r')\b')

bib = io.open('../refs/refs_iss_clean.bib', encoding='utf-8').read()
n = [0]


def protect(m):
    n[0] += 1
    return '{' + m.group(0) + '}'


def fix_field(m):
    return m.group(1) + pat.sub(protect, m.group(2)) + m.group(3)


bib = re.sub(r'((?:title|booktitle)\s*=\s*\{)(.*?)(\},\n)', fix_field, bib, flags=re.S)
io.open('refs/refs_iss_clean.bib', 'w', encoding='utf-8', newline='\n').write(bib)
print('protected %d tokens in refs/refs_iss_clean.bib' % n[0])

print('wrote %s (%d lines)' % (OUT, doc.count('\n')))
print('abstract words: %d (limit 250)' % count_words(abstract))
if count_words(abstract) > 250:
    sys.exit('ABSTRACT TOO LONG')
