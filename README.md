# Iss: a spectral audit of the substitution-saturation index

Code, results and manuscript for a methodological critique of the index of
substitution saturation,

> Xia, X., Xie, Z., Salemi, M., Chen, L. & Wang, Y. (2003). An index of
> substitution saturation and its application. *Molecular Phylogenetics and
> Evolution* **26**, 1–7.

$\mathrm{Iss}=\bar H/H_{\mathrm{FSS}}$ is computed from an alignment alone and
is compared against a tabulated critical value $\mathrm{Iss.c}$ to decide
whether a partition is too saturated for phylogenetics. This repository
contains a spectral treatment of the statistic (what it is a function of, and
when it loses that dependence), the simulations behind each claim, and the
figures the manuscript uses.

## What is here

| path | contents |
|---|---|
| `manuscript.tex` | the manuscript, built with `pdflatex` + `bibtex` |
| `figures/` | the plates in `manuscript.tex` |
| `results/` | every table the manuscript quotes |
| `src/issbench/` | the library: spectral decomposition, exact $\E[\mathrm{Iss}]$, simulation, inference wrappers |
| `scripts/` | one script per experiment, numbered in the order the paper reports them |
| `run/` | `nohup` launchers for the long stages |
| `docs/` | working notes and the audit trail |

## Reproducing

Simulation goes through AliSim, inside IQ-TREE 3:

```bash
export IQTREE=/path/to/iqtree3
make gate         # G0: AliSim must pass before anything else is run
make exp1 exp2 exp3 exp4 exp6
python scripts/s07_figures.py
```

Each stage writes `results/<experiment>/cells.csv` one replicate at a time and
resumes from what is already there, so a killed job restarts cleanly. Long
stages are launched with `bash run/<stage>.sh` and log to `logs/`.

The larger results are stored gzipped; `scripts/s20_pack_results.py --unpack`
restores the plain form and `--gzip` puts it back. Every reader goes through
`issbench.results.res`, so either form on disk works.

## Not in this repository

Third-party material is fetched rather than redistributed: the vertebrate COI
alignment (`data/`, described in `data/README.txt`), the Duchene et al.
supplementary files, and the concatenated genome-scale alignment from GigaDB
(`scripts/s09_fetch_gigadb.py`, 4.6 GB). Published PDFs are likewise not
shipped.

## Layout of the argument

The manuscript's Results section runs in the order the claims are established:
the statistic and its normalisation; what the spectral decomposition says it
is a function of; the exact and second-order closed forms; when the index
collapses onto $p$-distance; the dependence of the critical value on sequence
length; resolution collapse near saturation; then the simulations, the COI and
genome-scale applications, and the other tree-free indices ($\lambda$ and the
$C$ value).

## Licence

No licence has been chosen yet. The manuscript and figures are the author's
and all rights are reserved pending publication.
