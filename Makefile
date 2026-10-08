# Iss audit / critique project -- convenience targets.
# LaTeX is NOT installed on the current machine; the `ms` target is for machines
# that have pdflatex + bibtex.

PY      := python
IQTREE  ?= E:/Codes/YRTools/plugins/YR-MPE/YR_MPE/bin/iqtree3/iqtree3.exe

.PHONY: help gate design exp1 exp2 exp3 exp4 exp6 figs bib clean \
        ms ms-watch check

help:
	@echo "make gate     - G0 pipeline gate (AliSim must pass before anything else)"
	@echo "make bib      - regenerate refs/refs_iss_clean.bib from refs_iss.bib"
	@echo "make design   - Item 0: uniform-Iss sampling density"
	@echo "make exp1     - Item 1: 2-D maps of Iss location and scale"
	@echo "make exp2     - Item 2: Iss vs p-distance, closed-form overlay"
	@echo "make exp3     - Item 3: Felsenstein->Farris transfer"
	@echo "make exp4     - Item 4: one-at-a-time perturbation trends"
	@echo "make exp6     - Item 6: conditional independence Pr(Y|pbar,Iss) vs Pr(Y|pbar)"
	@echo "make figs     - Item 5: plates in the ref style"
	@echo "make check    - static checks on manuscript.tex (no LaTeX needed)"
	@echo "make ms       - build manuscript.pdf (requires pdflatex + bibtex)"

# ---------------------------------------------------------------- simulation
gate:
	$(PY) scripts/s00_gate.py

design:
	$(PY) scripts/s01_design.py

exp1:
	$(PY) scripts/s02_Lmap.py

exp2:
	$(PY) scripts/s03_iss_pdist.py

exp3:
	$(PY) scripts/s04_zone.py

exp4:
	$(PY) scripts/s05_perturb.py

exp6:
	$(PY) scripts/s06_condindep.py
	$(PY) scripts/s06_condindep.py --analyze

figs:
	$(PY) scripts/s07_figures.py

# ------------------------------------------------------------------ manuscript
bib:
	$(PY) scripts/make_clean_bib.py

check:
	$(PY) scripts/check_tex.py

ms: bib
	pdflatex -interaction=nonstopmode manuscript.tex
	bibtex manuscript
	pdflatex -interaction=nonstopmode manuscript.tex
	pdflatex -interaction=nonstopmode manuscript.tex

ms-watch:
	latexmk -pdf -pvc manuscript.tex

clean:
	rm -f manuscript.aux manuscript.bbl manuscript.blg manuscript.log \
	      manuscript.out manuscript.pdf manuscript.fls manuscript.fdb_latexmk
