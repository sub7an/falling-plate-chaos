# Thin convenience wrapper. Canonical entry point is run_all.py (portable).
PYTHON ?= ./.venv/bin/python

.PHONY: figures test list paper clean-paper

figures:
	$(PYTHON) run_all.py

list:
	$(PYTHON) run_all.py --list

test:
	$(PYTHON) -m pytest -q

# Build the manuscript PDF -> docs/paper.pdf. Figures are regenerated first.
# Prefers tectonic (self-contained, no full TeX install: `brew install tectonic`);
# falls back to latexmk/pdflatex if a system TeX is present.
paper: figures
	cd docs && ( tectonic paper.tex \
	  || latexmk -pdf -silent paper.tex \
	  || pdflatex -interaction=nonstopmode paper.tex \
	  || { echo "No LaTeX engine found. Install one, e.g.: brew install tectonic"; exit 1; } )

clean-paper:
	cd docs && rm -f paper.aux paper.log paper.out paper.fls paper.fdb_latexmk paper.pdf

# Assemble a flat arXiv upload: paper.tex + the 5 figures + a README, zipped.
FIGS = e1_regime_map e2_lyapunov e3_route e5_sensitivity e6_multistability
arxiv: figures
	rm -rf docs/arxiv docs/paper-arxiv.zip
	mkdir -p docs/arxiv
	cp docs/paper.tex docs/arxiv/
	$(foreach f,$(FIGS),cp figures/$(f).png docs/arxiv/;)
	printf '%s\n' "arXiv source: build with 'pdflatex paper.tex' (or tectonic). Figures are flat here." > docs/arxiv/00README.txt
	cd docs/arxiv && zip -q ../paper-arxiv.zip paper.tex 00README.txt *.png
	@echo "arXiv package -> docs/paper-arxiv.zip  (paper.tex + 5 figures, compiles flat)"
