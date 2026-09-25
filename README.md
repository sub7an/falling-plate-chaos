# constantofchange — chaos in the 2D quasi-steady falling-plate model

[![DOI](https://zenodo.org/badge/1388053043.svg)](https://doi.org/10.5281/zenodo.22967653)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Investigating whether the Andersen–Pesavento–Wang (2005) quasi-steady falling-plate model
has a chaotic regime as the dimensionless moment of inertia **I\*** varies — and if so,
where, by what route, and how sensitive landing displacement is to release conditions.

**Findings are properties of the model, not of real leaves.** All physics comes only from
`docs/model_spec.md`. See `docs/PLAN.md` for the full plan, `docs/paper.md` for the manuscript draft,
`docs/results.md` for the findings synthesis, and `docs/decisions.md` for the chronological
decision log.

## Status
Complete through **M8** (all experiments). Validation ladder (PLAN §4) passes; `python run_all.py`
regenerates the pipeline (E1, E2, E3, E5, E6; E4 skipped — see below). Headline findings, with
their nulls stated honestly (full detail in `docs/decisions.md`):

- **Chaos exists.** E2 gate PASSES: `λ_max` peaks at **+0.140 at I\*=2.2** and, cross-checked by
  the 0–1 test (median K), confirms a chaotic band **I\* ≈ [1.63, 2.83]**. Regimes run
  flutter → tumble (heteroclinic ≈1.2191) → period-2 → chaos → broadside flutter.
- **Route is NOT Feigenbaum (null).** Only 2 successive period-doublings (1→2 at ≈1.474, 2→4 at
  ≈1.599), then windowed chaos; power spectra rule out a torus. **E4 (Feigenbaum δ) is not run** —
  a δ from 2 doublings would overclaim.
- **Displacement sensitivity.** In the chaotic band σ(D)=std(Δx) grows exponentially; measured in
  *time*, the rate `d(log σ)/dt` equals **λ_max** within ~10% (estimation uncertainty) — the
  Lyapunov exponent sets the displacement-sensitivity growth rate. (The looser depth-based
  comparison to `λ_max/U` carries an extra depth↔time conversion; resolved 2026-09-25, see decisions
  log.) Regular bands are flat (insensitive), not linearly growing.
- **No multistability (null).** E6: no hysteresis and a single attractor at every I\* tested
  (only the trivial CW/CCW reflection-symmetry pair).
- **One anchor discrepancy:** at I\*=1.45 this model is period-**one** tumbling, not the paper's
  figure-read period-two (the first doubling sits at ≈1.474 here).

## Setup
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Reproduce
One command regenerates every cached figure from cache or from scratch:
```bash
python run_all.py
```
Build the manuscript PDF (regenerates figures first):
```bash
brew install tectonic   # one-time: self-contained LaTeX engine (no full MacTeX needed)
make paper              # -> docs/paper.pdf  (from docs/paper.tex, figures embedded)
```
`make paper` prefers `tectonic`, falling back to `latexmk`/`pdflatex` if a system TeX is present.
Build a flat arXiv upload with `make arxiv` (`-> docs/paper-arxiv.zip`). See
`docs/submission_checklist.md` and `docs/cover_letter.md` for the submission path.
The write-up exists in two forms: `docs/paper.md` (readable Markdown draft) and
`docs/paper.tex` (compilable LaTeX, `article` class, figures embedded).

## Layout
```
docs/         PLAN.md, decisions.md, model_spec.md (pending)
fallingplate/ package: model, integrate, regimes, lyapunov, chaos_tests,
              displacement, periodic_orbits, sweeps, io_cache, plotting
experiments/  e1..e6 drivers (one per experiment, YAML-configured)
config/       experiment configs (no magic numbers in code)
tests/        validation ladder
data/         cached arrays + provenance sidecars (gitignored)
figures/      generated figures (gitignored)
```

## Provenance
Every output carries a sidecar recording git hash, Python + SciPy versions, parameters,
and seed. The cache key is the parameter hash plus a hash of the physics/numerics source
files.

## Development note
Developed by the author with the assistance of general-purpose AI tools, under the author's
direction. The author defined the study, made all scientific decisions, independently verified the
validation anchors (fixed points, edge-on eigenvalues, the flutter/tumble transition) and the
Lyapunov code (on the logistic map and the Lorenz system), directed the rigor checks at every stage,
and is responsible for all results and their interpretation. Every decision and correction is logged
in `docs/decisions.md`.

## Citing this work
If you use this code or its results, please cite it via `CITATION.cff` (GitHub shows a "Cite this
repository" button) or the archived DOI once the release is minted (see below).

## License
Code is released under the MIT License (`LICENSE`); the manuscript text and figures
(`docs/paper.*`, `figures/*`) are under CC BY 4.0.
