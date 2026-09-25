# Submission checklist

A practical path from "finished draft" to "public, citable result." You are the author; the
manuscript is `docs/paper.tex` → `docs/paper.pdf` (build with `make paper`).

## 0. Two free things to do first

- [ ] **Get an ORCID** at <https://orcid.org> (2 minutes, free). It is a permanent researcher ID you
      will reuse for years. Put it in the author block and cover letter.
- [ ] **Keep the GitHub repo public** — it is your "code and data availability" backing. The paper
      already claims full reproducibility (`python run_all.py`, `make paper`, `pytest`).

## 1. Fill the author details (only things left blank)

- [ ] Email and ORCID in the author block of `paper.tex` (currently name + "Independent researcher").
- [ ] If you have an affiliation (e.g. once at university), replace "Independent researcher".

## 2. Choose a venue

| Venue | Good for | Notes |
|-------|----------|-------|
| **arXiv** (preprint) | A public, timestamped, citable version | Not peer-reviewed. Needs an *endorsement* for `nlin.CD` (see §3). Do this first regardless of journal. |
| **Journal of Emerging Investigators (JEI)** | High-school / pre-university research, real peer review + mentoring | Requires restructuring to JEI's format (Abstract, Introduction, Results, Discussion, Methods-at-end). Ask me to reformat. |
| A specialised journal later | A formal publication | Would want a sharper novelty case vs the 2005 paper and comparison with newer literature (e.g. Xu et al. 2021). |

## 3. arXiv submission (recommended first step)

- [ ] Create an account at <https://arxiv.org>.
- [ ] **Endorsement:** first-time submitters to a category usually need an endorser — a researcher
      who has posted several papers to `nlin` recently. Ask a mentor/teacher/contact, or use arXiv's
      endorsement request flow. (An academic email can sometimes auto-endorse.)
- [ ] Primary category **`nlin.CD`** (Chaotic Dynamics); cross-list **`physics.flu-dyn`**.
- [ ] Build the upload package: `make arxiv` → `docs/paper-arxiv.zip` (flat: `paper.tex` + the 5
      figure PNGs). arXiv compiles it with pdfLaTeX; the packages used (amsmath, graphicx, geometry,
      hyperref, booktabs, microtype, lmodern) are all standard on arXiv.
- [ ] Pick a license (CC BY 4.0 is a good open default).
- [ ] Paste the abstract; set the title and author exactly as in the paper.

## 4. Pre-submission checks (all currently pass)

- [x] Compiles clean (`make paper`): 11 pp, 0 undefined references, no overfull boxes.
- [x] All 5 figures embedded; 2 summary tables present.
- [x] No placeholders in the text; references complete.
- [x] AI-use disclosed per policy (Acknowledgements) — keep this; venues require it.
- [ ] Final proofread (spelling, a colleague's read-through).
- [ ] Verify reference details against originals (they are standard, but confirm volumes/pages).

## 5. After it is public

- [ ] Add the arXiv ID / DOI to the README.
- [ ] If submitting to a journal, use `docs/cover_letter.md` (fill the bracketed fields).
