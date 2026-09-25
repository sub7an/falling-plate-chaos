# References

PDFs are **not tracked in git** (copyrighted journal binaries; see `.gitignore`).
Place the source PDFs here locally.

## Primary source
- **[P91]** Andersen, A., Pesavento, U. & Wang, Z. J. (2005). *Analysis of transitions
  between fluttering, tumbling and steady descent of falling cards.* J. Fluid Mech. **541**,
  91–104. This is the source of every equation and coefficient in `docs/model_spec.md`.
  - Local file: `2005_JFM_pp91-104_transitions.pdf` (verified via PDF metadata and pages
    91–104). Originally named `..._b.pdf`: Cornell's a/b file-suffix convention is the reverse
    of the usual 2005a/2005b citation convention, so the content is correct — renamed by page
    range to avoid the ambiguity.

## Companion paper (not currently in repo)
- **[P65]** Andersen, A., Pesavento, U. & Wang, Z. J. (2005). *Unsteady aerodynamics of
  fluttering and tumbling plates.* J. Fluid Mech. **541**, 65–90. Dimensional formulation and
  experiments. Per `docs/model_spec.md` §0 it is **not** used for the equations.

## Second-reader verification
`docs/model_spec.md` was cross-checked against [P91] pages 98–102 on 2026-09-22; result and
recovered items (incl. the Hopf curve 6.3) are logged in `docs/decisions.md`.
