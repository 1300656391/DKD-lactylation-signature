# DKD-lactylation-signature

Analysis code for the manuscript *"A lactate–lactylation axis gene signature for diabetic
kidney disease: diagnostic model development, six-cohort evaluation, and immune
microenvironment correlates"* (Frontiers in Immunology, submitted).

## Repository layout

| Directory | Contents |
|---|---|
| [clean_pipeline/](clean_pipeline/) | End-to-end analysis pipeline (R scripts 00–13, `master.R`, Python figure assembly). Reproduces every main-text figure and table and all supplementary tables from public GEO data. |
| [additional_analyses/](additional_analyses/) | Standalone sensitivity and extension analyses (Tables S10, S11, S15a/b, S16–S18, head-to-head bootstrap). |

## Quick start

See **[clean_pipeline/README.md](clean_pipeline/README.md)** for the four-step guide:
install packages → download GEO inputs → run `master.R` → run `assemble_figures.py`.
All results derive from public GEO datasets under a fixed seed (42); the frozen model
coefficients reproduce all reported statistics exactly.

## License

MIT — see [clean_pipeline/LICENSE](clean_pipeline/LICENSE).
