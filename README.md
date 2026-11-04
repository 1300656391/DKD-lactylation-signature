# DKD-lactylation-signature

Analysis code for the PDK2/PDK4 lactate-axis signature study in diabetic kidney disease (six-cohort validation with orthogonal protein, DNA-methylation, and perturbation layers).

## Repository layout

| Directory | Contents |
|---|---|
| [pipeline/](pipeline/) | End-to-end analysis pipeline (R scripts 00–13, `master.R`, Python figure assembly). Reproduces all main-text figures and tables and the pipeline-generated supplementary tables (S5, S6, S12–S16, S18–S21; S2–S4, S7–S9 are direct copies of primary outputs) from public GEO data. |
| [additional_analyses/](additional_analyses/) | Standalone sensitivity and extension analyses (Tables S10, S11, S15a/b, S17; DeLong comparisons; head-to-head paired bootstrap; NNLS reference-rebuild sensitivity), shipped with their scripts and frozen outputs. |
| [orthogonal_layers/](orthogonal_layers/) | Orthogonal protein (two public DKD kidney proteomes), DNA-methylation (proximal-tubule methylome, primary plus replication), and perturbation-layer (iLINCS/L1000) analyses (Sections 2.9–2.10 and 3.9–3.11; Tables S22–S27), shipped with scripts and frozen outputs; `orthogonal_layers/assemble_esm_tables.py` rebuilds those tables byte-for-byte. |

## Quick start

See **[pipeline/README.md](pipeline/README.md)** for the four-step guide:
install packages → download GEO inputs → run `master.R` → run `assemble_figures.py`.
All results derive from public GEO datasets under a fixed seed (42); the frozen model
coefficients reproduce all reported statistics exactly.

## License

MIT — see [LICENSE](LICENSE) (also [pipeline/LICENSE](pipeline/LICENSE)).
