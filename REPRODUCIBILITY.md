# Reproducibility notes

## Required directory convention

The portable scripts resolve the project root from `SLAMF8_MASLD_ROOT`; when unset, they use the parent directory of `scripts/`. Large downloaded inputs are intentionally excluded from version control.

Expected working directories created or consumed by the pipeline include:

- `data_raw/`: downloaded raw or supplementary GEO files.
- `data_processed/`: normalized, mapped and lineage-specific intermediate objects.
- `results/`: frozen statistical outputs and audit tables.
- `figures/`: generated manuscript figures.
- `supplementary_tables/`: generated S1–S14 files.

## Determinism

The principal stochastic analyses use seed `20260920`; the paired bootstrap comparison uses seed `20260923`. Package versions are frozen under `environment/`. Exact numerical identity can still depend on BLAS, operating system and package compilation.

## Interpretation boundaries

- Pseudotime represents transcriptional state ordering, not lineage tracing or observed time.
- LIANA, CellChat and NicheNet provide expression-supported candidate communication, not protein-level or causal proof.
- The GRNBoost2/ridge counterfactual is a model-based prediction, not an experimental knockout.
- Visium correlations are spot-level co-expression audits and do not prove spatial proximity or signaling direction.
- Bulk AUCs are descriptive tissue-state discrimination estimates, not trained or calibrated clinical models.

## FARG95 provenance

The ferro-aging concept is traced to Liu et al., *Cell Metabolism* (2026), PMID 41819088, DOI 10.1016/j.cmet.2026.02.010. The 95-gene list was transcribed from Supplementary Table 1 of Lin et al., *Biology Direct* (2026), DOI 10.1186/s13062-026-00956-4. Only SLAMF8 among the prespecified five anchors is a member of FARG95; the anchor-excluded and SLAMF8-excluded scores are therefore identical.
