# Analysis script map

Run scripts from the repository root after setting `SLAMF8_MASLD_ROOT`. The numeric prefixes preserve the historical analysis order; later scripts supersede earlier exploratory implementations where their names overlap.

1. `001`–`010`: doublet detection, lineage-restricted reclustering, trajectory preparation, Monocle3/Slingshot, patient-level validation, communication inference and NicheNet.
2. `012`–`014` and `034`–`037`: communication summaries, SPP1 cross-cohort audits and state/network aggregation.
3. `015`–`025` and `028`–`032`: GEO bulk download, mapping, scoring, confounder checks and full-biopsy validation.
4. `038`–`045`: expanded public-cohort search, eligibility audit and robustness freeze.
5. `054`–`062`: FARG95 self-containment repair, continuous meta-analysis, donor-aware counterfactual re-audit, locked bulk models and specificity/bootstrap audits.
6. `065`–`067`: GSE192741 spatial audit and supplementary output assembly.
7. `033`, `046`–`053`, `057`–`059`, `066`–`067`, `071`: final figure, supplementary table and graphical workflow generation.

Local literature-extraction, manuscript-assembly and submission-packaging utilities are intentionally excluded because they do not generate scientific results and contained workstation-specific paths.

Raw GEO data are not bundled. See `data/metadata/cohorts.tsv` for accessions and `REPRODUCIBILITY.md` for expected folders and interpretation boundaries.
