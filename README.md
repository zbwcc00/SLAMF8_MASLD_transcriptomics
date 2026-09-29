# SLAMF8-marked macrophage program in human MASLD

This repository contains the analysis code, processed figure-level Source Data, Supplementary Tables and final vector figures supporting the manuscript:

> Cross-cohort transcriptomics identifies a SLAMF8-marked scar-associated macrophage program linked to hepatic stellate-cell activation and fibrosis in human MASLD

## Study overview

The workflow integrates three public human liver single-cell/single-nucleus cohorts, two primary independent bulk-biopsy validation cohorts, additional public validation cohorts and one Visium spatial transcriptomic audit. It evaluates a SLAMF8-marked scar-associated macrophage state, its ferro-aging-related transcriptional overlay, macrophage–HSC associations, SPP1-associated communication candidates, transcriptional ordering and donor-aware counterfactual network predictions.

The repository contains no new patient recruitment and no redistributed raw human sequencing data. Raw matrices remain available from NCBI GEO under the accession numbers in `data/metadata/cohorts.tsv`.

## Repository contents

- `scripts/`: 65 portable Python/R analysis and figure scripts.
- `data/source_data/`: figure-level processed Source Data.
- `data/supplementary_tables/`: Supplementary Tables S1–S14 in XLSX/TSV form with readme files.
- `data/reference/`: frozen FARG95 membership and provenance.
- `data/metadata/cohorts.tsv`: accession-level dataset inventory and stable GEO links.
- `figures/`: final PDF versions of Figure 1–5 and Figure S1–S12.
- `environment/`: pinned Python environment and exact R package-version inventory.
- `audit/`: consistency and figure-QA summaries.

## Quick start

```bash
git clone https://github.com/zbwcc00/SLAMF8_MASLD_transcriptomics.git
cd SLAMF8_MASLD_transcriptomics
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
python -m pip install -r environment/requirements.txt
python tools/validate_repository.py
```

Set `SLAMF8_MASLD_ROOT` to the repository root before running analysis scripts. Set `SLAMF8_MASLD_SOURCE_ROOT` only when raw GEO downloads are stored elsewhere. Scripts default to the repository root and expect large, non-redistributed inputs under `data_raw/` or the paths documented in `scripts/RUN_ORDER.md`.

## Reproducibility scope

The included Source Data and Supplementary Tables permit inspection of all reported values without downloading the large raw datasets. Full regeneration from raw matrices requires downloading the cited GEO records and can require substantial memory, storage and runtime. Communication and trajectory packages may also show platform-dependent numerical variation; preserve the stated seeds and versions.

## Availability

The versioned release is archived in Zenodo at https://doi.org/10.5281/zenodo.23042100 and mirrored at https://github.com/zbwcc00/SLAMF8_MASLD_transcriptomics/releases/tag/v1.0.0.

## Authors

Bowen Zheng, Guanghua Xie, Wangde Jin and Hao Li (corresponding author).  
Division of Hepatobiliary Pancreatic Surgery, The Affiliated Hospital of Yanbian University, Yanji 133000, China

## Licences

Analysis code is provided under the MIT License. Author-generated processed data and documentation are provided under CC BY 4.0 where the authors hold the necessary rights. Third-party raw datasets are not redistributed and remain subject to their source repository terms.
