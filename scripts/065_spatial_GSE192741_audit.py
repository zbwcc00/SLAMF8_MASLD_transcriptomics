#!/usr/bin/env python3
import os
"""Audit an independent human liver Visium cohort for the locked gene panels."""

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
DATA = ROOT / "data_processed" / "spatial_geo" / "GSE192741"
OUT = ROOT / "results" / "spatial_validation" / "GSE192741"
OUT.mkdir(parents=True, exist_ok=True)

FILES = {
    "GSM5764424": ("GSM5764424_filtered_feature_bc_matrix_JBO014.h5", "steatotic"),
    "GSM5764425": ("GSM5764425_filtered_feature_bc_matrix_JBO015.h5", "steatotic"),
    "GSM5764426": ("GSM5764426_filtered_feature_bc_matrix_JBO018.h5", "healthy"),
    "GSM5764427": ("GSM5764427_filtered_feature_bc_matrix_JBO019.h5", "steatotic"),
    "GSM5764428": ("GSM5764428_filtered_feature_bc_matrix_JBO022.h5", "healthy"),
}
PANELS = {
    "SAMac_like": ["TREM2", "GPNMB", "CD9", "SPP1", "LGALS3", "LPL", "APOC1", "FABP5", "CTSB", "CTSD"],
    "SAMac_like_no_SPP1": ["TREM2", "GPNMB", "CD9", "LGALS3", "LPL", "APOC1", "FABP5", "CTSB", "CTSD"],
    "HSC_fibrogenic": ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "DCN", "LUM", "ACTA2", "TAGLN", "PDGFRB"],
    "HSC_activation": ["YAP1", "NUAK2", "COL1A1", "COL1A2", "COL3A1", "LOX", "ACTA2", "TAGLN", "TGFB1", "PDGFRA", "PDGFRB"],
}
TARGETS = sorted(set(sum(PANELS.values(), []) + ["SLAMF8", "SPP1", "CD68", "CD163", "ITGAV", "ITGB5"]))


def score(data: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    found = [gene for gene in genes if gene in data.columns]
    if len(found) < 3:
        return pd.Series(np.nan, index=data.index), found
    values = data[found].astype(float)
    z = (values - values.mean()) / values.std(ddof=0).replace(0, np.nan)
    return z.mean(axis=1), found


def main() -> None:
    spot_rows, gene_rows, corr_rows = [], [], []
    for sample, (filename, condition) in FILES.items():
        adata = sc.read_10x_h5(DATA / filename)
        adata.var_names = adata.var_names.astype(str)
        adata.var_names_make_unique()
        counts = adata[:, [gene for gene in TARGETS if gene in adata.var_names]].to_df()
        log_expr = np.log1p(counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0) * 1e4)
        row = {"sample": sample, "condition": condition, "n_spots": len(log_expr)}
        for gene in TARGETS:
            row[f"detected_{gene}"] = int((log_expr[gene] > 0).sum()) if gene in log_expr else 0
            row[f"mean_{gene}"] = float(log_expr[gene].mean()) if gene in log_expr else np.nan
        for panel, genes in PANELS.items():
            values, found = score(log_expr, genes)
            row[f"mean_{panel}"] = float(values.mean())
            row[f"n_genes_{panel}"] = len(found)
        spot_rows.append(row)
        data = log_expr.copy()
        for panel, genes in PANELS.items():
            data[panel], _ = score(log_expr, genes)
        for gene in ["SLAMF8", "SPP1", "CD68", "CD163", "ITGAV", "ITGB5"]:
            for panel in ["SAMac_like", "SAMac_like_no_SPP1", "HSC_fibrogenic", "HSC_activation"]:
                if gene in data and panel in data:
                    paired = data[[gene, panel]].replace([np.inf, -np.inf], np.nan).dropna()
                    rho, pval = spearmanr(paired[gene], paired[panel]) if len(paired) >= 3 else (np.nan, np.nan)
                    corr_rows.append({"sample": sample, "condition": condition, "gene": gene, "panel": panel, "n_spots": len(data), "n_paired_spots": len(paired), "rho": rho, "p_value": pval})
        for spot, values in data.iterrows():
            spot_rows.append({"sample": sample, "condition": condition, "spot": spot, "SLAMF8": values.get("SLAMF8", np.nan), "SPP1": values.get("SPP1", np.nan), "SAMac_like": values.get("SAMac_like", np.nan), "HSC_fibrogenic": values.get("HSC_fibrogenic", np.nan), "HSC_activation": values.get("HSC_activation", np.nan)})
    summary = pd.DataFrame([row for row in spot_rows if "spot" not in row])
    summary.to_csv(OUT / "sample_level_spatial_summary.tsv", sep="\t", index=False)
    pd.DataFrame(corr_rows).to_csv(OUT / "spot_level_gene_panel_correlations.tsv", sep="\t", index=False)
    corr = pd.DataFrame(corr_rows)
    if not corr.empty:
        corr["direction"] = np.where(corr["rho"] >= 0, "positive", "negative")
    corr.to_csv(OUT / "spot_level_spatial_audit.tsv", sep="\t", index=False)
    (OUT / "README.md").write_text(
        "# GSE192741 spatial audit\n\n"
        "Five human Visium liver samples were audited at spot level. Three samples were labelled steatotic and two healthy in GEO metadata. The audit tests expression and relative panel co-variation; it is not a cell-resolved lineage or causal analysis. Because the data contain spot-level mixtures and no SLAMF8 protein assay, positive co-variation is spatial support for follow-up rather than mechanistic validation.\n",
        encoding="utf-8",
    )
    print(summary[["sample", "condition", "n_spots", "mean_SLAMF8", "mean_SPP1", "mean_SAMac_like", "mean_HSC_fibrogenic"]].to_string(index=False))


if __name__ == "__main__":
    main()
