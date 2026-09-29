from __future__ import annotations
import os

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu, norm, spearmanr


PROJECT_ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
PROCESSED = PROJECT_ROOT / "data_processed"
RESULTS = PROJECT_ROOT / "results"
FIGURES = PROJECT_ROOT / "figures"
DATASETS = ["GSE344087", "GSE298719", "GSE212837"]
COMPARISONS = {
    "GSE344087": ("fibrosis", "fibrosis", "no fibrosis"),
    "GSE298719": ("disease", "MASLD", "Healthy"),
    "GSE212837": ("disease", "NASH", "Control"),
}
TARGET_GENES = ["YAP1", "NUAK2", "COL1A1", "COL1A2", "COL3A1", "LOX", "ACTA2", "TAGLN", "TGFB1", "PDGFRA", "PDGFRB"]


def bh_adjust(values: pd.Series) -> pd.Series:
    order = np.argsort(values.to_numpy())
    ranked = values.to_numpy()[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1)
    return pd.Series(output, index=values.index)


def valid_hsc_mask(obs: pd.DataFrame) -> np.ndarray:
    return ~obs["subtype"].astype(str).str.contains("contaminant|low_quality|ambiguous", case=False, regex=True).to_numpy()


def sample_column(obs: pd.DataFrame) -> str:
    return "patient_id" if "patient_id" in obs.columns else "gsm"


def summarize_groups(adata: ad.AnnData, dataset: str) -> pd.DataFrame:
    adata = adata[valid_hsc_mask(adata.obs)].copy()
    sample_col = sample_column(adata.obs)
    condition_col, _, _ = COMPARISONS[dataset]
    genes = [gene for gene in TARGET_GENES if gene in adata.var_names]
    gene_positions = {gene: adata.var_names.get_loc(gene) for gene in genes}
    counts = adata.layers["counts"]
    module_columns = [column for column in adata.obs.columns if column.startswith("module_")]
    rows = []

    def add_group(sample: str, subtype: str, positions: np.ndarray) -> None:
        library_size = float(counts[positions, :].sum())
        row = {
            sample_col: sample,
            "subtype": subtype,
            "n_cells": len(positions),
            "library_size": library_size,
        }
        for gene, position in gene_positions.items():
            gene_count = float(counts[positions, position].sum())
            row[f"logCPM_{gene}"] = np.log1p(1e6 * gene_count / library_size) if library_size > 0 else np.nan
        for column in module_columns:
            row[column] = float(adata.obs.iloc[positions][column].mean())
        rows.append(row)

    for (sample, subtype), positions in adata.obs.groupby([sample_col, "subtype"], observed=True).indices.items():
        add_group(sample, str(subtype), np.asarray(positions))
    for sample, positions in adata.obs.groupby(sample_col, observed=True).indices.items():
        add_group(sample, "All_HSC", np.asarray(positions))

    values = pd.DataFrame(rows)
    sample_metadata = adata.obs[[sample_col, condition_col]].drop_duplicates(subset=[sample_col]).set_index(sample_col)
    values[condition_col] = values[sample_col].map(sample_metadata[condition_col])
    return values


def test_conditions(values: pd.DataFrame, dataset: str) -> pd.DataFrame:
    condition_col, case_label, control_label = COMPARISONS[dataset]
    features = [column for column in values.columns if column.startswith("logCPM_") or column == "module_HSC_activation"]
    rows = []
    for subtype, group in values.groupby("subtype"):
        eligible = group.loc[group["n_cells"] >= 20]
        case = eligible.loc[eligible[condition_col].astype(str).str.lower() == case_label.lower()]
        control = eligible.loc[eligible[condition_col].astype(str).str.lower() == control_label.lower()]
        if len(case) < 3 or len(control) < 3:
            continue
        for feature in features:
            case_values = case[feature].dropna().to_numpy(dtype=float)
            control_values = control[feature].dropna().to_numpy(dtype=float)
            if len(case_values) < 3 or len(control_values) < 3:
                continue
            statistic, pvalue = mannwhitneyu(case_values, control_values, alternative="two-sided")
            rows.append(
                {
                    "dataset": dataset,
                    "subtype": subtype,
                    "feature": feature,
                    "n_case": len(case_values),
                    "n_control": len(control_values),
                    "median_difference": np.median(case_values) - np.median(control_values),
                    "rank_biserial": 2 * statistic / (len(case_values) * len(control_values)) - 1,
                    "p_value": pvalue,
                }
            )
    output = pd.DataFrame(rows)
    if len(output):
        output["fdr"] = bh_adjust(output["p_value"])
    return output


def link_macrophage_hsc(dataset: str, hsc_values: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    macrophage = pd.read_csv(RESULTS / f"013_{dataset}_patient_subtype_values.tsv.gz", sep="\t")
    macrophage = macrophage.loc[macrophage["subtype"] == "All_macrophage"].copy()
    hsc = hsc_values.loc[hsc_values["subtype"] == "All_HSC"].copy()
    unit = "patient_id" if "patient_id" in hsc.columns else "gsm"
    merged = macrophage.merge(hsc, on=unit, suffixes=("_macrophage", "_HSC"))
    features = [column for column in hsc.columns if column.startswith("logCPM_") or column == "module_HSC_activation"]
    rows = []
    for feature in features:
        target_column = feature if feature in merged.columns else f"{feature}_HSC"
        x = merged["SLAMF8_logCPM"].to_numpy(dtype=float)
        y = merged[target_column].to_numpy(dtype=float)
        keep = np.isfinite(x) & np.isfinite(y)
        if keep.sum() < 5:
            continue
        correlation, pvalue = spearmanr(x[keep], y[keep])
        rows.append({"dataset": dataset, "feature": feature, "n_samples": int(keep.sum()), "spearman_rho": correlation, "p_value": pvalue})
    tests = pd.DataFrame(rows)
    if len(tests):
        tests["fdr"] = bh_adjust(tests["p_value"])
    return merged, tests


def main() -> None:
    correlation_frames = []
    scatter_frames = []
    for dataset in DATASETS:
        adata = ad.read_h5ad(PROCESSED / f"{dataset}_hsc_reclustered.h5ad")
        values = summarize_groups(adata, dataset)
        tests = test_conditions(values, dataset)
        merged, correlations = link_macrophage_hsc(dataset, values)
        values.to_csv(RESULTS / f"029_{dataset}_HSC_patient_values.tsv.gz", sep="\t", index=False, compression="gzip")
        tests.to_csv(RESULTS / f"030_{dataset}_HSC_patient_tests.tsv", sep="\t", index=False)
        merged.to_csv(RESULTS / f"031_{dataset}_macrophage_HSC_matched_values.tsv", sep="\t", index=False)
        correlations.to_csv(RESULTS / f"032_{dataset}_macrophage_HSC_correlations.tsv", sep="\t", index=False)
        correlation_frames.append(correlations)
        hsc_activation_column = "module_HSC_activation" if "module_HSC_activation" in merged else "module_HSC_activation_HSC"
        if hsc_activation_column in merged:
            scatter_frames.append(pd.DataFrame({"dataset": dataset, "SLAMF8_logCPM": merged["SLAMF8_logCPM"], "HSC_activation": merged[hsc_activation_column]}))

    combined = pd.concat(correlation_frames, ignore_index=True)
    combined.to_csv(RESULTS / "033_cross_dataset_macrophage_HSC_correlations.tsv", sep="\t", index=False)
    meta_rows = []
    for feature, group in combined.groupby("feature"):
        correlations = group["spearman_rho"].clip(-0.999999, 0.999999).to_numpy(dtype=float)
        weights = np.maximum(group["n_samples"].to_numpy(dtype=float) - 3, 1)
        pooled_z = np.sum(weights * np.arctanh(correlations)) / np.sum(weights)
        statistic = pooled_z * np.sqrt(np.sum(weights))
        meta_rows.append(
            {
                "feature": feature,
                "n_datasets": len(group),
                "n_samples_total": int(group["n_samples"].sum()),
                "positive_datasets": int((group["spearman_rho"] > 0).sum()),
                "pooled_spearman_rho": float(np.tanh(pooled_z)),
                "z_statistic": float(statistic),
                "p_value": float(2 * norm.sf(abs(statistic))),
            }
        )
    meta = pd.DataFrame(meta_rows)
    meta["fdr"] = bh_adjust(meta["p_value"])
    meta.sort_values("p_value").to_csv(RESULTS / "034_cross_dataset_macrophage_HSC_correlation_meta.tsv", sep="\t", index=False)
    if scatter_frames:
        plot_data = pd.concat(scatter_frames, ignore_index=True)
        sns.set_theme(style="whitegrid", context="notebook")
        grid = sns.lmplot(data=plot_data, x="SLAMF8_logCPM", y="HSC_activation", col="dataset", col_wrap=3, height=4, scatter_kws={"s": 35, "alpha": 0.8}, ci=None)
        grid.set_axis_labels("Macrophage SLAMF8 logCPM", "HSC activation module")
        grid.fig.suptitle("Matched-patient macrophage–HSC state association", y=1.04)
        grid.savefig(FIGURES / "010_cross_dataset_macrophage_SLAMF8_HSC_activation.png", dpi=240)
        plt.close("all")


if __name__ == "__main__":
    main()
