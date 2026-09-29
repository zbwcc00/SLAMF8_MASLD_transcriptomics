from __future__ import annotations
import os

import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse
from scipy.stats import mannwhitneyu, wilcoxon


PROJECT_ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
PROCESSED = PROJECT_ROOT / "data_processed"
RESULTS = PROJECT_ROOT / "results"
FIGURES = PROJECT_ROOT / "figures"

COMPARISONS = {
    "GSE344087": ("fibrosis", "fibrosis", "no fibrosis"),
    "GSE298719": ("disease", "MASLD", "Healthy"),
    "GSE212837": ("disease", "NASH", "Control"),
}


def bh_adjust(values: pd.Series) -> pd.Series:
    values = values.astype(float)
    order = np.argsort(values.to_numpy())
    ranked = values.to_numpy()[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1)
    return pd.Series(output, index=values.index)


def valid_lineage_mask(obs: pd.DataFrame) -> np.ndarray:
    excluded = "cDC|pDC|_DC|Cycling|ambiguous|contaminant|Neutrophil|low_quality"
    return ~obs["subtype"].astype(str).str.contains(excluded, case=False, regex=True).to_numpy()


def sample_column(obs: pd.DataFrame) -> str:
    return "patient_id" if "patient_id" in obs.columns else "gsm"


def pseudobulk_slamf8(adata: ad.AnnData, sample_col: str) -> pd.DataFrame:
    if "SLAMF8" not in adata.var_names:
        raise ValueError("SLAMF8 absent")
    gene_position = adata.var_names.get_loc("SLAMF8")
    counts = adata.layers["counts"]
    rows = []
    grouping = adata.obs.groupby([sample_col, "subtype"], observed=True).indices
    for (sample, subtype), positions in grouping.items():
        total_counts = float(counts[positions, :].sum())
        slamf8_counts = float(counts[positions, gene_position].sum())
        rows.append(
            {
                sample_col: sample,
                "subtype": str(subtype),
                "n_cells": len(positions),
                "library_size": total_counts,
                "SLAMF8_counts": slamf8_counts,
                "SLAMF8_logCPM": np.log1p(1e6 * slamf8_counts / total_counts) if total_counts > 0 else np.nan,
            }
        )
    return pd.DataFrame(rows)


def build_sample_values(adata: ad.AnnData, dataset: str) -> pd.DataFrame:
    adata = adata[valid_lineage_mask(adata.obs)].copy()
    sample_col = sample_column(adata.obs)
    condition_col, _, _ = COMPARISONS[dataset]
    module_columns = [column for column in adata.obs.columns if column.startswith("module_")]
    group_columns = [sample_col, "subtype"]
    means = adata.obs.groupby(group_columns, observed=True)[module_columns].mean().reset_index()
    expression = pseudobulk_slamf8(adata, sample_col)
    values = expression.merge(means, on=group_columns, how="left")

    original_subtypes = adata.obs["subtype"].copy()
    adata.obs["subtype"] = "All_macrophage"
    pooled_expression = pseudobulk_slamf8(adata, sample_col)
    pooled_means = adata.obs.groupby(group_columns, observed=True)[module_columns].mean().reset_index()
    pooled_values = pooled_expression.merge(pooled_means, on=group_columns, how="left")
    adata.obs["subtype"] = original_subtypes
    values = pd.concat([values, pooled_values], ignore_index=True)

    sample_metadata = adata.obs[[sample_col, condition_col]].drop_duplicates(subset=[sample_col]).set_index(sample_col)
    values[condition_col] = values[sample_col].map(sample_metadata[condition_col])
    total_cells = adata.obs.groupby(sample_col, observed=True).size()
    values["myeloid_fraction"] = values["n_cells"] / values[sample_col].map(total_cells)
    return values


def condition_tests(values: pd.DataFrame, dataset: str) -> pd.DataFrame:
    condition_col, case_label, control_label = COMPARISONS[dataset]
    features = ["myeloid_fraction", "SLAMF8_logCPM"] + [column for column in values.columns if column.startswith("module_")]
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
                    "case": case_label,
                    "control": control_label,
                    "n_case": len(case_values),
                    "n_control": len(control_values),
                    "median_case": np.median(case_values),
                    "median_control": np.median(control_values),
                    "median_difference": np.median(case_values) - np.median(control_values),
                    "rank_biserial": 2 * statistic / (len(case_values) * len(control_values)) - 1,
                    "p_value": pvalue,
                }
            )
    output = pd.DataFrame(rows)
    if len(output):
        output["fdr"] = bh_adjust(output["p_value"])
    return output


def slamf8_high_analysis(adata: ad.AnnData, dataset: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    adata = adata[valid_lineage_mask(adata.obs)].copy()
    sample_col = sample_column(adata.obs)
    if "SLAMF8" not in adata.raw.var_names:
        raise ValueError("SLAMF8 absent from normalized expression")
    genes = [gene for gene in ["SLAMF8", "TREM2", "CD9", "GPNMB", "SPP1", "OLR1", "IL1B", "LGALS3", "HMOX1", "FTH1", "FTL", "PTGS2", "ACSL4"] if gene in adata.raw.var_names]
    expression = adata.raw[:, genes].X
    if sparse.issparse(expression):
        expression = expression.toarray()
    expression = pd.DataFrame(expression, index=adata.obs_names, columns=genes)
    cell_values = adata.obs[[sample_col, "subtype"] + [column for column in adata.obs.columns if column.startswith("module_")]].copy()
    for gene in genes:
        cell_values[f"gene_{gene}"] = expression[gene].to_numpy()
    cell_values["SLAMF8_expression"] = expression["SLAMF8"].to_numpy()

    rows = []
    minimum_high = 5 if dataset == "GSE212837" else 20
    feature_columns = [column for column in cell_values.columns if column.startswith("module_") or column.startswith("gene_")]
    for sample, group in cell_values.groupby(sample_col):
        threshold = float(group["SLAMF8_expression"].quantile(0.75))
        group = group.copy()
        group["SLAMF8_state"] = np.where(group["SLAMF8_expression"] > threshold, "high", "low")
        counts = group["SLAMF8_state"].value_counts()
        if counts.get("high", 0) < minimum_high or counts.get("low", 0) < 20:
            continue
        state_means = group.groupby("SLAMF8_state")[feature_columns].mean()
        for feature in feature_columns:
            rows.append(
                {
                    "dataset": dataset,
                    sample_col: sample,
                    "feature": feature,
                    "n_high": int(counts["high"]),
                    "n_low": int(counts["low"]),
                    "mean_high": state_means.loc["high", feature],
                    "mean_low": state_means.loc["low", feature],
                    "paired_difference": state_means.loc["high", feature] - state_means.loc["low", feature],
                }
            )
    differences = pd.DataFrame(rows)
    tests = []
    if len(differences):
        for feature, group in differences.groupby("feature"):
            delta = group["paired_difference"].dropna().to_numpy(dtype=float)
            if len(delta) < 5 or np.allclose(delta, 0):
                continue
            statistic, pvalue = wilcoxon(delta, alternative="two-sided", zero_method="wilcox")
            tests.append(
                {
                    "dataset": dataset,
                    "feature": feature,
                    "n_samples": len(delta),
                    "median_paired_difference": np.median(delta),
                    "positive_fraction": np.mean(delta > 0),
                    "wilcoxon_statistic": statistic,
                    "p_value": pvalue,
                }
            )
    tests = pd.DataFrame(tests)
    if len(tests):
        tests["fdr"] = bh_adjust(tests["p_value"])
    return differences, tests


def save_figures(dataset: str, tests: pd.DataFrame, high_tests: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="notebook")
    if len(tests):
        selected = tests.loc[tests["feature"].isin(["myeloid_fraction", "SLAMF8_logCPM", "module_FARG95", "module_Iron_homeostasis", "module_Lipid_peroxidation", "module_Senescence_SASP", "module_Inflammation"])]
        if len(selected):
            heatmap = selected.pivot(index="subtype", columns="feature", values="rank_biserial")
            fig, ax = plt.subplots(figsize=(11, max(4, 0.45 * len(heatmap))))
            sns.heatmap(heatmap, cmap="vlag", center=0, vmin=-1, vmax=1, annot=True, fmt=".2f", ax=ax)
            ax.set_title(f"{dataset}: patient-level disease effects")
            fig.tight_layout()
            fig.savefig(FIGURES / f"006_{dataset}_patient_level_subtype_effects.png", dpi=240)
            plt.close(fig)
    if len(high_tests):
        selected = high_tests.sort_values("median_paired_difference")
        fig, ax = plt.subplots(figsize=(9, max(5, 0.3 * len(selected))))
        colors = ["#B2182B" if value > 0 else "#2166AC" for value in selected["median_paired_difference"]]
        ax.barh(selected["feature"], selected["median_paired_difference"], color=colors)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_title(f"{dataset}: SLAMF8-high vs low, paired by patient/sample")
        ax.set_xlabel("Median paired difference")
        fig.tight_layout()
        fig.savefig(FIGURES / f"007_{dataset}_SLAMF8_high_paired_effects.png", dpi=240)
        plt.close(fig)


def run_dataset(dataset: str) -> None:
    adata = ad.read_h5ad(PROCESSED / f"{dataset}_myeloid_reclustered.h5ad")
    values = build_sample_values(adata, dataset)
    tests = condition_tests(values, dataset)
    differences, high_tests = slamf8_high_analysis(adata, dataset)
    values.to_csv(RESULTS / f"013_{dataset}_patient_subtype_values.tsv.gz", sep="\t", index=False, compression="gzip")
    tests.to_csv(RESULTS / f"014_{dataset}_patient_subtype_tests.tsv", sep="\t", index=False)
    differences.to_csv(RESULTS / f"015_{dataset}_SLAMF8_high_paired_differences.tsv.gz", sep="\t", index=False, compression="gzip")
    high_tests.to_csv(RESULTS / f"016_{dataset}_SLAMF8_high_paired_tests.tsv", sep="\t", index=False)
    save_figures(dataset, tests, high_tests)


def main() -> None:
    datasets = sys.argv[1:] or list(COMPARISONS)
    for dataset in datasets:
        run_dataset(dataset)


if __name__ == "__main__":
    main()
