#!/usr/bin/env python
"""Top-journal redraw of Figure 2 and missing state-association panels."""

from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "cache" / "matplotlib"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(ROOT / "cache" / "numba"))

import anndata  # noqa: F401
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns


RESULTS = ROOT / "results"
FIG = ROOT / "figures"
DATASETS = ["GSE344087", "GSE298719", "GSE212837"]
DATASET_COLORS = {"GSE344087": "#0072B2", "GSE298719": "#D55E00", "GSE212837": "#009E73"}
FEATURE_LABELS = {
    "gene_SLAMF8": "SLAMF8",
    "gene_GPNMB": "GPNMB",
    "gene_TREM2": "TREM2",
    "gene_CD9": "CD9",
    "gene_SPP1": "SPP1",
    "module_FARG95": "FARG95",
    "module_Iron_homeostasis": "Iron homeostasis",
    "module_Senescence_SASP": "Senescence/SASP",
    "module_Inflammation": "Inflammation",
    "module_Lipid_peroxidation": "Lipid peroxidation",
    "module_HSC_activation": "HSC activation",
    "module_PUFA_ACSL4": "PUFA–ACSL4",
}


def star(value: float) -> str:
    if not np.isfinite(value):
        return ""
    if value < 0.001:
        return "***"
    if value < 0.01:
        return "**"
    if value < 0.05:
        return "*"
    return ""


def read_state_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    combined = pd.read_csv(RESULTS / "021_cross_dataset_SLAMF8_high_tests.tsv", sep="\t")
    meta = pd.read_csv(RESULTS / "022_cross_dataset_SLAMF8_high_meta.tsv", sep="\t")
    return combined, meta


def read_disease_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    combined = pd.read_csv(RESULTS / "023_cross_dataset_disease_macrophage_tests.tsv", sep="\t")
    meta = pd.read_csv(RESULTS / "024_cross_dataset_disease_macrophage_meta.tsv", sep="\t")
    return combined, meta


def styled_heatmap(axis, values: pd.DataFrame, fdr: pd.DataFrame, title: str, xlabel: str):
    finite = values.to_numpy(dtype=float)
    limit = max(0.05, float(np.nanpercentile(np.abs(finite), 97)))
    sns.heatmap(values, ax=axis, cmap="vlag", center=0, vmin=-limit, vmax=limit,
                linewidths=0, cbar_kws={"pad": 0.12},
                annot=False, xticklabels=True, yticklabels=True)
    axis.grid(False)
    colorbar = axis.collections[0].colorbar
    colorbar.set_label("")
    colorbar.set_ticks([0.0])
    for i, row in enumerate(values.index):
        for j, col in enumerate(values.columns):
            value = values.loc[row, col]
            q = fdr.loc[row, col] if row in fdr.index and col in fdr.columns else np.nan
            if np.isfinite(value):
                axis.text(j + 0.5, i + 0.5, f"{value:.2f}", ha="center", va="center", fontsize=7,
                          color="#17202A")
    axis.set_title(title, fontsize=8.5, fontweight="bold", pad=4, loc="left")
    axis.set_xlabel("")
    axis.set_ylabel("")
    axis.tick_params(axis="x", rotation=0, labelsize=7)
    axis.tick_params(axis="y", rotation=0, labelsize=7)


def state_heatmap(axis):
    combined, meta = read_state_tables()
    features = ["gene_SLAMF8", "module_FARG95", "module_Iron_homeostasis", "module_Senescence_SASP",
                "module_Inflammation", "gene_GPNMB", "gene_TREM2", "gene_CD9", "gene_SPP1"]
    effects = combined.loc[combined["feature"].isin(features)].pivot(index="feature", columns="dataset", values="median_paired_difference")
    qvalues = combined.loc[combined["feature"].isin(features)].pivot(index="feature", columns="dataset", values="fdr")
    effects = effects.reindex(index=features, columns=DATASETS)
    qvalues = qvalues.reindex(index=features, columns=DATASETS)
    compact_names = {
        "gene_SLAMF8": "SLAMF8",
        "module_FARG95": "FARG95",
        "module_Iron_homeostasis": "Iron",
        "module_Senescence_SASP": "SASP",
        "module_Inflammation": "Inflam.",
        "gene_GPNMB": "GPNMB",
        "gene_TREM2": "TREM2",
        "gene_CD9": "CD9",
        "gene_SPP1": "SPP1",
    }
    labels = []
    for feature in features:
        row = meta.loc[meta["feature"] == feature]
        q = float(row["signed_stouffer_fdr"].iloc[0]) if len(row) else np.nan
        labels.append(f"{compact_names[feature]} [q={q:.2g}]")
    effects.index = labels
    qvalues.index = labels
    styled_heatmap(axis, effects, qvalues, "SLAMF8-high versus low macrophage state", "Difference")
    axis.set_xticklabels(["GSE344087", "GSE298719", "GSE212837"], fontsize=7)


def disease_heatmap(axis):
    combined, meta = read_disease_tables()
    features = ["SLAMF8_logCPM", "module_FARG95", "module_Iron_homeostasis", "module_Senescence_SASP",
                "module_Inflammation", "module_HSC_activation"]
    effects = combined.loc[combined["feature"].isin(features)].pivot(index="feature", columns="dataset", values="median_difference")
    qvalues = combined.loc[combined["feature"].isin(features)].pivot(index="feature", columns="dataset", values="fdr")
    effects = effects.reindex(index=features, columns=DATASETS)
    qvalues = qvalues.reindex(index=features, columns=DATASETS)
    compact_names = {
        "SLAMF8_logCPM": "SLAMF8",
        "module_FARG95": "FARG95",
        "module_Iron_homeostasis": "Iron",
        "module_Senescence_SASP": "SASP",
        "module_Inflammation": "Inflam.",
        "module_HSC_activation": "HSC activation",
    }
    labels = []
    for feature in features:
        row = meta.loc[meta["feature"] == feature]
        q = float(row["signed_stouffer_fdr"].iloc[0]) if len(row) else np.nan
        labels.append(compact_names[feature])
    effects.index = labels
    qvalues.index = labels
    styled_heatmap(axis, effects, qvalues, "Disease/fibrosis differences", "Difference")
    axis.set_xticklabels(["GSE344087\nfibrosis", "GSE298719\nMASLD", "GSE212837\nNASH"], fontsize=6.3)


def hsc_scatter(outer_axis):
    outer_axis.axis("off")
    axes = outer_axis.get_subplotspec().subgridspec(1, 3, wspace=0.42)
    corr = pd.read_csv(RESULTS / "033_cross_dataset_macrophage_HSC_correlations.tsv", sep="\t")
    for index, dataset in enumerate(DATASETS):
        axis = outer_axis.figure.add_subplot(axes[0, index])
        values = pd.read_csv(RESULTS / f"031_{dataset}_macrophage_HSC_matched_values.tsv", sep="\t")
        x = values["SLAMF8_logCPM"].astype(float)
        y = values["module_HSC_activation_HSC"].astype(float)
        keep = x.notna() & y.notna()
        axis.scatter(x[keep], y[keep], s=24, color=DATASET_COLORS[dataset], alpha=0.82,
                     edgecolor="white", linewidth=0.35)
        if keep.sum() >= 3:
            coefficient = np.polyfit(x[keep], y[keep], 1)
            grid = np.linspace(float(x[keep].min()), float(x[keep].max()), 60)
            axis.plot(grid, coefficient[0] * grid + coefficient[1], color="#34495E", linewidth=1.4)
        axis.set_title(dataset, fontsize=7.5, fontweight="bold", pad=3)
        axis.set_xlabel("SLAMF8\n(logCPM)", fontsize=6.5)
        if index == 0:
            axis.set_ylabel("HSC activation", fontsize=6.5)
        else:
            axis.set_ylabel("")
        axis.tick_params(axis="x", labelsize=6.5)
        axis.tick_params(axis="y", labelsize=6.5, pad=4)
        axis.grid(False)
        for spine in axis.spines.values():
            spine.set_color("#AAB7B8")


def receptor_summary(axis):
    meta = pd.read_csv(RESULTS / "034_cross_dataset_macrophage_HSC_correlation_meta.tsv", sep="\t")
    features = ["logCPM_COL1A2", "logCPM_PDGFRA", "logCPM_PDGFRB", "logCPM_COL1A1", "logCPM_COL3A1",
                "module_HSC_activation", "logCPM_TGFB1", "logCPM_ACTA2"]
    table = meta.loc[meta["feature"].isin(features)].copy()
    table["label"] = table["feature"].map({
        "logCPM_COL1A2": "COL1A2", "logCPM_PDGFRA": "PDGFRA", "logCPM_PDGFRB": "PDGFRB",
        "logCPM_COL1A1": "COL1A1", "logCPM_COL3A1": "COL3A1", "module_HSC_activation": "HSC score",
        "logCPM_TGFB1": "TGFB1", "logCPM_ACTA2": "ACTA2",
    })
    table = table.set_index("label").reindex(["COL1A2", "PDGFRA", "PDGFRB", "COL1A1", "COL3A1", "HSC score", "TGFB1", "ACTA2"]).dropna(subset=["pooled_spearman_rho"])
    y = np.arange(len(table))
    colors = ["#D55E00" if q < 0.05 else "#7F8C8D" for q in table["fdr"]]
    axis.axvline(0, color="#17202A", linewidth=0.8)
    axis.hlines(y, 0, table["pooled_spearman_rho"], color=colors, linewidth=2.2)
    axis.scatter(table["pooled_spearman_rho"], y, s=48, color=colors, edgecolor="white", linewidth=0.5, zorder=3)
    axis.set_yticks(y, table.index, fontsize=7)
    axis.set_xlabel("Cross-cohort pooled Spearman rho", fontsize=7)
    axis.set_title("HSC receptor/program association", fontsize=8.5, fontweight="bold", loc="left", pad=4)
    axis.grid(False)
    axis.tick_params(axis="x", labelsize=7)
    for y_position, (_, row) in zip(y, table.iterrows()):
        rho = float(row["pooled_spearman_rho"])
        if rho > 0.36:
            axis.text(rho - 0.015, y_position + 0.08, f"q={row['fdr']:.3f}",
                      ha="right", va="bottom", fontsize=6.3, clip_on=False)
        elif y_position == y.max():
            axis.text(rho + 0.015, y_position - 0.08, f"q={row['fdr']:.3f}",
                      ha="left", va="top", fontsize=6.3, clip_on=False)
        else:
            axis.text(rho + 0.015, y_position + 0.08, f"q={row['fdr']:.3f}",
                      ha="left", va="bottom", fontsize=6.3, clip_on=False)


def bootstrap_median(values: np.ndarray, seed: int = 20260920, n_boot: int = 2000) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    samples = rng.choice(values, size=(n_boot, len(values)), replace=True)
    medians = np.median(samples, axis=1)
    return float(np.median(values)), float(np.quantile(medians, 0.025)), float(np.quantile(medians, 0.975))


def save_bootstrap_supplement():
    combined, _ = read_state_tables()
    features = ["gene_SLAMF8", "module_FARG95", "module_Iron_homeostasis", "module_Senescence_SASP",
                "module_Inflammation", "gene_GPNMB", "gene_TREM2", "gene_CD9", "gene_SPP1"]
    records = []
    for feature in features:
        for dataset in DATASETS:
            path = RESULTS / f"015_{dataset}_SLAMF8_high_paired_differences.tsv.gz"
            values = pd.read_csv(path, sep="\t")
            values = values.loc[values["feature"] == feature, "paired_difference"].to_numpy(dtype=float)
            median, lower, upper = bootstrap_median(values, seed=20260920 + features.index(feature) * 19 + DATASETS.index(dataset))
            row = combined.loc[(combined["dataset"] == dataset) & (combined["feature"] == feature)]
            records.append({"feature": feature, "dataset": dataset, "median": median, "lower": lower, "upper": upper,
                            "fdr": float(row["fdr"].iloc[0]) if len(row) else np.nan})
    table = pd.DataFrame(records)
    table.to_csv(RESULTS / "035_Figure2_state_bootstrap_ci.tsv", sep="\t", index=False)
    sns.set_theme(style="whitegrid", context="paper", font="DejaVu Sans")
    fig, axis = plt.subplots(figsize=(7.2, 5.6), dpi=300)
    y_base = np.arange(len(features))[::-1]
    offsets = {"GSE344087": -0.18, "GSE298719": 0.0, "GSE212837": 0.18}
    for dataset in DATASETS:
        subset = table[table["dataset"] == dataset].set_index("feature").reindex(features)
        y = y_base + offsets[dataset]
        axis.errorbar(subset["median"], y,
                      xerr=[subset["median"] - subset["lower"], subset["upper"] - subset["median"]],
                      fmt="o", color=DATASET_COLORS[dataset], markersize=5, capsize=2.5, linewidth=1.2,
                      label=dataset)
    axis.axvline(0, color="#17202A", linewidth=0.8)
    axis.set_yticks(y_base, [FEATURE_LABELS[f] for f in features], fontsize=7)
    axis.set_xlabel("Median paired difference (SLAMF8-high minus low)", fontsize=7.5)
    axis.set_ylabel("")
    axis.set_title("Bootstrap uncertainty of within-unit SLAMF8 state effects", fontsize=9.5, fontweight="bold", loc="left")
    if axis.get_legend() is not None:
        axis.get_legend().remove()
    axis.grid(axis="y", visible=False)
    axis.text(0.0, -0.11, "Points are medians; bars are deterministic percentile bootstrap 95% intervals within cohort.",
              transform=axis.transAxes, fontsize=6.5, color="#566573")
    fig.tight_layout()
    fig.savefig(FIG / "Supplementary_Figure2_SLAMF8_state_bootstrap.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / "Supplementary_Figure2_SLAMF8_state_bootstrap.svg", bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / "Supplementary_Figure2_SLAMF8_state_bootstrap.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / "Supplementary_Figure2_SLAMF8_state_bootstrap.tiff", dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def save_main():
    sns.set_theme(style="whitegrid", context="paper", font="DejaVu Sans")
    fig = plt.figure(figsize=(7.2, 5.8), dpi=300)
    grid = fig.add_gridspec(2, 2, width_ratios=[1.15, 1.0], height_ratios=[1.0, 1.02], wspace=0.26, hspace=0.36)
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]
    state_heatmap(axes[0])
    disease_heatmap(axes[1])
    hsc_scatter(axes[2])
    receptor_summary(axes[3])
    for label, axis in zip("ABCD", axes):
        axis.text(-0.10, 1.08, label, transform=axis.transAxes, fontsize=8, fontweight="bold", va="top", clip_on=False)
    fig.suptitle("SLAMF8-high macrophage state and matched HSC activation in MASLD/NASH\n"
                 "Cohort colors: GSE344087 blue · GSE298719 orange · GSE212837 green",
                 fontsize=9.2, fontweight="bold", y=0.985)
    fig.text(0.02, 0.012, "Effects are patient/sample-level associations; no panel establishes causality.",
             fontsize=6.5, color="#566573")
    fig.savefig(FIG / "Figure2_SLAMF8_macrophage_state_v2.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / "Figure2_SLAMF8_macrophage_state_v2.svg", bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / "Figure2_SLAMF8_macrophage_state_v2.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / "Figure2_SLAMF8_macrophage_state_v2.tiff", dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    save_main()
    save_bootstrap_supplement()
    print("Figure 2 v2 and bootstrap supplement written to", FIG)
