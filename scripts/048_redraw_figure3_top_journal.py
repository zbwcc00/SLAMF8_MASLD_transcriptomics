#!/usr/bin/env python
"""Top-journal redraw of Monocle3/Slingshot trajectory outputs."""

from __future__ import annotations

import gzip
import os
from pathlib import Path

ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "cache" / "matplotlib"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(ROOT / "cache" / "numba"))

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import io, sparse
from scipy.stats import spearmanr


RESULTS = ROOT / "results"
DATA = ROOT / "data_processed"
FIG = ROOT / "figures"
REVISION_FIG = FIG / "revision_p0"
SOURCE_DATA = ROOT / "source_data"
REVISION_FIG.mkdir(parents=True, exist_ok=True)
SOURCE_DATA.mkdir(parents=True, exist_ok=True)
DATASETS = ["GSE344087", "GSE298719", "GSE212837"]
COHORT_COLORS = {"GSE344087": "#0072B2", "GSE298719": "#D55E00", "GSE212837": "#009E73"}
GENES = ["FCN1", "VCAN", "TREM2", "CD9", "GPNMB", "SPP1", "SLAMF8", "HMOX1", "LGALS3", "FTH1", "FTL", "OLR1", "IL1B"]


def load_trajectory(dataset: str) -> pd.DataFrame:
    umap = pd.read_csv(DATA / f"{dataset}_trajectory_umap.tsv.gz", sep="\t")
    metadata = pd.read_csv(DATA / f"{dataset}_trajectory_cell_metadata.tsv.gz", sep="\t")
    monocle = pd.read_csv(RESULTS / f"007_{dataset}_monocle3_pseudotime.tsv.gz", sep="\t")
    slingshot = pd.read_csv(RESULTS / f"008_{dataset}_slingshot_pseudotime.tsv.gz", sep="\t")
    table = umap.merge(metadata[["barcode", "subtype"]], on="barcode", how="left")
    table = table.merge(monocle, on="barcode", how="left").merge(slingshot, on="barcode", how="left")
    return table


def finite_pseudotime(table: pd.DataFrame) -> pd.Series:
    return np.isfinite(table["monocle3_pseudotime"].to_numpy())


def subtype_color(subtype: str) -> str:
    value = str(subtype).lower()
    if "contaminant" in value or "low_quality" in value or "ambiguous" in value:
        return "#BDBDBD"
    if "monocyte" in value or "fcn1" in value or "fcgr3a" in value:
        return "#E69F00"
    if "kupffer" in value or "marco" in value:
        return "#0072B2"
    if "trem2" in value or "spp1" in value or "gpnmb" in value or "macrophage" in value or "samac" in value:
        return "#D55E00"
    if "dendritic" in value or "cdc" in value or "pdc" in value or "lamp3" in value:
        return "#009E73"
    if "cycling" in value:
        return "#CC79A7"
    return "#666666"


def clean_axis(axis):
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_facecolor("white")
    for spine in axis.spines.values():
        spine.set_visible(False)


def plot_monocle_maps(axes):
    for axis, dataset in zip(axes, DATASETS):
        table = load_trajectory(dataset)
        finite = finite_pseudotime(table)
        values = table["monocle3_pseudotime"].to_numpy(dtype=float)
        axis.scatter(table.loc[~finite, "umap_1"], table.loc[~finite, "umap_2"], s=2.0, c="#D5D8DC", alpha=0.45,
                     linewidths=0, rasterized=True)
        if finite.any():
            plot = axis.scatter(table.loc[finite, "umap_1"], table.loc[finite, "umap_2"], c=values[finite],
                                cmap="viridis", s=2.4, alpha=0.75, linewidths=0, rasterized=True)
        else:
            plot = None
        edge_path = RESULTS / f"037_{dataset}_monocle3_graph_edges.tsv"
        if edge_path.exists():
            edges = pd.read_csv(edge_path, sep="\t")
            for _, edge in edges.iterrows():
                axis.plot([edge["x1"], edge["x2"]], [edge["y1"], edge["y2"]], color="#17202A", linewidth=0.55,
                          alpha=0.85, zorder=4)
        clean_axis(axis)
        coverage = 100 * finite.mean()
        axis.set_title(dataset, fontsize=10, fontweight="bold", pad=4)
        axis.text(0.02, 0.04, f"finite pseudotime: {coverage:.1f}%\nroot: FCN1/VCAN state",
                  transform=axis.transAxes, fontsize=7.5, color="#34495E")
        if plot is not None:
            axis.figure.colorbar(plot, ax=axis, fraction=0.035, pad=0.02, label="Monocle3 pseudotime")


def plot_slingshot_maps(axes):
    for axis, dataset in zip(axes, DATASETS):
        table = load_trajectory(dataset)
        column = "Lineage1"
        values = table[column].to_numpy(dtype=float)
        finite = np.isfinite(values)
        normalized = values.copy()
        if finite.any():
            lo, hi = np.nanpercentile(values[finite], [1, 99])
            normalized[finite] = np.clip((values[finite] - lo) / max(hi - lo, 1e-9), 0, 1)
            plot = axis.scatter(table.loc[finite, "umap_1"], table.loc[finite, "umap_2"], c=normalized[finite],
                                cmap="plasma", s=2.4, alpha=0.75, linewidths=0, rasterized=True)
        else:
            plot = None
        clean_axis(axis)
        n_lineages = len([column_name for column_name in table.columns if column_name.startswith("Lineage")])
        axis.set_title(dataset, fontsize=10, fontweight="bold", pad=4)
        axis.text(0.02, 0.04, f"representative Lineage 1\nSlingshot lineages: {n_lineages}",
                  transform=axis.transAxes, fontsize=7.5, color="#34495E")
        if plot is not None:
            axis.figure.colorbar(plot, ax=axis, fraction=0.035, pad=0.02, label="Slingshot Lineage 1")


def plot_concordance(outer_axis):
    outer_axis.axis("off")
    subgrid = outer_axis.get_subplotspec().subgridspec(1, 3, wspace=0.28)
    for index, dataset in enumerate(DATASETS):
        axis = outer_axis.figure.add_subplot(subgrid[0, index])
        table = load_trajectory(dataset)
        finite = finite_pseudotime(table) & np.isfinite(table["Lineage1"].to_numpy())
        x = table.loc[finite, "monocle3_pseudotime"].to_numpy(dtype=float)
        y = table.loc[finite, "Lineage1"].to_numpy(dtype=float)
        if len(x):
            x = (x - x.min()) / max(x.max() - x.min(), 1e-9)
            y = (y - y.min()) / max(y.max() - y.min(), 1e-9)
            axis.scatter(x, y, s=9, color=COHORT_COLORS[dataset], alpha=0.55, edgecolor="white", linewidth=0.2)
            rho, pvalue = spearmanr(x, y)
            axis.text(0.05, 0.95, f"rho={rho:.2f}\nn={len(x):,}", transform=axis.transAxes, va="top", fontsize=7)
        axis.set_title(dataset, fontsize=8.5, fontweight="bold")
        axis.set_xlabel("Monocle3", fontsize=7)
        if index == 0:
            axis.set_ylabel("Slingshot L1", fontsize=7)
        else:
            axis.set_ylabel("")
        axis.tick_params(labelsize=6.5)
        axis.grid(True, color="#E5E7E9", linewidth=0.5)
    outer_axis.figure.text(0.03, 0.365, "Monocle3–Slingshot concordance", fontsize=10.5, fontweight="bold", ha="left")


def load_expression_bins(dataset: str, n_bins: int = 20) -> tuple[np.ndarray, list[str]]:
    pseudo = pd.read_csv(RESULTS / f"007_{dataset}_monocle3_pseudotime.tsv.gz", sep="\t")
    metadata = pd.read_csv(DATA / f"{dataset}_trajectory_cell_metadata.tsv.gz", sep="\t")
    features = pd.read_csv(DATA / f"{dataset}_trajectory_features.tsv.gz", sep="\t")
    barcodes = metadata["barcode"].astype(str).to_numpy()
    pseudo = pseudo.set_index("barcode").reindex(barcodes)["monocle3_pseudotime"].to_numpy(dtype=float)
    keep_cells = np.isfinite(pseudo)
    counts_path = DATA / f"{dataset}_trajectory_counts.mtx.gz"
    with gzip.open(counts_path, "rb") as handle:
        counts = io.mmread(handle).tocsr()
    gene_names = features["var_name"].astype(str).to_numpy()
    gene_to_idx = {gene: index for index, gene in enumerate(gene_names)}
    selected = [gene for gene in GENES if gene in gene_to_idx]
    positions = np.where(keep_cells)[0]
    order = np.argsort(pseudo[positions])
    positions = positions[order]
    n_bins = min(n_bins, max(5, len(positions) // 20))
    bins = np.array_split(np.arange(len(positions)), n_bins)
    library = np.asarray(counts[:, positions].sum(axis=0)).ravel()
    library[library <= 0] = 1
    matrix = np.zeros((len(selected), len(bins)), dtype=float)
    for gene_index, gene in enumerate(selected):
        row = counts[gene_to_idx[gene], positions].toarray().ravel()
        log_cpm = np.log1p(row / library * 1e6)
        for bin_index, relative in enumerate(bins):
            matrix[gene_index, bin_index] = float(np.mean(log_cpm[relative]))
        row_values = matrix[gene_index]
        matrix[gene_index] = (row_values - row_values.mean()) / max(row_values.std(), 1e-9)
    return matrix, selected


def plot_dynamic_heatmap(axis):
    pieces = []
    selected_genes = None
    for dataset in DATASETS:
        matrix, genes = load_expression_bins(dataset)
        selected_genes = genes if selected_genes is None else [gene for gene in selected_genes if gene in genes]
        pieces.append((dataset, matrix, genes))
    common_genes = [gene for gene in GENES if all(gene in genes for _, _, genes in pieces)]
    blocks = []
    for dataset, matrix, genes in pieces:
        indices = [genes.index(gene) for gene in common_genes]
        blocks.append(matrix[indices])
    heat = np.concatenate(blocks, axis=1)
    image = axis.imshow(heat, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2, interpolation="nearest")
    axis.set_yticks(np.arange(len(common_genes)), common_genes, fontsize=7.5)
    boundaries = np.cumsum([block.shape[1] for block in blocks])
    centers = [0] + [int((boundaries[i - 1] + boundaries[i]) / 2) for i in range(1, len(boundaries))]
    axis.set_xticks(centers, DATASETS, fontsize=7)
    for boundary in boundaries[:-1]:
        axis.axvline(boundary - 0.5, color="white", linewidth=1.2)
    axis.set_title("Candidate gene dynamics along Monocle3 pseudotime", fontsize=10.5, fontweight="bold", loc="left", pad=6)
    axis.set_xlabel("Early → late pseudotime bins; each cohort standardized separately", fontsize=7.5)
    axis.set_ylabel("")
    axis.figure.colorbar(image, ax=axis, fraction=0.04, pad=0.02, label="Within-cohort z-score")
    axis.text(0.0, -0.18, "Binned expression is descriptive; it does not replace a fitted differential-expression test.",
              transform=axis.transAxes, fontsize=7.5, color="#566573")


def save_main():
    sns.set_theme(style="white", context="paper", font="DejaVu Sans")
    fig = plt.figure(figsize=(16.0, 13.5), dpi=300)
    grid = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 1.25], hspace=0.36)
    top = grid[0].subgridspec(1, 3, wspace=0.12)
    middle = grid[1].subgridspec(1, 3, wspace=0.12)
    top_axes = [fig.add_subplot(top[0, i]) for i in range(3)]
    middle_axes = [fig.add_subplot(middle[0, i]) for i in range(3)]
    plot_monocle_maps(top_axes)
    plot_slingshot_maps(middle_axes)
    bottom = grid[2].subgridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.24)
    concordance_axis = fig.add_subplot(bottom[0, 0])
    heatmap_axis = fig.add_subplot(bottom[0, 1])
    plot_concordance(concordance_axis)
    plot_dynamic_heatmap(heatmap_axis)
    fig.text(0.015, 0.675, "Monocle3", rotation=90, rotation_mode="anchor", va="center", ha="left", fontsize=10.5, fontweight="bold")
    fig.text(0.015, 0.405, "Slingshot", rotation=90, rotation_mode="anchor", va="center", ha="left", fontsize=10.5, fontweight="bold")
    fig.text(0.015, 0.105, "Cross-method audit", rotation=90, rotation_mode="anchor", va="center", ha="left", fontsize=10.5, fontweight="bold")
    for label, x, y in [("A", 0.025, 0.985), ("B", 0.025, 0.66), ("C", 0.025, 0.34)]:
        fig.text(x, y, label, fontsize=14, fontweight="bold", va="top")
    fig.suptitle("Macrophage trajectory inference across three liver single-cell cohorts", fontsize=15, fontweight="bold", y=0.995)
    fig.text(0.03, 0.012, "Monocle3 disconnected cells are shown in gray; trajectory inference is state ordering, not lineage tracing or temporal proof.",
             fontsize=8.5, color="#566573")
    fig.savefig(FIG / "Figure3_macrophage_trajectories_v2.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / "Figure3_macrophage_trajectories_v2.png", dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def save_audit_summary():
    rows = []
    for dataset in DATASETS:
        table = load_trajectory(dataset)
        finite = finite_pseudotime(table)
        lineages = [column for column in table.columns if column.startswith("Lineage")]
        correlations = []
        for lineage in lineages:
            keep = finite & np.isfinite(table[lineage].to_numpy())
            if keep.sum() >= 20:
                correlations.append(spearmanr(table.loc[keep, "monocle3_pseudotime"], table.loc[keep, lineage]).statistic)
        rows.append({"dataset": dataset, "n_cells": len(table), "monocle3_finite_fraction": finite.mean(),
                     "n_slingshot_lineages": len(lineages), "median_monocle3_slingshot_rho": np.nanmedian(correlations)})
    summary = pd.DataFrame(rows)
    summary.to_csv(RESULTS / "036_Figure3_trajectory_audit_summary.tsv", sep="\t", index=False)
    sns.set_theme(style="whitegrid", context="paper", font="DejaVu Sans")
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.3), dpi=300, constrained_layout=True)
    axes[0].bar(summary["dataset"], summary["monocle3_finite_fraction"], color=[COHORT_COLORS[d] for d in DATASETS])
    axes[0].set_ylim(0, 1)
    axes[0].set_ylabel("Fraction finite")
    axes[0].set_title("Monocle3 finite pseudotime", fontweight="bold")
    axes[1].bar(summary["dataset"], summary["n_slingshot_lineages"], color=[COHORT_COLORS[d] for d in DATASETS])
    axes[1].set_ylabel("Number of inferred lineages")
    axes[1].set_title("Slingshot branching output", fontweight="bold")
    axes[2].bar(summary["dataset"], summary["median_monocle3_slingshot_rho"], color=[COHORT_COLORS[d] for d in DATASETS])
    axes[2].set_ylim(0, 1)
    axes[2].set_ylabel("Median Spearman rho")
    axes[2].set_title("Cross-method concordance", fontweight="bold")
    for axis in axes:
        axis.tick_params(axis="x", rotation=0, labelsize=8)
        axis.grid(axis="x", visible=False)
    fig.suptitle("Trajectory robustness audit", fontsize=12, fontweight="bold")
    stems = [FIG / "Supplementary_Figure3_trajectory_audit",
             REVISION_FIG / "Supplementary_Figure5_trajectory_audit"]
    for stem in stems:
        fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
        fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
        fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")
        fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    summary.to_csv(SOURCE_DATA / "Source_Data_Supplementary_Figure5_trajectory_audit.tsv.gz",
                  sep="\t", index=False, compression="gzip")
    plt.close(fig)


if __name__ == "__main__":
    save_main()
    save_audit_summary()
    print("Figure 3 v2 and trajectory audit written to", FIG)
