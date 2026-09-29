#!/usr/bin/env python
"""Submission redraw of Monocle3 and Slingshot macrophage trajectories."""

from __future__ import annotations
import os

import gzip
from pathlib import Path

import matplotlib

ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
FIG = ROOT / "figures" / "revision_p0"
RESULTS = ROOT / "results"
DATA = ROOT / "data_processed"
SOURCE_DATA = ROOT / "source_data"
FIG.mkdir(parents=True, exist_ok=True)
SOURCE_DATA.mkdir(parents=True, exist_ok=True)

matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.linewidth": 0.7,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from scipy import io
from scipy.stats import spearmanr


DATASETS = ["GSE344087", "GSE298719", "GSE212837"]
COHORT_COLORS = {"GSE344087": "#0072B2", "GSE298719": "#D55E00", "GSE212837": "#009E73"}
GENES = ["FCN1", "VCAN", "TREM2", "CD9", "GPNMB", "SPP1", "SLAMF8", "HMOX1", "LGALS3", "FTH1", "FTL", "OLR1", "IL1B"]


def load_trajectory(dataset: str) -> pd.DataFrame:
    umap = pd.read_csv(DATA / f"{dataset}_trajectory_umap.tsv.gz", sep="\t")
    metadata = pd.read_csv(DATA / f"{dataset}_trajectory_cell_metadata.tsv.gz", sep="\t")
    monocle = pd.read_csv(RESULTS / f"007_{dataset}_monocle3_pseudotime.tsv.gz", sep="\t")
    slingshot = pd.read_csv(RESULTS / f"008_{dataset}_slingshot_pseudotime.tsv.gz", sep="\t")
    table = umap.merge(metadata[["barcode", "subtype"]], on="barcode", how="left")
    return table.merge(monocle, on="barcode", how="left").merge(slingshot, on="barcode", how="left")


def finite_pseudotime(table: pd.DataFrame) -> np.ndarray:
    return np.isfinite(table["monocle3_pseudotime"].to_numpy(dtype=float))


def normalize_order(values: np.ndarray, robust: bool = False) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    output = np.full(values.shape, np.nan, dtype=float)
    if not finite.any():
        return output
    if robust:
        lower, upper = np.nanpercentile(values[finite], [1, 99])
    else:
        lower, upper = np.nanmin(values[finite]), np.nanmax(values[finite])
    output[finite] = np.clip((values[finite] - lower) / max(upper - lower, 1e-12), 0, 1)
    return output


def clean_map_axis(axis: plt.Axes) -> None:
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_facecolor("white")
    for spine in axis.spines.values():
        spine.set_visible(False)


def add_panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(-0.08, 1.12, label, transform=axis.transAxes, fontsize=8.5,
              fontweight="bold", va="top", ha="left", clip_on=False)


def plot_monocle_maps(fig: plt.Figure, axes: list[plt.Axes]) -> None:
    mappable = None
    norm = Normalize(0, 1)
    for axis, dataset in zip(axes, DATASETS):
        table = load_trajectory(dataset)
        finite = finite_pseudotime(table)
        values = normalize_order(table["monocle3_pseudotime"].to_numpy(dtype=float))
        axis.scatter(table.loc[~finite, "umap_1"], table.loc[~finite, "umap_2"], s=1.7,
                     c="#D5D8DC", alpha=0.48, linewidths=0, rasterized=True)
        if finite.any():
            mappable = axis.scatter(table.loc[finite, "umap_1"], table.loc[finite, "umap_2"],
                                    c=values[finite], cmap="viridis", norm=norm, s=2.0,
                                    alpha=0.76, linewidths=0, rasterized=True)
        edge_path = RESULTS / f"037_{dataset}_monocle3_graph_edges.tsv"
        if edge_path.exists():
            edges = pd.read_csv(edge_path, sep="\t")
            for _, edge in edges.iterrows():
                axis.plot([edge["x1"], edge["x2"]], [edge["y1"], edge["y2"]],
                          color="#17202A", linewidth=0.45, alpha=0.75, zorder=4)
        clean_map_axis(axis)
        coverage = 100 * finite.mean()
        axis.set_title(f"{dataset}\nfinite {coverage:.1f}% · root FCN1/VCAN",
                       fontsize=7.2, fontweight="bold", pad=3)
    if mappable is not None:
        colorbar = fig.colorbar(mappable, ax=axes, fraction=0.018, pad=0.015)
        colorbar.set_ticks([0, 1])
        colorbar.set_ticklabels(["early", "late"])
        colorbar.ax.tick_params(labelsize=6.2, length=2)
        colorbar.set_label("Within-cohort order", fontsize=6.2, labelpad=3)


def plot_slingshot_maps(fig: plt.Figure, axes: list[plt.Axes]) -> None:
    mappable = None
    norm = Normalize(0, 1)
    for axis, dataset in zip(axes, DATASETS):
        table = load_trajectory(dataset)
        values = table["Lineage1"].to_numpy(dtype=float)
        normalized = normalize_order(values, robust=True)
        finite = np.isfinite(normalized)
        if finite.any():
            mappable = axis.scatter(table.loc[finite, "umap_1"], table.loc[finite, "umap_2"],
                                    c=normalized[finite], cmap="plasma", norm=norm, s=2.0,
                                    alpha=0.76, linewidths=0, rasterized=True)
        clean_map_axis(axis)
        n_lineages = sum(column.startswith("Lineage") for column in table.columns)
        axis.set_title(f"{dataset}\nLineage 1 · {n_lineages} inferred lineages",
                       fontsize=7.2, fontweight="bold", pad=3)
    if mappable is not None:
        colorbar = fig.colorbar(mappable, ax=axes, fraction=0.018, pad=0.015)
        colorbar.set_ticks([0, 1])
        colorbar.set_ticklabels(["early", "late"])
        colorbar.ax.tick_params(labelsize=6.2, length=2)
        colorbar.set_label("Within-cohort order", fontsize=6.2, labelpad=3)


def plot_concordance(outer_axis: plt.Axes) -> None:
    outer_axis.axis("off")
    subgrid = outer_axis.get_subplotspec().subgridspec(1, 3, wspace=0.42)
    for index, dataset in enumerate(DATASETS):
        axis = outer_axis.figure.add_subplot(subgrid[0, index])
        table = load_trajectory(dataset)
        finite = finite_pseudotime(table) & np.isfinite(table["Lineage1"].to_numpy(dtype=float))
        x = table.loc[finite, "monocle3_pseudotime"].to_numpy(dtype=float)
        y = table.loc[finite, "Lineage1"].to_numpy(dtype=float)
        if len(x):
            x = normalize_order(x)
            y = normalize_order(y, robust=True)
            rho = float(spearmanr(x, y).statistic)
            axis.scatter(x, y, s=5.5, color=COHORT_COLORS[dataset], alpha=0.48,
                         edgecolor="none", rasterized=True)
            title = dataset
        else:
            title = f"{dataset}\nno finite pairs"
        axis.set_title(title, fontsize=6.8, fontweight="bold", pad=3)
        axis.set_xlabel("Monocle3 order", fontsize=6.3)
        axis.set_ylabel("Slingshot L1 order" if index == 0 else "", fontsize=6.3)
        axis.tick_params(labelsize=6.0, pad=2)
        axis.grid(False)
        for spine in axis.spines.values():
            spine.set_color("#AAB7B8")
    add_panel_label(outer_axis, "C")


def load_expression_bins(dataset: str, n_bins: int = 20) -> tuple[np.ndarray, list[str]]:
    pseudo = pd.read_csv(RESULTS / f"007_{dataset}_monocle3_pseudotime.tsv.gz", sep="\t")
    metadata = pd.read_csv(DATA / f"{dataset}_trajectory_cell_metadata.tsv.gz", sep="\t")
    features = pd.read_csv(DATA / f"{dataset}_trajectory_features.tsv.gz", sep="\t")
    barcodes = metadata["barcode"].astype(str).to_numpy()
    pseudo_values = pseudo.set_index("barcode").reindex(barcodes)["monocle3_pseudotime"].to_numpy(dtype=float)
    positions = np.where(np.isfinite(pseudo_values))[0]
    positions = positions[np.argsort(pseudo_values[positions])]
    counts_path = DATA / f"{dataset}_trajectory_counts.mtx.gz"
    with gzip.open(counts_path, "rb") as handle:
        counts = io.mmread(handle).tocsr()
    gene_names = features["var_name"].astype(str).to_numpy()
    gene_to_idx = {gene: index for index, gene in enumerate(gene_names)}
    selected = [gene for gene in GENES if gene in gene_to_idx]
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


def plot_dynamic_heatmap(axis: plt.Axes) -> None:
    pieces = [(dataset, *load_expression_bins(dataset)) for dataset in DATASETS]
    common_genes = [gene for gene in GENES if all(gene in genes for _, _, genes in pieces)]
    blocks = [matrix[[genes.index(gene) for gene in common_genes]] for _, matrix, genes in pieces]
    heat = np.concatenate(blocks, axis=1)
    image = axis.imshow(heat, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2,
                        interpolation="nearest", rasterized=True)
    axis.set_yticks(np.arange(len(common_genes)), common_genes, fontsize=6.5)
    boundaries = np.cumsum([block.shape[1] for block in blocks])
    centers = [0] + [int((boundaries[i - 1] + boundaries[i]) / 2) for i in range(1, len(boundaries))]
    axis.set_xticks(centers, DATASETS, fontsize=6.2)
    for boundary in boundaries[:-1]:
        axis.axvline(boundary - 0.5, color="white", linewidth=0.9)
    axis.set_title("Candidate gene dynamics along Monocle3 pseudotime", fontsize=8.0,
                   fontweight="bold", loc="left", pad=3)
    axis.set_xlabel("Early → late bins; each cohort standardized separately", fontsize=6.4, labelpad=3)
    axis.set_ylabel("")
    colorbar = axis.figure.colorbar(image, ax=axis, fraction=0.035, pad=0.018)
    colorbar.set_ticks([-1.5, 0, 1.5])
    colorbar.set_label("Within-cohort z-score", fontsize=6.2, labelpad=3)
    colorbar.ax.tick_params(labelsize=6.0, length=2)
    add_panel_label(axis, "D")


def build_trajectory_source_data() -> None:
    rows = []
    summary_rows = []
    for dataset in DATASETS:
        table = load_trajectory(dataset).copy()
        table.insert(0, "dataset", dataset)
        table.insert(1, "panel", "trajectory_cells")
        table["monocle3_finite"] = np.isfinite(table["monocle3_pseudotime"].to_numpy(dtype=float))
        keep = ["panel", "dataset", "barcode", "umap_1", "umap_2", "subtype",
                "monocle3_pseudotime", "monocle3_finite"]
        keep += [column for column in table.columns if column.startswith("Lineage")]
        rows.append(table[keep])
        finite = finite_pseudotime(table) & np.isfinite(table["Lineage1"].to_numpy(dtype=float))
        rho = float(spearmanr(table.loc[finite, "monocle3_pseudotime"], table.loc[finite, "Lineage1"]).statistic)
        summary_rows.append({"panel": "concordance_summary", "dataset": dataset,
                             "spearman_rho_lineage1": rho, "n_pairs": int(finite.sum()),
                             "monocle3_finite_fraction": float(finite_pseudotime(table).mean()),
                             "n_slingshot_lineages": int(sum(column.startswith("Lineage") for column in table.columns))})
    dynamic_rows = []
    for dataset in DATASETS:
        matrix, genes = load_expression_bins(dataset)
        for gene_index, gene in enumerate(genes):
            for bin_index, value in enumerate(matrix[gene_index], start=1):
                dynamic_rows.append({"panel": "dynamic_bins", "dataset": dataset, "gene": gene,
                                     "bin": bin_index, "z_score": float(value)})
    source = pd.concat(rows, ignore_index=True, sort=False)
    source = pd.concat([source, pd.DataFrame(dynamic_rows)], ignore_index=True, sort=False)
    source = pd.concat([source, pd.DataFrame(summary_rows)], ignore_index=True, sort=False)
    source.to_csv(SOURCE_DATA / "Source_Data_Figure3_trajectory.tsv.gz", sep="\t", index=False, compression="gzip")


def save_main() -> None:
    fig = plt.figure(figsize=(7.2, 8.0), dpi=300)
    grid = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 1.18], hspace=0.48)
    top = grid[0].subgridspec(1, 3, wspace=0.10)
    middle = grid[1].subgridspec(1, 3, wspace=0.10)
    top_axes = [fig.add_subplot(top[0, index]) for index in range(3)]
    middle_axes = [fig.add_subplot(middle[0, index]) for index in range(3)]
    plot_monocle_maps(fig, top_axes)
    plot_slingshot_maps(fig, middle_axes)
    bottom = grid[2].subgridspec(1, 2, width_ratios=[1.0, 1.32], wspace=0.26)
    concordance_axis = fig.add_subplot(bottom[0, 0])
    heatmap_axis = fig.add_subplot(bottom[0, 1])
    plot_concordance(concordance_axis)
    plot_dynamic_heatmap(heatmap_axis)
    add_panel_label(top_axes[0], "A")
    add_panel_label(middle_axes[0], "B")
    fig.suptitle("Cross-method inference of macrophage state ordering", fontsize=10.5,
                 fontweight="bold", y=0.995)
    fig.text(0.03, 0.012,
             "Monocle3-disconnected cells are gray; trajectory inference represents state ordering, not lineage tracing or temporal proof.",
             fontsize=6.4, color="#566573")
    stem = FIG / "Figure3_macrophage_trajectories_v3_corrected"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    save_main()
    build_trajectory_source_data()
    print("Figure 3 submission redraw written to", FIG)
