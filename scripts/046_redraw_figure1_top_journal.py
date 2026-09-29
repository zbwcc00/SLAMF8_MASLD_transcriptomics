#!/usr/bin/env python
"""Redraw the single-cell atlas and produce missing scRNA QC/annotation panels.

All outputs are written to the project directory on D:. The script only uses
the already curated AnnData objects and does not alter biological annotations.
"""

from __future__ import annotations

from pathlib import Path
import os

ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "cache" / "matplotlib"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(ROOT / "cache" / "numba"))

import anndata as ad
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import pandas as pd
import seaborn as sns


DATA = ROOT / "data_processed"
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
SOURCE_DATA = ROOT / "source_data"
REVISION_FIG = FIG / "revision_p0"
REVISION_FIG.mkdir(parents=True, exist_ok=True)
SOURCE_DATA.mkdir(parents=True, exist_ok=True)

DATASETS = ["GSE344087", "GSE298719", "GSE212837"]
LINEAGES = ["myeloid", "hsc"]
DISPLAY = {"GSE344087": "GSE344087", "GSE298719": "GSE298719", "GSE212837": "GSE212837"}

MYELOID_COLORS = {
    "Kupffer": "#0072B2",
    "Macrophage": "#D55E00",
    "Monocyte": "#E69F00",
    "Dendritic": "#009E73",
    "Cycling": "#CC79A7",
    "Contaminant": "#BDBDBD",
    "Other myeloid": "#666666",
}
HSC_COLORS = {
    "Quiescent HSC": "#0072B2",
    "Activated HSC": "#D55E00",
    "Myofibroblast": "#E69F00",
    "Portal fibroblast": "#009E73",
    "Contaminant": "#BDBDBD",
    "Other HSC": "#666666",
}

MARKERS = {
    "myeloid": {
        "Kupffer": ["C1QA", "MARCO", "VSIG4"],
        "Macrophage": ["TREM2", "SPP1", "GPNMB"],
        "Monocyte": ["FCN1", "VCAN", "S100A8"],
        "Dendritic": ["FCER1A", "CD1C", "CLEC9A"],
    },
    "hsc": {
        "Quiescent HSC": ["RBP1", "LRAT", "CYGB"],
        "Activated HSC": ["COL1A1", "COL3A1", "PDGFRA"],
        "Myofibroblast": ["ACTA2", "TAGLN", "MYL9"],
        "Portal fibroblast": ["COL15A1", "PI16", "DPT"],
    },
}


def broad_label(subtype: str, lineage: str) -> str:
    label = str(subtype)
    low = label.lower()
    if "contaminant" in low or "low_quality" in low or "ambiguous" in low:
        return "Contaminant"
    if lineage == "myeloid":
        if "cycling" in low:
            return "Cycling"
        if any(token in low for token in ("cdc", "pdc", "dc", "lamp3")):
            return "Dendritic"
        if "monocyte" in low or "fcn1" in low or "fcgr3a" in low:
            return "Monocyte"
        if "kupffer" in low or "marco" in low:
            return "Kupffer"
        if any(token in low for token in ("macrophage", "samac", "spp1", "trem2", "gpnmb")):
            return "Macrophage"
        return "Other myeloid"
    if "quiescent" in low:
        return "Quiescent HSC"
    if "activated" in low:
        return "Activated HSC"
    if "myofibroblast" in low:
        return "Myofibroblast"
    if "portal" in low:
        return "Portal fibroblast"
    return "Other HSC"


def load_object(dataset: str, lineage: str):
    path = DATA / f"{dataset}_{lineage}_reclustered.h5ad"
    return ad.read_h5ad(path, backed="r")


def save_figure1_source_data():
    """Export the cell-level coordinates and labels used by the atlas."""
    records = []
    for dataset in DATASETS:
        for lineage in LINEAGES:
            adata = load_object(dataset, lineage)
            coordinates = np.asarray(adata.obsm["X_umap"])
            observations = adata.obs.copy()
            donor_column = "patient_id" if "patient_id" in observations else ("sample_id" if "sample_id" in observations else "gsm")
            for index, cell_id in enumerate(adata.obs_names.astype(str)):
                row = {
                    "dataset": dataset,
                    "lineage": lineage,
                    "cell_id": cell_id,
                    "UMAP1": float(coordinates[index, 0]),
                    "UMAP2": float(coordinates[index, 1]),
                    "subtype": str(observations.iloc[index]["subtype"]),
                    "broad_annotation": broad_label(str(observations.iloc[index]["subtype"]), lineage),
                    "donor_id": str(observations.iloc[index][donor_column]),
                }
                for column in ("disease", "fibrosis", "sample_id", "patient_id", "gsm"):
                    if column in observations:
                        row[column] = str(observations.iloc[index][column])
                records.append(row)
            adata.file.close()
    source_data = pd.DataFrame(records)
    source_data.to_csv(SOURCE_DATA / "Source_Data_Figure1_single_cell_atlas.tsv.gz", sep="\t", index=False, compression="gzip")


def style_axes(axis):
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(False)
    axis.set_facecolor("white")


def plot_umap(axis, adata, lineage: str):
    coordinates = np.asarray(adata.obsm["X_umap"])
    labels = pd.Series(adata.obs["subtype"].astype(str).map(lambda value: broad_label(value, lineage)).to_numpy())
    palette = MYELOID_COLORS if lineage == "myeloid" else HSC_COLORS
    order = [key for key in palette if key in set(labels)]
    for label in order:
        mask = labels.eq(label).to_numpy()
        axis.scatter(coordinates[mask, 0], coordinates[mask, 1], s=3.2, alpha=0.72,
                     c=palette[label], linewidths=0, rasterized=True)
    style_axes(axis)
    return labels


def save_main_figure():
    sns.set_theme(style="white", context="paper", font="DejaVu Sans")
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.6), dpi=300,
                             gridspec_kw={"wspace": 0.05, "hspace": 0.14})
    panel_labels = list("ABCDEF")
    legend_handles = {}
    for col, dataset in enumerate(DATASETS):
        for row, lineage in enumerate(LINEAGES):
            adata = load_object(dataset, lineage)
            labels = plot_umap(axes[row, col], adata, lineage)
            axes[row, col].text(0.0, 1.02, panel_labels[row * 3 + col], transform=axes[row, col].transAxes,
                                ha="left", va="bottom", fontsize=8, fontweight="bold", clip_on=False)
            if row == 0:
                axes[row, col].set_title(DISPLAY[dataset], fontsize=9, fontweight="bold", pad=3)
            palette = MYELOID_COLORS if lineage == "myeloid" else HSC_COLORS
            for label in labels.unique():
                legend_handles.setdefault((lineage, label),
                    Line2D([0], [0], marker="o", linestyle="", markersize=5,
                           markerfacecolor=palette.get(label, "#666666"), markeredgewidth=0, label=label))
            adata.file.close()

    fig.text(0.01, 0.72, "Myeloid compartment", rotation=90, rotation_mode="anchor", va="center", ha="left",
             fontsize=8.5, fontweight="bold")
    fig.text(0.01, 0.28, "HSC/fibroblast compartment", rotation=90, rotation_mode="anchor", va="center", ha="left",
             fontsize=8.5, fontweight="bold")
    myeloid_legend = [legend_handles[("myeloid", key)] for key in MYELOID_COLORS if ("myeloid", key) in legend_handles]
    hsc_legend = [legend_handles[("hsc", key)] for key in HSC_COLORS if ("hsc", key) in legend_handles]
    fig.subplots_adjust(left=0.055, right=0.995, top=0.88, bottom=0.20, wspace=0.05, hspace=0.14)
    fig.legend(handles=myeloid_legend, loc="lower center", bbox_to_anchor=(0.50, 0.090), ncol=4,
               frameon=False, fontsize=7, title="Myeloid annotation", title_fontsize=7.2, handletextpad=0.3,
               columnspacing=1.0)
    fig.legend(handles=hsc_legend, loc="lower center", bbox_to_anchor=(0.50, 0.012), ncol=5,
               frameon=False, fontsize=7, title="HSC/fibroblast annotation", title_fontsize=7.2, handletextpad=0.3,
               columnspacing=1.0)
    fig.suptitle("Cross-cohort single-cell atlas of myeloid and HSC/fibroblast compartments",
                 fontsize=10.5, fontweight="bold", y=0.985)
    fig.savefig(FIG / "Figure1_myeloid_HSC_atlas_v2.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / "Figure1_myeloid_HSC_atlas_v2.svg", bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / "Figure1_myeloid_HSC_atlas_v2.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / "Figure1_myeloid_HSC_atlas_v2.tiff", dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _long_metadata():
    rows = []
    for dataset in DATASETS:
        for lineage in LINEAGES:
            adata = load_object(dataset, lineage)
            obs = adata.obs.copy()
            obs["dataset"] = dataset
            obs["lineage"] = lineage
            obs["broad_label"] = obs["subtype"].astype(str).map(lambda value: broad_label(value, lineage))
            donor_column = "patient_id" if "patient_id" in obs else ("sample_id" if "sample_id" in obs else "gsm")
            obs["donor_id"] = obs[donor_column].astype(str)
            keep = [column for column in ["dataset", "lineage", "broad_label", "disease", "total_counts", "n_genes", "pct_mito", "scDblFinder_class", "donor_id"] if column in obs]
            rows.append(obs[keep].reset_index(drop=True))
            adata.file.close()
    return pd.concat(rows, ignore_index=True)


def save_qc_composition():
    meta = _long_metadata()
    sns.set_theme(style="whitegrid", context="paper", font="DejaVu Sans")
    fig, axes = plt.subplots(2, 3, figsize=(14.2, 8.2), dpi=300, constrained_layout=True)
    qc_columns = [("n_genes", "Detected genes"), ("total_counts", "UMI counts"), ("pct_mito", "Mitochondrial fraction")]
    for index, (column, title) in enumerate(qc_columns):
        axis = axes[0, index]
        sns.boxenplot(data=meta, x="dataset", y=column, hue="lineage", order=DATASETS,
                      palette={"myeloid": "#0072B2", "hsc": "#D55E00"}, showfliers=False, linewidth=0.5, ax=axis)
        axis.set_title(title, fontsize=10, fontweight="bold")
        axis.set_xlabel("")
        axis.set_ylabel("")
        axis.tick_params(axis="x", labelsize=8)
        plt.setp(axis.get_xticklabels(), rotation=25, ha="right", rotation_mode="anchor")
        if axis.get_legend() is not None:
            axis.get_legend().remove()
    comp = (meta.groupby(["dataset", "lineage", "disease", "broad_label"], observed=True)
            .size().rename("n").reset_index())
    comp["fraction"] = comp["n"] / comp.groupby(["dataset", "lineage", "disease"], observed=True)["n"].transform("sum")
    for index, lineage in enumerate(LINEAGES):
        axis = axes[1, index]
        subset = comp[comp["lineage"] == lineage]
        pivot = subset.pivot_table(index=["dataset", "disease"], columns="broad_label", values="fraction", fill_value=0)
        palette = MYELOID_COLORS if lineage == "myeloid" else HSC_COLORS
        pivot = pivot[[key for key in palette if key in pivot.columns]]
        pivot.plot(kind="bar", stacked=True, ax=axis, color=[palette[key] for key in pivot.columns], width=0.82,
                   legend=False, edgecolor="white", linewidth=0.25)
        axis.set_title(f"{lineage.capitalize()} composition by disease", fontsize=10, fontweight="bold")
        axis.set_xlabel("")
        axis.set_ylabel("Cell fraction")
        axis.set_xticklabels([f"{dataset}\n{disease}" for dataset, disease in pivot.index], rotation=0,
                             ha="center", fontsize=7)
        axis.set_ylim(0, 1)
    axis = axes[1, 2]
    summary = (meta.groupby(["dataset", "lineage"], observed=True)
               .agg(n_cells=("broad_label", "size"), n_units=("donor_id", "nunique"))
               .reset_index())
    summary["label"] = summary["dataset"] + "\n" + summary["lineage"].str.capitalize()
    bars = axis.bar(summary["label"], summary["n_cells"],
                    color=summary["lineage"].map({"myeloid": "#0072B2", "hsc": "#D55E00"}),
                    width=0.72, edgecolor="white", linewidth=0.4)
    axis.set_ylim(0, float(summary["n_cells"].max()) * 1.18)
    axis.set_title("Cells and donors retained for atlas", fontsize=10, fontweight="bold", pad=8)
    axis.set_xlabel("")
    axis.set_ylabel("Retained cells")
    axis.tick_params(axis="x", labelsize=7)
    axis.grid(axis="x", visible=False)
    axis.grid(axis="y", visible=False)
    fig.suptitle("Single-cell quality control and compartment composition", fontsize=13, fontweight="bold")
    stems = [FIG / "Supplementary_Figure_scRNA_QC_composition",
             REVISION_FIG / "Supplementary_Figure2_scRNA_QC_composition"]
    for stem in stems:
        fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
        fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
        fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")
        fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    meta.to_csv(SOURCE_DATA / "Source_Data_Supplementary_Figure2_scRNA_QC_composition.tsv.gz",
                sep="\t", index=False, compression="gzip")
    plt.close(fig)


def save_marker_validation():
    sns.set_theme(style="white", context="paper", font="DejaVu Sans")
    fig, axes = plt.subplots(1, 2, figsize=(14.8, 5.3), dpi=300, constrained_layout=True)
    source_tables = []
    for axis, lineage in zip(axes, LINEAGES):
        records = []
        groups = list(MARKERS[lineage])
        genes = [gene for group in groups for gene in MARKERS[lineage][group]]
        for dataset in DATASETS:
            adata = load_object(dataset, lineage)
            present = [gene for gene in genes if gene in adata.var_names]
            if present:
                matrix = adata[:, present].to_memory().X
                if hasattr(matrix, "toarray"):
                    matrix = matrix.toarray()
                matrix = np.asarray(matrix)
                expr = pd.DataFrame(matrix, columns=present)
                expr["broad_label"] = adata.obs["subtype"].astype(str).map(lambda value: broad_label(value, lineage)).to_numpy()
                for broad in groups:
                    subset = expr[expr["broad_label"] == broad]
                    if subset.empty:
                        continue
                    for gene in present:
                        records.append({"dataset": dataset, "group": broad, "gene": gene,
                                        "mean": float(subset[gene].mean()),
                                        "fraction": float((subset[gene] > 0).mean())})
            adata.file.close()
        table = pd.DataFrame(records)
        if table.empty:
            axis.set_axis_off()
            continue
        table["gene"] = pd.Categorical(table["gene"], categories=genes, ordered=True)
        table["group"] = pd.Categorical(table["group"], categories=groups, ordered=True)
        summary = table.groupby(["group", "gene"], observed=True).agg(mean=("mean", "mean"), fraction=("fraction", "mean")).reset_index()
        summary.insert(0, "lineage", lineage)
        source_tables.append(summary.copy())
        matrix = summary.pivot(index="group", columns="gene", values="mean").reindex(index=groups, columns=genes)
        size = summary.pivot(index="group", columns="gene", values="fraction").reindex(index=groups, columns=genes)
        x_positions, y_positions = np.meshgrid(np.arange(len(genes)), np.arange(len(groups)))
        plot = axis.scatter(x_positions.ravel(), y_positions.ravel(), s=(size.to_numpy().ravel() * 220) + 5,
                            c=matrix.to_numpy().ravel(), cmap="magma", vmin=0, vmax=max(1, float(np.nanpercentile(matrix, 98))),
                            edgecolors="white", linewidths=0.3)
        axis.set_xticks(np.arange(len(genes)), genes, fontsize=7.2)
        plt.setp(axis.get_xticklabels(), rotation=55, ha="right", rotation_mode="anchor")
        axis.set_yticks(np.arange(len(groups)), groups, fontsize=8)
        axis.set_title(f"{lineage.capitalize()} marker validation", fontsize=11, fontweight="bold")
        axis.set_xlabel("Marker gene")
        axis.set_ylabel("Broad annotation")
        fig.colorbar(plot, ax=axis, fraction=0.035, pad=0.02, label="Mean log-normalized expression")
    fig.suptitle("Marker-panel support for curated cell annotations", fontsize=13, fontweight="bold")
    stems = [FIG / "Supplementary_Figure_scRNA_marker_validation",
             REVISION_FIG / "Supplementary_Figure3_scRNA_marker_validation"]
    for stem in stems:
        fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
        fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
        fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")
        fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    if source_tables:
        pd.concat(source_tables, ignore_index=True).to_csv(
            SOURCE_DATA / "Source_Data_Supplementary_Figure3_scRNA_marker_validation.tsv.gz",
            sep="\t", index=False, compression="gzip")
    plt.close(fig)


def save_workflow():
    """Create a concise workflow schematic for the actual scRNA analysis chain."""
    fig, axis = plt.subplots(figsize=(14.0, 4.6), dpi=300)
    axis.set_xlim(0, 14)
    axis.set_ylim(0, 4.6)
    axis.axis("off")
    blocks = [
        (0.25, 2.65, 1.55, 1.05, "Public\nscRNA cohorts", "#EAF2F8"),
        (2.05, 2.65, 1.55, 1.05, "QC +\nscDblFinder", "#E8F6F3"),
        (3.85, 2.65, 1.55, 1.05, "Myeloid / HSC\nsubset", "#FEF5E7"),
        (5.65, 2.65, 1.55, 1.05, "Per-cohort\nreclustering", "#FDEDEC"),
        (7.45, 2.65, 1.75, 1.05, "Marker-guided\nannotation", "#F4ECF7"),
        (9.45, 2.65, 1.75, 1.05, "Ferro-aging +\nstate scores", "#E8F8F5"),
        (11.45, 2.65, 2.15, 1.05, "Patient-level\nassociation", "#FDF2E9"),
    ]
    for x, y, width, height, label, color in blocks:
        patch = FancyBboxPatch((x, y), width, height, boxstyle="round,pad=0.03,rounding_size=0.08",
                               linewidth=1.2, edgecolor="#34495E", facecolor=color)
        axis.add_patch(patch)
        axis.text(x + width / 2, y + height / 2, label, ha="center", va="center",
                  fontsize=10, fontweight="bold", color="#25313B")
    for left, right in zip(blocks[:-1], blocks[1:]):
        x1 = left[0] + left[2] + 0.03
        x2 = right[0] - 0.03
        axis.add_patch(FancyArrowPatch((x1, left[1] + left[3] / 2), (x2, right[1] + right[3] / 2),
                                       arrowstyle="-|>", mutation_scale=13, linewidth=1.2, color="#566573"))
    downstream = [
        (4.05, 0.68, 2.05, 0.84, "Monocle3 / Slingshot\ntrajectory", "#F5EEF8"),
        (6.65, 0.68, 2.05, 0.84, "LIANA / CellChat\ncommunication", "#EBDEF0"),
        (9.25, 0.68, 2.05, 0.84, "NicheNet + virtual\nnetwork perturbation", "#D6EAF8"),
        (11.85, 0.68, 1.8, 0.84, "Bulk\nvalidation", "#D5F5E3"),
    ]
    for x, y, width, height, label, color in downstream:
        patch = FancyBboxPatch((x, y), width, height, boxstyle="round,pad=0.03,rounding_size=0.08",
                               linewidth=1.0, edgecolor="#5D6D7E", facecolor=color)
        axis.add_patch(patch)
        axis.text(x + width / 2, y + height / 2, label, ha="center", va="center", fontsize=8.8, color="#25313B")
    connections = [
        ((6.43, 2.64), (5.08, 1.58), 0.10),
        ((10.30, 2.64), (7.68, 1.58), 0.12),
        ((12.45, 2.64), (10.28, 1.58), 0.10),
        ((12.45, 2.64), (12.75, 1.58), -0.25),
    ]
    for start, target, curvature in connections:
        axis.add_patch(FancyArrowPatch(start, target, connectionstyle=f"arc3,rad={curvature}",
                                       arrowstyle="-|>", mutation_scale=11, linewidth=1.0, color="#566573"))
    axis.text(0.25, 4.24, "Actual single-cell analysis workflow", fontsize=14, fontweight="bold", color="#17202A")
    axis.text(0.25, 3.96, "Three public liver cohorts were processed independently before patient-level integration.",
              fontsize=9.3, color="#566573")
    axis.text(0.25, 0.18, "Arrows indicate analysis order; downstream outputs are exploratory unless independently validated.",
              fontsize=8.3, color="#7B7D7D")
    legacy_stem = FIG / "Figure1_scRNA_workflow"
    final_stem = REVISION_FIG / "Supplementary_Figure1_scRNA_workflow"
    for stem in (legacy_stem, final_stem):
        fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
        fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
        fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")
        fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    workflow_rows = [
        {"order": 1, "branch": "main", "step": "Public scRNA cohorts", "description": "Three public human liver single-cell cohorts."},
        {"order": 2, "branch": "main", "step": "QC + scDblFinder", "description": "Library-level quality control and doublet identification."},
        {"order": 3, "branch": "main", "step": "Myeloid / HSC subset", "description": "Lineage-restricted subset extraction."},
        {"order": 4, "branch": "main", "step": "Per-cohort reclustering", "description": "Independent cohort-specific reclustering."},
        {"order": 5, "branch": "main", "step": "Marker-guided annotation", "description": "Canonical marker-supported subtype annotation."},
        {"order": 6, "branch": "main", "step": "Ferro-aging + state scores", "description": "State and pathway scoring, including iron-homeostasis and senescence-related modules."},
        {"order": 7, "branch": "main", "step": "Patient-level association", "description": "Patient/sample-level aggregation and matched analyses."},
        {"order": 8, "branch": "downstream", "step": "Monocle3 / Slingshot trajectory", "description": "Transcriptional state ordering."},
        {"order": 9, "branch": "downstream", "step": "LIANA / CellChat communication", "description": "Expression-based candidate communication inference."},
        {"order": 10, "branch": "downstream", "step": "NicheNet + virtual network perturbation", "description": "Exploratory ligand prioritization and counterfactual network prediction."},
        {"order": 11, "branch": "downstream", "step": "Bulk validation", "description": "Independent complete-bulk validation."},
    ]
    pd.DataFrame(workflow_rows).to_csv(SOURCE_DATA / "Source_Data_Supplementary_Figure1_workflow.tsv.gz", sep="\t", index=False, compression="gzip")
    plt.close(fig)


if __name__ == "__main__":
    save_main_figure()
    save_qc_composition()
    save_marker_validation()
    save_workflow()
    save_figure1_source_data()
    print("Figure 1 v2 and supplementary scRNA panels written to", FIG)
