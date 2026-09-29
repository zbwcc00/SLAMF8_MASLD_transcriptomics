#!/usr/bin/env python3
import os
"""Submission redraw and source-data export for Figure 4 communication audit."""

from pathlib import Path
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
COMM = ROOT / "results" / "communication"
AUDIT = COMM / "audit"
SPP1 = COMM / "cross_cohort_spp1"
FIG = ROOT / "figures" / "revision_p0"
SOURCE_DATA = ROOT / "source_data"
FIG.mkdir(parents=True, exist_ok=True)
SOURCE_DATA.mkdir(parents=True, exist_ok=True)

OKABE = {
    "blue": "#0072B2", "sky": "#56B4E9", "green": "#009E73",
    "orange": "#E69F00", "vermillion": "#D55E00", "gray": "#7A7A7A",
    "ink": "#263238",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.titleweight": "bold",
    "axes.labelsize": 7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.7,
    "xtick.color": OKABE["ink"],
    "ytick.color": OKABE["ink"],
    "text.color": OKABE["ink"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIG / f"{stem}.svg", bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"{stem}.tiff", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"{stem}.png", dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fisher_ci(rho: float, n: int) -> tuple[float, float]:
    if pd.isna(rho) or pd.isna(n) or n <= 4 or abs(rho) >= 1:
        return np.nan, np.nan
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    se = 1 / math.sqrt(n - 3)
    return tuple(np.tanh(z + np.array([-1.96, 1.96]) * se))


def load_tables():
    liana = pd.read_csv(AUDIT / "liana_candidate_audit.tsv", sep="\t")
    cellchat = pd.read_csv(AUDIT / "cellchat_candidate_audit.tsv", sep="\t")
    liana_path = pd.read_csv(COMM / "021_GSE298719_LIANA_Healthy_vs_MASLD_pathway_summary.tsv", sep="\t")
    cellchat_path = pd.read_csv(COMM / "022_GSE298719_CellChat_Healthy_vs_MASLD_pathway_summary.tsv", sep="\t")
    corr = pd.read_csv(SPP1 / "cross_cohort_SPP1_correlations.tsv", sep="\t")
    meta = pd.read_csv(SPP1 / "cross_cohort_SPP1_partial_condition_meta.tsv", sep="\t")
    expr_tests = pd.read_csv(AUDIT / "patient_level_ligand_receptor_expression_tests.tsv", sep="\t")
    return liana, cellchat, liana_path, cellchat_path, corr, meta, expr_tests


def panel_a(ax, liana, cellchat):
    ligands = ["SPP1", "PDGFC", "PDGFB", "TGFB1", "CXCL16"]
    columns = ["Healthy\nLIANA", "MASLD\nLIANA", "Healthy\nCellChat", "MASLD\nCellChat"]
    values = pd.DataFrame(index=ligands, columns=columns, dtype=float)
    labels = pd.DataFrame(index=ligands, columns=columns, dtype=object)
    for ligand in ligands:
        for condition in ["Healthy", "MASLD"]:
            row = liana[(liana.ligand == ligand) & (liana.condition == condition)]
            if len(row):
                row = row.iloc[0]
                den = float(row.n_inferred_rows) if pd.notna(row.n_inferred_rows) else 0.0
                num = float(row.n_active) if pd.notna(row.n_active) else 0.0
                values.loc[ligand, f"{condition}\nLIANA"] = num / den if den else np.nan
                labels.loc[ligand, f"{condition}\nLIANA"] = f"{int(num)}/{int(den)}" if den else "—"
            row = cellchat[(cellchat.ligand == ligand) & (cellchat.condition == condition)]
            if len(row):
                row = row.iloc[0]
                den = float(row.n_pairs) if pd.notna(row.n_pairs) else 0.0
                num = float(row.n_significant) if pd.notna(row.n_significant) else 0.0
                values.loc[ligand, f"{condition}\nCellChat"] = num / den if den else np.nan
                labels.loc[ligand, f"{condition}\nCellChat"] = f"{int(num)}/{int(den)}" if den else "—"
    sns.heatmap(values, ax=ax, cmap="YlGnBu", vmin=0, vmax=1, linewidths=0.8,
                linecolor="white", annot=labels, fmt="", annot_kws={"fontsize": 6.2},
                cbar_kws={"label": "Supported fraction", "shrink": 0.76, "pad": 0.03})
    for row_index in range(values.shape[0]):
        for column_index in range(values.shape[1]):
            if pd.isna(values.iloc[row_index, column_index]):
                ax.text(column_index + 0.5, row_index + 0.5, "—", ha="center", va="center",
                        fontsize=6.2, color=OKABE["ink"])
    ax.set_title("A  Candidate communication support", loc="left", pad=7, fontsize=7.5)
    ax.set_xlabel("")
    ax.set_ylabel("Candidate ligand", fontsize=6.5)
    ax.tick_params(axis="x", length=0, labelsize=6.2)
    ax.tick_params(axis="y", length=0, labelsize=6.5)
    ax.text(0, -0.20, "Cells show supported / inferred pairs; — indicates no inferred pairs.",
            transform=ax.transAxes, fontsize=6.1, color=OKABE["gray"])


def panel_b(subspec, liana_path, cellchat_path):
    inner = GridSpecFromSubplotSpec(1, 2, subspec, wspace=0.42, width_ratios=[1, 1])
    ax1, ax2 = plt.subplot(inner[0]), plt.subplot(inner[1])
    pathways = ["SPP1", "NOTCH", "TGFb", "PDGF", "TNF", "IL1"]
    lp = liana_path.set_index("focus_pathway").reindex(pathways)
    cp = cellchat_path.set_index("focus_pathway").reindex(["SPP1", "NOTCH", "PDGF", "Other"])
    y, h = np.arange(len(pathways)), 0.32
    ax1.barh(y - h / 2, lp["Healthy_active_fraction"], h, color=OKABE["sky"], label="Healthy")
    ax1.barh(y + h / 2, lp["MASLD_active_fraction"], h, color=OKABE["vermillion"], label="MASLD")
    ax1.set_yticks(y, pathways, fontsize=6.2); ax1.invert_yaxis(); ax1.set_xlim(0, 0.40)
    ax1.set_xlabel("LIANA active fraction", fontsize=6.3)
    ax1.set_title("LIANA pathway audit", loc="left", fontsize=7.0, pad=5)
    ax1.grid(axis="x", alpha=0.18, linewidth=0.5)
    ax1.legend(frameon=False, fontsize=5.6, loc="lower right", handlelength=1.3)
    y2 = np.arange(len(cp.index))
    ax2.barh(y2 - h / 2, cp["Healthy_interactions"].fillna(0), h, color=OKABE["sky"])
    ax2.barh(y2 + h / 2, cp["MASLD_interactions"].fillna(0), h, color=OKABE["vermillion"])
    ax2.set_yticks(y2, cp.index, fontsize=6.2); ax2.invert_yaxis(); ax2.set_xlim(0, 40)
    ax2.set_xlabel("CellChat significant interactions", fontsize=6.3)
    ax2.set_title("CellChat pathway audit", loc="left", fontsize=7.0, pad=5)
    ax2.grid(axis="x", alpha=0.18, linewidth=0.5)
    for axis in (ax1, ax2):
        axis.tick_params(axis="y", length=0, pad=2)
    ax2.text(-0.06, -0.25, "Healthy and MASLD are compared within each algorithm; scales are not shared.",
             transform=ax2.transAxes, fontsize=5.9, color=OKABE["gray"])
    ax1.text(-0.18, 1.12, "B", transform=ax1.transAxes, fontsize=8.5, fontweight="bold", va="top")


def panel_c(ax, corr, meta):
    features = ["SLAMF8_vs_SPP1_in_SAMac", "SPP1_source_vs_HSC_activation", "SPP1_source_vs_HSC_ITGAV"]
    labels = ["SLAMF8 ↔ SPP1 in SAMac", "SPP1-source → HSC activation", "SPP1-source → HSC ITGAV"]
    corr = corr[(corr.scope == "partial_condition") & corr.feature.isin(features)].copy()
    meta = meta[meta.feature.isin(features)].set_index("feature")
    rows = []
    for fi, feature in enumerate(features):
        for _, row in corr[corr.feature == feature].sort_values("dataset").iterrows():
            rows.append((fi, row.dataset, row.rho, row.n, False, row.fdr_within_scope))
        if feature in meta.index:
            row = meta.loc[feature]
            rows.append((fi, "Pooled", row.pooled_rho, row.n_total, True, row.fdr))
    positions, cursor = [], 0
    for fi in range(len(features)):
        block = [r for r in rows if r[0] == fi]
        for row in block:
            positions.append((row, cursor)); cursor += 1
        cursor += 0.7
    colors = {"GSE212837": OKABE["green"], "GSE298719": OKABE["blue"], "GSE344087": OKABE["orange"], "Pooled": OKABE["vermillion"]}
    for row, y in positions:
        _, dataset, rho, n, pooled, _ = row
        lo, hi = fisher_ci(rho, n); color = colors.get(dataset, OKABE["gray"])
        if pooled:
            ax.plot([lo, hi], [y, y], color=color, linewidth=1.8, zorder=2)
            ax.scatter(rho, y, marker="D", s=27, color=color, edgecolor="white", linewidth=0.4, zorder=3)
        else:
            ax.plot([lo, hi], [y, y], color=color, linewidth=1.0, zorder=1)
            ax.scatter(rho, y, s=18, color=color, edgecolor="white", linewidth=0.35, zorder=3)
    ax.axvline(0, color=OKABE["ink"], linewidth=0.7); ax.set_xlim(-1, 1)
    ax.set_xlabel("Partial Spearman ρ (95% CI)", fontsize=6.3)
    ax.set_yticks([y for _, y in positions], [row[1] for row, _ in positions], fontsize=5.8)
    ax.invert_yaxis()
    for fi, feature_label in enumerate(labels):
        ys = [y for row, y in positions if row[0] == fi]
        if ys:
            ax.text(-0.30, np.mean(ys), feature_label, transform=ax.get_yaxis_transform(),
                    ha="right", va="center", fontsize=5.9, fontweight="bold", clip_on=False)
            if fi < len(labels) - 1:
                ax.axhline(max(ys) + 0.85, color="#D9DDE0", linewidth=0.6)
    ax.set_title("C  Patient-level cross-cohort association audit", loc="left", pad=6, fontsize=7.3)
    ax.grid(axis="x", alpha=0.15, linewidth=0.5)
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=colors[k], markeredgecolor="white", markersize=3.8, label=k) for k in ["GSE344087", "GSE298719", "GSE212837"]]
    handles.append(Line2D([0], [0], marker="D", color="none", markerfacecolor=OKABE["vermillion"], markeredgecolor="white", markersize=3.8, label="Pooled"))
    ax.legend(handles=handles, frameon=False, fontsize=5.2, ncol=2, loc="lower left",
              bbox_to_anchor=(0.30, 1.05), handletextpad=0.3, columnspacing=0.8)


def panel_d(ax, meta):
    order = ["SPP1_source_vs_HSC_ITGB5", "SPP1_source_vs_HSC_ITGAV", "SPP1_source_vs_HSC_activation", "SPP1_source_vs_HSC_PDGFRA", "SPP1_source_vs_HSC_ENG", "SPP1_source_vs_HSC_ITGB1", "SPP1_source_vs_HSC_PDGFRB", "SPP1_source_vs_HSC_CD44", "SPP1_source_vs_HSC_ITGA9", "SPP1_source_vs_HSC_TGFBR1", "SPP1_source_vs_HSC_TGFBR2"]
    labels = [x.replace("SPP1_source_vs_HSC_", "") for x in order]
    d = meta.set_index("feature").reindex(order)
    y = np.arange(len(order)); sig = d["fdr"] < 0.05
    colors = [OKABE["vermillion"] if x else OKABE["gray"] for x in sig.fillna(False)]
    ax.axvline(0, color=OKABE["ink"], linewidth=0.7)
    ax.hlines(y, 0, d["pooled_rho"], color=colors, linewidth=2.3, alpha=0.58)
    ax.scatter(d["pooled_rho"], y, s=24, color=colors, edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_yticks(y, labels, fontsize=5.8); ax.invert_yaxis(); ax.set_xlim(-0.35, 0.91)
    ax.set_xlabel("Pooled partial Spearman ρ", fontsize=6.3)
    ax.set_title("D  SPP1-source receptor/program audit", loc="left", pad=6, fontsize=7.3)
    ax.grid(axis="x", alpha=0.15, linewidth=0.5)
    for yi, (_, row) in enumerate(d.iterrows()):
        if pd.notna(row.pooled_rho):
            ax.text(1.02, yi, f"ρ={row.pooled_rho:+.2f}\nFDR={row.fdr:.3f}",
                    transform=ax.get_yaxis_transform(), va="center", fontsize=5.2, clip_on=False)
    ax.text(0, -0.14, "Red: pooled FDR < 0.05; gray: weak or non-significant evidence.", transform=ax.transAxes, fontsize=5.9, color=OKABE["gray"])


def make_main(liana, cellchat, liana_path, cellchat_path, corr, meta):
    fig = plt.figure(figsize=(7.2, 7.9), facecolor="white")
    grid = GridSpec(2, 2, figure=fig, height_ratios=[0.94, 1.22], hspace=0.58, wspace=0.42,
                    left=0.18, right=0.98, top=0.93, bottom=0.09)
    panel_a(fig.add_subplot(grid[0, 0]), liana, cellchat)
    panel_b(grid[0, 1], liana_path, cellchat_path)
    panel_c(fig.add_subplot(grid[1, 0]), corr, meta)
    panel_d(fig.add_subplot(grid[1, 1]), meta)
    fig.suptitle("Macrophage–HSC communication audit in MASLD", fontsize=9.5, fontweight="bold", x=0.18, ha="left", y=0.985)
    fig.text(0.18, 0.018, "Expression-based inference and patient-level associations support a context-dependent SPP1 candidate; no causal signal is established.", fontsize=5.9, color=OKABE["gray"])
    save_figure(fig, "Figure4_macrophage_HSC_communication_v3_corrected")


def make_supplement(expr_tests):
    selected = ["mac_SPP1", "hsc_SPP1", "mac_SLAMF8", "hsc_ITGAV", "hsc_ITGB5", "hsc_HSC_activation"]
    d = expr_tests[expr_tests.feature.isin(selected)].copy()
    if d.empty:
        d = expr_tests.head(12).copy()
    d["display"] = d.feature.str.replace("mac_", "macrophage: ", regex=False).str.replace("hsc_", "HSC: ", regex=False)
    d = d.iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.2, 3.6), facecolor="white")
    y = np.arange(len(d)); delta = d.masld_minus_healthy.to_numpy(float)
    colors = [OKABE["vermillion"] if x >= 0 else OKABE["blue"] for x in delta]
    ax.axvline(0, color=OKABE["ink"], linewidth=0.7)
    ax.hlines(y, 0, delta, color=colors, linewidth=2.5, alpha=0.60)
    ax.scatter(delta, y, s=26, color=colors, edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_yticks(y, d.display, fontsize=6.1); ax.set_xlabel("MASLD − Healthy median expression", fontsize=6.3)
    ax.set_title("Supplementary Figure 8  Patient-level ligand/receptor expression audit", loc="left", pad=6, fontweight="bold", fontsize=8)
    ax.grid(axis="x", alpha=0.15, linewidth=0.5)
    for yi, (_, row) in enumerate(d.iterrows()):
        ax.text(max(delta.max(), 0.1) * 1.03, yi, f"P={row.p_value:.3f}", va="center", fontsize=5.2)
    ax.text(0, -0.18, "Permutation-based patient-level comparisons; medians are descriptive.", transform=ax.transAxes, fontsize=5.9, color=OKABE["gray"])
    save_figure(fig, "Supplementary_Figure8_patient_level_communication_audit")


def build_source_data(liana, cellchat, liana_path, cellchat_path, corr, meta, expr_tests):
    blocks = []
    for panel, table, source_file in [
        ("A_LIANA_candidate_support", liana, "liana_candidate_audit.tsv"),
        ("A_CellChat_candidate_support", cellchat, "cellchat_candidate_audit.tsv"),
        ("B_LIANA_pathways", liana_path, "021_GSE298719_LIANA_Healthy_vs_MASLD_pathway_summary.tsv"),
        ("B_CellChat_pathways", cellchat_path, "022_GSE298719_CellChat_Healthy_vs_MASLD_pathway_summary.tsv"),
        ("C_cross_cohort_correlations", corr, "cross_cohort_SPP1_correlations.tsv"),
        ("C_pooled_meta", meta, "cross_cohort_SPP1_partial_condition_meta.tsv"),
        ("D_pooled_receptor_program", meta, "cross_cohort_SPP1_partial_condition_meta.tsv"),
        ("S8_patient_expression_audit", expr_tests, "patient_level_ligand_receptor_expression_tests.tsv"),
    ]:
        frame = table.copy()
        frame.insert(0, "panel", panel)
        frame.insert(1, "source_file", source_file)
        blocks.append(frame)
    pd.concat(blocks, ignore_index=True, sort=False).to_csv(
        SOURCE_DATA / "Source_Data_Figure4_communication.tsv.gz", sep="\t", index=False, compression="gzip"
    )


def main():
    tables = load_tables()
    make_main(*tables[:5], tables[5])
    make_supplement(tables[6])
    build_source_data(*tables)


if __name__ == "__main__":
    main()
