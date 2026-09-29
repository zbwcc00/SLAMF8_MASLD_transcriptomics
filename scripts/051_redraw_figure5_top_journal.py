#!/usr/bin/env python3
import os
"""Submission redraw and source-data export for Figure 5 bulk validation."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
RES = ROOT / "results" / "bulk_validation"
EXP = RES / "expanded_public"
FIG = ROOT / "figures" / "revision_p0"
SOURCE_DATA = ROOT / "source_data"
FIG.mkdir(parents=True, exist_ok=True)
SOURCE_DATA.mkdir(parents=True, exist_ok=True)

COLORS = {
    "GSE162694": "#0072B2", "GSE130970": "#D55E00", "gray": "#6F7A80",
    "green": "#009E73", "orange": "#E69F00", "purple": "#CC79A7", "ink": "#263238",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.titleweight": "bold",
    "axes.labelsize": 6.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.7,
    "xtick.color": COLORS["ink"],
    "ytick.color": COLORS["ink"],
    "text.color": COLORS["ink"],
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


def load():
    fibrosis = pd.read_csv(RES / "full_bulk_fibrosis_validation.tsv", sep="\t")
    modules = pd.read_csv(RES / "full_bulk_SLAMF8_module_correlations.tsv", sep="\t")
    cov = pd.read_csv(RES / "full_bulk_covariate_sensitivity.tsv", sep="\t")
    scores = {ds: pd.read_csv(RES / f"{ds}_full_bulk_scores.tsv", sep="\t") for ds in ["GSE162694", "GSE130970"]}
    exp = pd.read_csv(EXP / "expanded_public_validation.tsv", sep="\t")
    add = pd.read_csv(EXP / "additional_human_validation.tsv", sep="\t")
    return fibrosis, modules, cov, scores, exp, add


def panel_a_b(ax, data, ds, result, panel_label):
    d = data.loc[~data.normal_histology].copy()
    d["fibrosis"] = d.fibrosis.astype(int)
    order = sorted(d.fibrosis.unique())
    sns.boxplot(data=d, x="fibrosis", y="SLAMF8", order=order, ax=ax,
                color="#DCEAF3" if ds == "GSE162694" else "#F6DDD3", width=0.62,
                fliersize=0, linewidth=0.75)
    sns.stripplot(data=d, x="fibrosis", y="SLAMF8", order=order, ax=ax,
                  color=COLORS[ds], size=2.2, alpha=0.62, jitter=0.16, linewidth=0)
    ax.text(-0.11, 1.10, panel_label, transform=ax.transAxes, fontsize=8.5,
            fontweight="bold", va="top", ha="left", clip_on=False)
    ax.set_title(
        f"{ds}  |  NAFLD biopsies, n={int(result.n)}\n"
        f"ρ={result.rho:.3f} · BH-FDR={result.fdr_within_dataset_scope:.2g} · F≥2 AUC={result.auc_f_ge2:.3f}",
        loc="left", fontsize=6.8, pad=5, linespacing=1.05,
    )
    ax.set_xlabel("Fibrosis stage (F0–F4)", fontsize=6.2)
    ax.set_ylabel("SLAMF8 expression" if ds == "GSE162694" else "", fontsize=6.2)
    ax.tick_params(labelsize=5.8, length=2, pad=2)
    ax.grid(axis="y", alpha=0.15, linewidth=0.5)


def panel_c(ax, fibrosis):
    features = ["SLAMF8", "HSC_activation", "FARG95", "Iron_homeostasis", "Senescence_SASP", "SPP1", "GPNMB"]
    display = {"HSC_activation": "HSC activation", "Iron_homeostasis": "Iron homeostasis", "Senescence_SASP": "Senescence/SASP"}
    d = fibrosis[(fibrosis.scope == "NAFLD_only") & fibrosis.feature.isin(features)]
    mat = d.pivot(index="feature", columns="dataset", values="rho").reindex(index=features, columns=["GSE162694", "GSE130970"])
    fdr = d.pivot(index="feature", columns="dataset", values="fdr_within_dataset_scope").reindex(index=features, columns=["GSE162694", "GSE130970"])
    labels = mat.copy().astype(object)
    for feature in mat.index:
        for dataset in mat.columns:
            value = mat.loc[feature, dataset]
            q = fdr.loc[feature, dataset]
            labels.loc[feature, dataset] = "" if pd.isna(value) else f"{value:+.2f}" + ("*" if q < 0.05 else "")
    sns.heatmap(mat, ax=ax, cmap="RdBu_r", center=0, vmin=-0.65, vmax=0.65, annot=labels, fmt="",
                linewidths=0.7, linecolor="white", cbar_kws={"label": "Spearman ρ", "shrink": 0.72, "pad": 0.03},
                annot_kws={"fontsize": 5.8})
    ax.text(-0.11, 1.10, "C", transform=ax.transAxes, fontsize=8.5, fontweight="bold", va="top", clip_on=False)
    ax.set_title("Within-NAFLD fibrosis associations", loc="left", pad=5, fontsize=7.2)
    ax.set_xlabel(""); ax.set_ylabel("")
    ax.set_yticklabels([display.get(label, label) for label in features], fontsize=5.8)
    ax.set_xticklabels(["GSE162694", "GSE130970"], fontsize=5.8)
    ax.tick_params(length=0)
    ax.text(0, -0.18, "* BH-FDR < 0.05 within cohort; modules are not uniformly replicated.", transform=ax.transAxes, fontsize=5.8, color=COLORS["gray"])


def panel_d(subspec, scores, modules):
    inner = GridSpecFromSubplotSpec(1, 2, subspec, wspace=0.42)
    axes = [plt.subplot(inner[0]), plt.subplot(inner[1])]
    for ax, ds in zip(axes, ["GSE162694", "GSE130970"]):
        d = scores[ds].loc[~scores[ds].normal_histology].copy()
        ax.scatter(d.SLAMF8, d.HSC_activation, s=10, alpha=0.70, color=COLORS[ds], edgecolor="white", linewidth=0.2, rasterized=True)
        coef = np.polyfit(d.SLAMF8, d.HSC_activation, 1)
        xx = np.linspace(d.SLAMF8.min(), d.SLAMF8.max(), 50)
        ax.plot(xx, coef[0] * xx + coef[1], color=COLORS["ink"], linewidth=1.0)
        row = modules[(modules.dataset == ds) & (modules.feature == "HSC_activation")].iloc[0]
        ax.set_title("")
        ax.text(0.00, 1.10, ds, transform=ax.transAxes, va="bottom", ha="left", fontsize=6.8, fontweight="bold", clip_on=False)
        ax.text(0.00, 1.03, f"ρ={row.rho_with_SLAMF8:.3f} · BH-FDR={row.fdr_within_dataset:.2g}", transform=ax.transAxes,
                va="bottom", ha="left", fontsize=5.2, clip_on=False)
        ax.set_xlabel("SLAMF8 expression", fontsize=5.8)
        ax.set_ylabel("HSC activation" if ds == "GSE162694" else "", fontsize=5.8)
        ax.tick_params(labelsize=5.4, length=2, pad=2)
        ax.grid(alpha=0.15, linewidth=0.5)
    axes[0].text(-0.11, 1.12, "D", transform=axes[0].transAxes, fontsize=8.5, fontweight="bold", va="top", clip_on=False)
    axes[1].text(0.0, -0.22, "Bulk-tissue association; cell–cell causality is not inferred.", transform=axes[1].transAxes, fontsize=5.8, color=COLORS["gray"])


def make_main(fibrosis, modules, scores):
    fig = plt.figure(figsize=(7.2, 7.7), facecolor="white")
    grid = GridSpec(2, 2, figure=fig, height_ratios=[1.02, 1.02], hspace=0.58, wspace=0.47,
                    left=0.16, right=0.95, top=0.93, bottom=0.09)
    for ax, ds, label in zip([fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])], ["GSE162694", "GSE130970"], ["A", "B"]):
        row = fibrosis[(fibrosis.dataset == ds) & (fibrosis.scope == "NAFLD_only") & (fibrosis.feature == "SLAMF8")].iloc[0]
        panel_a_b(ax, scores[ds], ds, row, label)
    panel_c(fig.add_subplot(grid[1, 0]), fibrosis)
    panel_d(grid[1, 1], scores, modules)
    fig.suptitle("Independent complete-bulk validation of the SLAMF8–HSC axis", fontsize=9.2, fontweight="bold", x=0.16, ha="left", y=0.985)
    fig.text(0.16, 0.018, "Expression units are cohort-specific (CPM+1 versus TPM+1). F≥2 AUC is descriptive and untrained, not a clinical model.", fontsize=5.8, color=COLORS["gray"])
    save_figure(fig, "Figure5_independent_bulk_validation_v3_corrected")


def make_covariate(cov):
    features = ["SLAMF8", "HSC_activation", "FARG95", "Iron_homeostasis"]
    cov = cov[(cov.scope == "NAFLD_only") & cov.feature.isin(features)].copy()
    cov["label"] = cov.dataset + " | " + cov.feature
    order = [f"{ds} | {feat}" for ds in ["GSE162694", "GSE130970"] for feat in features]
    cov["label"] = pd.Categorical(cov.label, categories=order, ordered=True); cov = cov.sort_values("label")
    y = np.arange(len(order))
    fig, ax = plt.subplots(figsize=(7.2, 4.9), facecolor="white")
    for label, block in cov.groupby("label", observed=True):
        block = block.set_index("adjustment").reindex(["age_sex", "age_sex_NAS"])
        yi = order.index(str(label))
        ax.plot(block["beta_fibrosis_stage_per_sd"], [yi - 0.13, yi + 0.13], color=COLORS["gray"], linewidth=1.1)
        ax.scatter(block.loc["age_sex", "beta_fibrosis_stage_per_sd"], yi - 0.13, s=20, color=COLORS["GSE162694"], zorder=3)
        ax.scatter(block.loc["age_sex_NAS", "beta_fibrosis_stage_per_sd"], yi + 0.13, s=20, color=COLORS["orange"], zorder=3)
    ax.axvline(0, color=COLORS["ink"], linewidth=0.7)
    cov_display = {"SLAMF8": "SLAMF8", "HSC_activation": "HSC activation", "FARG95": "FARG95", "Iron_homeostasis": "Iron homeostasis"}
    ax.set_yticks(y, ["  |  ".join([parts[0], cov_display.get(parts[1], parts[1])]) for parts in [str(x).split(" | ") for x in order]], fontsize=5.7)
    ax.set_xlabel("β for fibrosis stage per 1 SD feature increase", fontsize=6.2)
    ax.set_title("Supplementary Figure 9  Covariate sensitivity of bulk associations", loc="left", fontweight="bold", pad=6, fontsize=8)
    ax.tick_params(axis="x", labelsize=5.8, length=2); ax.grid(axis="x", alpha=0.15, linewidth=0.5)
    ax.text(0, -0.16, "Blue: age + sex; orange: age + sex + NAS. NAS attenuates several associations and reduces complete-case sample size.", transform=ax.transAxes, fontsize=5.8, color=COLORS["gray"])
    save_figure(fig, "Supplementary_Figure9_bulk_covariate_sensitivity")


def make_extended(exp, add):
    rows = []
    for _, row in exp[(exp.feature.isin(["SLAMF8", "HSC_activation"])) & exp.analysis.isin(["NASH_vs_NAFL", "NASH_vs_nonNASH"])].iterrows():
        rows.append({"label": f"{row.accession} | {row.analysis} | {row.feature}", "auc": row.auc})
    for _, row in add[(add.feature.isin(["SLAMF8", "HSC_activation"])) & add.analysis.isin(["MASH_vs_MASL", "NASH_vs_NO_NASH"])].iterrows():
        rows.append({"label": f"{row.accession} | {row.analysis} | {row.feature}", "auc": row.auc})
    d = pd.DataFrame(rows).drop_duplicates("label").sort_values("auc")
    d["display_label"] = d["label"].str.replace("_", " ", regex=False).str.replace("nonNASH", "non-NASH", regex=False).str.replace("NO NASH", "no-NASH", regex=False).str.replace("HSC activation", "HSC activation", regex=False)
    fig, ax = plt.subplots(figsize=(7.2, max(4.8, 0.30 * len(d))), facecolor="white")
    y = np.arange(len(d)); colors = [COLORS["GSE162694"] if "SLAMF8" in x else COLORS["purple"] for x in d.label]
    ax.scatter(d.auc, y, s=24, c=colors, edgecolor="white", linewidth=0.4)
    ax.axvline(0.5, color=COLORS["ink"], linewidth=0.7, linestyle="--")
    ax.set_yticks(y, d.display_label, fontsize=5.7); ax.set_xlim(0.25, 1.0)
    ax.set_xlabel("AUC (point estimate; extended public cohorts)", fontsize=6.2)
    ax.set_title("Supplementary Figure 10  Extended public-cohort directional validation", loc="left", fontweight="bold", pad=6, fontsize=8)
    ax.tick_params(axis="x", labelsize=5.8, length=2); ax.grid(axis="x", alpha=0.15, linewidth=0.5)
    ax.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=COLORS["GSE162694"], markersize=4, label="SLAMF8"), Line2D([0], [0], marker="o", color="none", markerfacecolor=COLORS["purple"], markersize=4, label="HSC activation")], frameon=False, fontsize=5.6, loc="lower right")
    ax.text(0, -0.12, "Point estimates are exploratory and heterogeneous; this panel does not establish a clinical classifier.", transform=ax.transAxes, fontsize=5.8, color=COLORS["gray"])
    save_figure(fig, "Supplementary_Figure10_extended_public_validation")


def build_source_data(fibrosis, modules, cov, scores, exp, add):
    blocks = []
    for panel, table, source_file in [
        ("A_GSE162694_fibrosis_stage", scores["GSE162694"], "GSE162694_full_bulk_scores.tsv"),
        ("B_GSE130970_fibrosis_stage", scores["GSE130970"], "GSE130970_full_bulk_scores.tsv"),
        ("A_B_fibrosis_statistics", fibrosis[fibrosis.scope == "NAFLD_only"], "full_bulk_fibrosis_validation.tsv"),
        ("C_fibrosis_associations", fibrosis[(fibrosis.scope == "NAFLD_only") & fibrosis.feature.isin(["SLAMF8", "HSC_activation", "FARG95", "Iron_homeostasis", "Senescence_SASP", "SPP1", "GPNMB"])], "full_bulk_fibrosis_validation.tsv"),
        ("D_SLAMF8_HSC_associations", modules[modules.feature == "HSC_activation"], "full_bulk_SLAMF8_module_correlations.tsv"),
        ("S9_covariate_sensitivity", cov[(cov.scope == "NAFLD_only") & cov.feature.isin(["SLAMF8", "HSC_activation", "FARG95", "Iron_homeostasis"])], "full_bulk_covariate_sensitivity.tsv"),
        ("S10_expanded_public", exp[exp.feature.isin(["SLAMF8", "HSC_activation"])], "expanded_public_validation.tsv"),
        ("S10_additional_human", add[add.feature.isin(["SLAMF8", "HSC_activation"])], "additional_human_validation.tsv"),
    ]:
        frame = table.copy(); frame.insert(0, "panel", panel); frame.insert(1, "source_file", source_file); blocks.append(frame)
    pd.concat(blocks, ignore_index=True, sort=False).to_csv(SOURCE_DATA / "Source_Data_Figure5_bulk_validation.tsv.gz", sep="\t", index=False, compression="gzip")


def main():
    fibrosis, modules, cov, scores, exp, add = load()
    make_main(fibrosis, modules, scores)
    make_covariate(cov)
    make_extended(exp, add)
    build_source_data(fibrosis, modules, cov, scores, exp, add)


if __name__ == "__main__":
    main()
