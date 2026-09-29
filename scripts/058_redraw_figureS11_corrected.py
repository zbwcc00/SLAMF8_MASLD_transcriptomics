#!/usr/bin/env python3
import os
"""Redraw the virtual SLAMF8 perturbation audit with donor-level uncertainty.

The figure is deliberately descriptive: it does not display legacy Fisher
P values and does not label the counterfactual as an experimental knockout.
"""
from pathlib import Path

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
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
IN_DIR = ROOT / "results" / "revision_p0" / "virtual_perturbation"
OUT_DIR = ROOT / "figures" / "revision_p0"
SOURCE_DIR = ROOT / "source_data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
SOURCE_DIR.mkdir(parents=True, exist_ok=True)

DATASET_ORDER = ["GSE344087", "GSE298719", "GSE212837"]
DATASET_COLORS = {"GSE344087": "#0072B2", "GSE298719": "#D55E00", "GSE212837": "#009E73"}
TARGETS = [
    "FARG95", "Iron_homeostasis", "Senescence_SASP", "Inflammation",
    "ligand_SPP1", "ligand_PDGFB", "ligand_TGFB1", "ligand_CXCL16",
]
TARGET_LABELS = {
    "FARG95": "FARG95",
    "Iron_homeostasis": "Iron homeostasis",
    "Senescence_SASP": "Senescence/SASP",
    "Inflammation": "Inflammation",
    "ligand_SPP1": "Ligand: SPP1",
    "ligand_PDGFB": "Ligand: PDGFB",
    "ligand_TGFB1": "Ligand: TGFB1",
    "ligand_CXCL16": "Ligand: CXCL16",
}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "axes.labelsize": 8.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "savefig.dpi": 600,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def save(fig: plt.Figure) -> None:
    stem = OUT_DIR / "FigureS11_SLAMF8_counterfactual_audit_corrected"
    fig.savefig(f"{stem}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(f"{stem}.svg", bbox_inches="tight", facecolor="white")
    fig.savefig(f"{stem}.tiff", dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    audit = pd.read_csv(IN_DIR / "donor_level_counterfactual_audit.tsv", sep="\t")
    cross = pd.read_csv(IN_DIR / "cross_dataset_counterfactual_audit.tsv", sep="\t")
    audit = audit[audit.target.isin(TARGETS)].copy()
    audit["target"] = pd.Categorical(audit["target"], categories=TARGETS, ordered=True)
    audit["dataset"] = pd.Categorical(audit["dataset"], categories=DATASET_ORDER, ordered=True)
    audit = audit.sort_values(["target", "dataset"])
    cross = cross[cross.target.isin(TARGETS)].copy()
    cross["target"] = pd.Categorical(cross["target"], categories=TARGETS, ordered=True)
    cross = cross.sort_values("target")

    fig = plt.figure(figsize=(14.2, 12.0), facecolor="white")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.65, 1.0], width_ratios=[1.8, 1.0],
                          hspace=0.42, wspace=0.28)

    # A: cohort-specific donor-level deltas, with a separate x scale per target.
    ax = fig.add_subplot(gs[0, :])
    y_positions = np.arange(len(TARGETS))[::-1]
    offsets = {"GSE344087": 0.22, "GSE298719": 0.0, "GSE212837": -0.22}
    for yi, target in zip(y_positions, TARGETS):
        block = audit[audit.target == target]
        max_abs = np.nanmax(np.abs(block[["mean_bootstrap_ci_low", "mean_bootstrap_ci_high", "mean_delta"]].to_numpy(float)))
        pad = max(max_abs * 0.18, 1e-5)
        for _, row in block.iterrows():
            y = yi + offsets[str(row.dataset)]
            color = DATASET_COLORS[str(row.dataset)]
            ax.errorbar(row.mean_delta, y,
                        xerr=[[row.mean_delta - row.mean_bootstrap_ci_low],
                              [row.mean_bootstrap_ci_high - row.mean_delta]],
                        fmt="o", color=color, markerfacecolor=color, markeredgecolor="white",
                        markeredgewidth=0.55, markersize=5.5, elinewidth=1.15, capsize=2.2,
                        zorder=3)
        ax.axhline(yi - 0.5, color="#E5E9EC", linewidth=0.65, zorder=0)
    ax.axvline(0, color="#263238", linewidth=0.85, zorder=1)
    ax.set_yticks(y_positions, [TARGET_LABELS[t] for t in TARGETS])
    ax.set_xlabel("Predicted donor-level mean change after model-based SLAMF8 reduction (raw score units)")
    ax.set_title("A  Donor-level counterfactual changes with bootstrap 95% CIs", loc="left", pad=9)
    ax.grid(axis="x", alpha=0.18, linewidth=0.6)
    ax.text(0.0, -0.16, "Intervals are donor bootstrap intervals for the predicted mean delta; scales are raw and target-specific.",
            transform=ax.transAxes, fontsize=7.5, color="#6F7A80")

    # B: direction consistency and standardized effect size.
    ax_b = fig.add_subplot(gs[1, 0])
    direction = cross.set_index("target")["median_standardized_mean_delta"].reindex(TARGETS)
    consistency = cross.set_index("target").reindex(TARGETS)
    colors = ["#B2182B" if v < 0 else "#2166AC" for v in direction]
    y = np.arange(len(TARGETS))
    ax_b.axvline(0, color="#263238", linewidth=0.8)
    ax_b.barh(y, direction.values, color=colors, alpha=0.86, height=0.62)
    ax_b.set_yticks(y, [TARGET_LABELS[t] for t in TARGETS])
    ax_b.set_xlabel("Median standardized mean delta")
    ax_b.set_title("B  Cross-cohort direction summary", loc="left", pad=9)
    ax_b.grid(axis="x", alpha=0.18, linewidth=0.6)
    ax_b.text(0.0, -0.19, "Direction counts summarize cohort-level predicted deltas; they are not causal tests.",
              transform=ax_b.transAxes, fontsize=7.5, color="#6F7A80")

    # C: group-CV R2.
    ax_c = fig.add_subplot(gs[1, 1])
    cv = audit.pivot(index="target", columns="dataset", values="group_cv_r2_mean").reindex(TARGETS, columns=DATASET_ORDER)
    sns.heatmap(cv, ax=ax_c, annot=True, fmt=".2f", cmap="YlGnBu", vmin=-0.9, vmax=0.9,
                linewidths=0.8, linecolor="white", cbar_kws={"label": "Group-CV R²", "shrink": 0.82, "pad": 0.03},
                annot_kws={"fontsize": 7.4})
    ax_c.set_yticklabels([TARGET_LABELS[t] for t in TARGETS], rotation=0)
    ax_c.set_xlabel("")
    ax_c.set_ylabel("")
    ax_c.tick_params(axis="x", rotation=0, pad=8, length=0)
    for tick in ax_c.get_xticklabels():
        tick.set_rotation_mode("anchor")
    ax_c.tick_params(axis="y", length=0)
    ax_c.set_title("C  Donor-aware predictive fit", loc="left", pad=9)
    ax_c.text(0.0, -0.19, "R² is a predictive-fit diagnostic under group cross-validation.",
              transform=ax_c.transAxes, fontsize=7.5, color="#6F7A80")

    fig.suptitle("Figure S11  SLAMF8 counterfactual network audit", fontsize=15, fontweight="bold",
                 x=0.02, ha="left", y=0.995)
    fig.text(0.02, 0.008,
             "Counterfactual network prediction only; not experimental SLAMF8 knockout. Legacy Fisher P values are intentionally not shown.",
             fontsize=8.2, color="#6F7A80")
    save(fig)
    source = pd.concat([
        audit.assign(panel="A_donor_level_counterfactual"),
        cross.assign(panel="B_cross_cohort_direction_summary"),
    ], ignore_index=True, sort=False)
    source.to_csv(SOURCE_DIR / "Source_Data_Supplementary_Figure11_SLAMF8_counterfactual_audit.tsv.gz",
                  sep="\t", index=False, compression="gzip")
    print(f"Wrote {OUT_DIR / 'FigureS11_SLAMF8_counterfactual_audit_corrected.pdf'}")


if __name__ == "__main__":
    main()
