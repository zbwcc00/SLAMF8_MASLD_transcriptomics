#!/usr/bin/env python3
import os
"""Plot a bounded, non-causal audit of the in-silico SLAMF8 perturbation."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
RES = ROOT / "results" / "perturbation"
FIG = ROOT / "figures"
COLORS = {"GSE344087": "#0072B2", "GSE298719": "#D55E00", "GSE212837": "#009E73"}
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10,
    "axes.titleweight": "bold", "axes.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "savefig.dpi": 600,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def save(fig, stem):
    fig.savefig(FIG / f"{stem}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    summary = pd.read_csv(RES / "033_cross_dataset_SLAMF8_network_perturbation_summary.tsv", sep="\t")
    focus = ["FARG95", "Iron_homeostasis", "Senescence_SASP", "Inflammation", "ligand_SPP1", "ligand_PDGFB", "ligand_TGFB1", "ligand_CXCL16"]
    summary = summary[summary.target.isin(focus)].copy()
    summary["target"] = pd.Categorical(summary.target, categories=focus, ordered=True)
    summary = summary.sort_values(["target", "dataset"])
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.7), gridspec_kw={"width_ratios": [1.25, 1]}, facecolor="white")
    ax = axes[0]
    targets = focus
    for yi, target in enumerate(targets):
        block = summary[summary.target == target]
        for _, row in block.iterrows():
            ax.scatter(row.mean_delta, yi, s=38, color=COLORS.get(row.dataset, "#6F7A80"), edgecolor="white", linewidth=0.45, zorder=3)
        ax.plot([block.mean_delta.min(), block.mean_delta.max()], [yi, yi], color="#C7CDD1", linewidth=2, zorder=1)
    ax.axvline(0, color="#263238", linewidth=0.8)
    ax.set_yticks(np.arange(len(targets)), targets)
    ax.invert_yaxis()
    ax.set_xlabel("Predicted mean change after SLAMF8 reduction")
    ax.set_title("A  Cross-cohort perturbation direction", loc="left", pad=10)
    ax.grid(axis="x", alpha=0.15, linewidth=0.6)
    ax.legend(handles=[plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=COLORS[d], markersize=6, label=d) for d in ["GSE344087", "GSE298719", "GSE212837"]], frameon=False, fontsize=7.5, loc="lower right")
    ax.text(0, -0.17, "Dots are donor-level mean predicted deltas; near-zero effects and direction discordance are retained.", transform=ax.transAxes, fontsize=7.5, color="#6F7A80")

    ax = axes[1]
    cv = summary.pivot(index="target", columns="dataset", values="group_cv_r2_mean").reindex(targets)
    sns.heatmap(cv, ax=ax, annot=True, fmt=".2f", cmap="YlGnBu", vmin=-0.35, vmax=0.85, linewidths=1, linecolor="white", cbar_kws={"label": "Group-CV R²", "shrink": 0.78, "pad": 0.03}, annot_kws={"fontsize": 8})
    ax.set_title("B  Donor-aware prediction audit", loc="left", pad=10)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=20, length=0)
    ax.tick_params(axis="y", rotation=0, length=0)
    ax.text(0, -0.17, "R² summarizes predictive fit under group cross-validation; it is not evidence of biological causality.", transform=ax.transAxes, fontsize=7.5, color="#6F7A80")
    fig.suptitle("Supplementary Figure 7  In-silico SLAMF8 network perturbation audit", fontsize=14, fontweight="bold", x=0.02, ha="left", y=0.995)
    fig.text(0.02, 0.006, "Counterfactual network prediction only; not experimental knockout. Fisher P values from the legacy aggregation are intentionally not shown.", fontsize=8.2, color="#6F7A80")
    save(fig, "Supplementary_Figure7_SLAMF8_perturbation_audit")


if __name__ == "__main__":
    main()
