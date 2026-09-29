#!/usr/bin/env python
import os
"""Plot the final condition-adjusted cross-cohort SPP1 review."""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
TABLE = ROOT / "results" / "communication" / "cross_cohort_spp1"
FIG = ROOT / "figures"
DATASETS = ["GSE344087", "GSE298719", "GSE212837"]


def main():
    tests = pd.read_csv(TABLE / "cross_cohort_SPP1_correlations.tsv", sep="\t")
    meta = pd.read_csv(TABLE / "cross_cohort_SPP1_partial_condition_meta.tsv", sep="\t")
    tests = tests[tests.scope.eq("partial_condition")]
    selected = ["SPP1_source_vs_HSC_activation", "SPP1_source_vs_HSC_ITGAV", "SPP1_source_vs_HSC_ITGB5", "SPP1_source_vs_HSC_PDGFRA", "SLAMF8_vs_SPP1_in_SAMac"]
    display = tests[tests.feature.isin(selected)].pivot(index="feature", columns="dataset", values="rho").reindex(index=selected, columns=DATASETS)
    fdr = tests[tests.feature.isin(selected)].pivot(index="feature", columns="dataset", values="fdr_within_scope").reindex(index=selected, columns=DATASETS)
    labels = display.copy().astype(object)
    for row in display.index:
        for col in display.columns:
            value = display.loc[row, col]
            labels.loc[row, col] = "" if pd.isna(value) else f"{value:+.2f}" + ("*" if fdr.loc[row, col] < 0.05 else "")

    sns.set_theme(style="ticks", context="paper", font_scale=1.1)
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), gridspec_kw={"width_ratios": [1.55, 1, 1]})
    sns.heatmap(display, ax=axes[0], annot=labels, fmt="", center=0, vmin=-1, vmax=1, cmap="RdBu_r", linewidths=1, linecolor="white", cbar_kws={"label": "Condition-adjusted Spearman ρ"})
    axes[0].set_title("Condition-adjusted patient-level associations\n* within-dataset FDR < 0.05", loc="left")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("")
    axes[0].tick_params(axis="x", rotation=20)
    axes[0].tick_params(axis="y", rotation=0)

    for ax, feature, title in zip(axes[1:], ["SPP1_source_vs_HSC_activation", "SLAMF8_vs_SPP1_in_SAMac"], ["SPP1-source state → HSC activation", "SLAMF8 ↔ SPP1 within SAMac"]):
        values = tests[tests.feature.eq(feature)].set_index("dataset").reindex(DATASETS)
        pooled = meta[meta.feature.eq(feature)]
        ax.axvline(0, color="#555555", linewidth=0.8)
        ax.errorbar(values.rho, np.arange(len(DATASETS)), xerr=None, fmt="o", color="#b2182b", markersize=6)
        if len(pooled):
            ax.scatter(pooled.iloc[0].pooled_rho, len(DATASETS), marker="D", s=48, color="#2166ac", zorder=3)
            ax.text(pooled.iloc[0].pooled_rho, len(DATASETS) + 0.25, f"pooled={pooled.iloc[0].pooled_rho:+.2f}\nFDR={pooled.iloc[0].fdr:.3f}", ha="center", va="bottom", fontsize=9)
        ax.set_yticks(np.arange(len(DATASETS) + 1))
        ax.set_yticklabels(DATASETS + ["pooled"])
        ax.set_ylim(len(DATASETS) + 0.6, -0.6)
        ax.set_xlim(-1, 1)
        ax.set_xlabel("Partial Spearman ρ")
        ax.set_title(title, loc="left")
        sns.despine(ax=ax)
    fig.suptitle("Final cross-cohort review of the SLAMF8–SPP1–HSC hypothesis", x=0.01, ha="left", y=1.02, fontsize=16)
    fig.tight_layout(w_pad=2)
    for suffix in ["png", "pdf"]:
        fig.savefig(FIG / f"016_cross_cohort_SPP1_review.{suffix}", dpi=320, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
