#!/usr/bin/env python
import os
"""Plot the frozen expanded public-cohort validation and robustness results."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
FREEZE = ROOT / "results" / "final_freeze"
RESULTS = ROOT / "results" / "bulk_validation" / "expanded_public"
FIGURES = ROOT / "figures"
FEATURES = ["SLAMF8", "SPP1", "FARG95", "HSC_activation"]
COLORS = {"SLAMF8": "#2b6ca3", "SPP1": "#c44e52", "FARG95": "#55a868", "HSC_activation": "#8172b3"}


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="ticks", context="paper", font_scale=1.05)
    fig, axes = plt.subplots(1, 3, figsize=(16.2, 5.1), gridspec_kw={"width_ratios": [1.35, 1.0, 1.0]})

    auc = pd.read_csv(FREEZE / "external_validation_bootstrap_auc.tsv", sep="\t")
    keep_analyses = ["NASH_vs_NAFL", "NASH_vs_NO_NASH", "MASH_vs_MASL"]
    auc = auc[auc.analysis.isin(keep_analyses) & auc.feature.isin(FEATURES)].copy()
    labels = {"NASH_vs_NAFL": "GSE167523\nNASH vs NAFL", "NASH_vs_NO_NASH": "GSE83452\nNASH vs no NASH", "MASH_vs_MASL": "GSE281797\nMASH vs MASL"}
    auc["label"] = auc.analysis.map(labels)
    rows = []
    for label in [labels[x] for x in keep_analyses]:
        for feature in FEATURES:
            rows.append((label, feature))
    y = np.arange(len(rows))[::-1]
    for yi, (label, feature) in zip(y, rows):
        row = auc[(auc.label == label) & (auc.feature == feature)].iloc[0]
        axes[0].errorbar(row.auc, yi, xerr=[[row.auc - row.auc_ci_low], [row.auc_ci_high - row.auc]], fmt="o", color=COLORS[feature], markersize=5, capsize=2)
    axes[0].axvline(0.5, color="grey", linestyle="--", linewidth=0.9)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([f"{label} | {feature}" for label, feature in rows], fontsize=8)
    axes[0].set_xlim(0.15, 1.0)
    axes[0].set_xlabel("AUC (bootstrap 95% CI)")
    axes[0].set_title("Independent group validation", loc="left")
    sns.despine(ax=axes[0])

    follow = pd.read_csv(RESULTS / "GSE83452_paired_followup.tsv", sep="\t")
    follow = follow[follow.feature.isin(FEATURES)].set_index("feature").reindex(FEATURES).reset_index()
    positions = np.arange(len(FEATURES))
    axes[1].axhline(0, color="grey", linewidth=0.8)
    axes[1].bar(positions, follow.median_delta, color=[COLORS[x] for x in FEATURES], alpha=0.82, width=0.68)
    axes[1].set_xticks(positions)
    axes[1].set_xticklabels(FEATURES, rotation=35, ha="right")
    axes[1].set_ylabel("Follow-up − baseline score")
    axes[1].set_title("GSE83452 paired follow-up\n60 patient pairs", loc="left")
    for pos, value, fdr in zip(positions, follow.median_delta, follow.fdr_within_analysis):
        axes[1].text(pos, value + 0.018, f"FDR={fdr:.2g}", ha="center", va="bottom", fontsize=7)
    sns.despine(ax=axes[1])

    paired = pd.read_csv(RESULTS / "GSE193066_HUn164_paired.tsv", sep="\t")
    paired = paired[paired.feature.isin(FEATURES)].set_index("feature").reindex(FEATURES).reset_index()
    axes[2].axvline(0, color="grey", linewidth=0.8)
    for pos, (feature, rho, fdr) in enumerate(zip(paired.feature, paired.delta_fibrosis_rho, paired.delta_fibrosis_fdr)):
        axes[2].errorbar(rho, pos, fmt="o", color=COLORS[feature], markersize=5)
        axes[2].text(rho + (0.025 if rho >= 0 else -0.025), pos, f"FDR={fdr:.2g}", va="center", ha="left" if rho >= 0 else "right", fontsize=7)
    axes[2].set_yticks(range(len(FEATURES)))
    axes[2].set_yticklabels(FEATURES)
    axes[2].set_xlim(-0.55, 0.55)
    axes[2].set_xlabel("Spearman rho: Δfeature vs Δfibrosis")
    axes[2].set_title("GSE193066 paired biopsies\n58 patient pairs", loc="left")
    sns.despine(ax=axes[2])

    fig.suptitle("Expanded public-cohort validation of the SLAMF8/SAMac–HSC axis", x=0.03, ha="left", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94], w_pad=2.3)
    for suffix, dpi in [("png", 320), ("pdf", 320)]:
        fig.savefig(FIGURES / f"017_expanded_public_validation.{suffix}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print("017_expanded_public_validation", flush=True)


if __name__ == "__main__":
    main()
