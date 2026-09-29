#!/usr/bin/env python
import os
"""Plot independent NAFLD liver-biopsy transcriptome validation."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
TABLES = ROOT / "results" / "bulk_validation"
FIGURES = ROOT / "figures"
COHORTS = ["GSE162694", "GSE130970"]
FEATURES = ["SLAMF8", "GPNMB", "FARG95", "Iron_homeostasis", "HSC_activation"]


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    results = pd.read_csv(TABLES / "full_bulk_fibrosis_validation.tsv", sep="\t")
    sns.set_theme(style="ticks", context="paper", font_scale=1.15)
    fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.8), gridspec_kw={"width_ratios": [1, 1, 1.35]})
    for ax, accession in zip(axes[:2], COHORTS):
        data = pd.read_csv(TABLES / f"{accession}_full_bulk_scores.tsv", sep="\t")
        data = data.loc[~data.normal_histology].copy()
        data["fibrosis"] = data.fibrosis.astype(int)
        orders = sorted(data.fibrosis.unique())
        sns.boxplot(data=data, x="fibrosis", y="SLAMF8", order=orders, ax=ax,
                    color="#92bad5", width=0.65, fliersize=0, linewidth=1)
        sns.stripplot(data=data, x="fibrosis", y="SLAMF8", order=orders, ax=ax,
                      color="#173e61", size=3, alpha=0.56, jitter=0.18)
        result = results.loc[(results.dataset == accession) &
                             (results.scope == "NAFLD_only") &
                             (results.feature == "SLAMF8")].iloc[0]
        ax.set_title(f"{accession}  |  n={int(result.n)}\n"
                     f"Spearman rho={result.rho:.3f}, P={result.p_value:.2g}", loc="left")
        ax.set_xlabel("Fibrosis stage")
        ax.set_ylabel("SLAMF8  log2(CPM + 1)" if accession == COHORTS[0]
                      else "SLAMF8  log2(TPM + 1)")
        sns.despine(ax=ax)

    selected = results[(results.scope == "NAFLD_only") &
                       results.feature.isin(FEATURES)].pivot(index="feature", columns="dataset", values="rho")
    selected = selected.reindex(index=FEATURES, columns=COHORTS)
    sig = results[(results.scope == "NAFLD_only") &
                  results.feature.isin(FEATURES)].pivot(index="feature", columns="dataset", values="fdr_within_dataset_scope")
    sig = sig.reindex(index=FEATURES, columns=COHORTS)
    labels = np.array([[f"{selected.loc[feature, ds]:+.2f}" +
                        ("*" if sig.loc[feature, ds] < 0.05 else "")
                        for ds in COHORTS] for feature in FEATURES])
    sns.heatmap(selected, ax=axes[2], cmap="RdBu_r", center=0, vmin=-0.65, vmax=0.65,
                annot=labels, fmt="", cbar_kws={"label": "Spearman rho", "shrink": 0.8},
                linewidths=1, linecolor="white")
    axes[2].set_title("Within-NAFLD fibrosis association\n* BH-FDR < 0.05", loc="left")
    axes[2].set_xlabel("")
    axes[2].set_ylabel("")
    axes[2].tick_params(axis="y", rotation=0)
    axes[2].tick_params(axis="x", rotation=20)
    fig.tight_layout(w_pad=2.0)
    for suffix in ["png", "pdf"]:
        fig.savefig(FIGURES / f"014_full_bulk_independent_validation.{suffix}", dpi=320, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
