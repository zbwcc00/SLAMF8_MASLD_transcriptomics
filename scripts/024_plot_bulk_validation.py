#!/usr/bin/env python
import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc

PROJECT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
ROOT = PROJECT / "results" / "bulk_validation"
FIG = PROJECT / "figures"


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(ROOT / "GSE48452_bulk_scores.tsv", sep="\t")
    fibrosis = pd.to_numeric(data["fibrosis"], errors="coerce")
    data["fibrosis_group"] = fibrosis.map(lambda x: "F0" if x == 0 else "F≥1" if pd.notna(x) and x >= 1 else "NA")
    plot_data = data[data.fibrosis_group != "NA"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4))
    sns.boxplot(data=plot_data, x="fibrosis_group", y="SLAMF8", ax=axes[0], color="#8fb9d4")
    sns.stripplot(data=plot_data, x="fibrosis_group", y="SLAMF8", ax=axes[0], color="black", size=3, alpha=0.55)
    axes[0].set_title("GSE48452: SLAMF8 by fibrosis group")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Bulk expression")
    sns.regplot(data=data, x="SLAMF8", y="HSC_activation", ax=axes[1], scatter_kws={"s": 25, "alpha": 0.65}, line_kws={"color": "#b33b3b"})
    axes[1].set_title("SLAMF8 versus HSC activation")
    plt.tight_layout()
    plt.savefig(FIG / "012_bulk_SLAMF8_fibrosis_validation.png", dpi=300)
    plt.close()

    perf = pd.read_csv(ROOT / "bulk_composite_score_performance.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for axis, dataset, outcome in [(axes[0], "GSE48452", "fibrosis_ge1"), (axes[1], "GSE63067", "NASH_vs_other")]:
        score_file = ROOT / f"{dataset}_{outcome}_composite_scores.tsv"
        scores = pd.read_csv(score_file, sep="\t")
        for name in ["SLAMF8_only", "SLAMF8_FARG95_Iron", "SLAMF8_FARG95_Iron_HSC"]:
            fpr, tpr, _ = roc_curve(scores.outcome, scores[name])
            value = auc(fpr, tpr)
            axis.plot(fpr, tpr, label=f"{name} (AUC={value:.2f})")
        axis.plot([0, 1], [0, 1], linestyle="--", color="grey")
        axis.set_title(f"{dataset}: {outcome}")
        axis.set_xlabel("1 - specificity")
        axis.set_ylabel("sensitivity")
        axis.legend(fontsize=7, loc="lower right")
    plt.tight_layout()
    plt.savefig(FIG / "013_bulk_composite_score_ROC.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
