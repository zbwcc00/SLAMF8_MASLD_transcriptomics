from __future__ import annotations
import os

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import chi2, norm


PROJECT_ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
RESULTS = PROJECT_ROOT / "results"
FIGURES = PROJECT_ROOT / "figures"
DATASETS = ["GSE344087", "GSE298719", "GSE212837"]


def bh_adjust(values: pd.Series) -> pd.Series:
    order = np.argsort(values.to_numpy())
    ranked = values.to_numpy()[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1)
    return pd.Series(output, index=values.index)


def integrate_state_tests() -> tuple[pd.DataFrame, pd.DataFrame]:
    frames = []
    for dataset in DATASETS:
        frame = pd.read_csv(RESULTS / f"016_{dataset}_SLAMF8_high_paired_tests.tsv", sep="\t")
        frames.append(frame)
    combined = pd.concat(frames, ignore_index=True)
    rows = []
    for feature, group in combined.groupby("feature"):
        pvalues = group["p_value"].clip(lower=np.finfo(float).tiny).to_numpy()
        signs = np.sign(group["median_paired_difference"].to_numpy())
        weights = np.sqrt(group["n_samples"].to_numpy(dtype=float))
        zscores = signs * norm.isf(pvalues / 2)
        signed_z = float(np.sum(weights * zscores) / np.sqrt(np.sum(weights**2)))
        fisher_statistic = float(-2 * np.log(pvalues).sum())
        rows.append(
            {
                "feature": feature,
                "n_datasets": len(group),
                "n_samples_total": int(group["n_samples"].sum()),
                "positive_datasets": int((group["median_paired_difference"] > 0).sum()),
                "median_dataset_effect": float(group["median_paired_difference"].median()),
                "signed_stouffer_z": signed_z,
                "signed_stouffer_p": float(2 * norm.sf(abs(signed_z))),
                "fisher_statistic": fisher_statistic,
                "fisher_p": float(chi2.sf(fisher_statistic, 2 * len(group))),
            }
        )
    meta = pd.DataFrame(rows)
    meta["signed_stouffer_fdr"] = bh_adjust(meta["signed_stouffer_p"])
    meta["fisher_fdr"] = bh_adjust(meta["fisher_p"])
    combined.to_csv(RESULTS / "021_cross_dataset_SLAMF8_high_tests.tsv", sep="\t", index=False)
    meta.sort_values("signed_stouffer_p").to_csv(RESULTS / "022_cross_dataset_SLAMF8_high_meta.tsv", sep="\t", index=False)
    return combined, meta


def integrate_disease_tests() -> pd.DataFrame:
    frames = []
    for dataset in DATASETS:
        frame = pd.read_csv(RESULTS / f"014_{dataset}_patient_subtype_tests.tsv", sep="\t")
        frame = frame.loc[frame["subtype"] == "All_macrophage"].copy()
        frames.append(frame)
    combined = pd.concat(frames, ignore_index=True)
    rows = []
    for feature, group in combined.groupby("feature"):
        pvalues = group["p_value"].clip(lower=np.finfo(float).tiny).to_numpy()
        signs = np.sign(group["rank_biserial"].to_numpy())
        weights = np.sqrt((group["n_case"] + group["n_control"]).to_numpy(dtype=float))
        zscores = signs * norm.isf(pvalues / 2)
        signed_z = float(np.sum(weights * zscores) / np.sqrt(np.sum(weights**2)))
        rows.append(
            {
                "feature": feature,
                "n_datasets": len(group),
                "positive_datasets": int((group["rank_biserial"] > 0).sum()),
                "median_rank_biserial": float(group["rank_biserial"].median()),
                "signed_stouffer_z": signed_z,
                "signed_stouffer_p": float(2 * norm.sf(abs(signed_z))),
            }
        )
    output = pd.DataFrame(rows)
    output["signed_stouffer_fdr"] = bh_adjust(output["signed_stouffer_p"])
    combined.to_csv(RESULTS / "023_cross_dataset_disease_macrophage_tests.tsv", sep="\t", index=False)
    output.sort_values("signed_stouffer_p").to_csv(RESULTS / "024_cross_dataset_disease_macrophage_meta.tsv", sep="\t", index=False)
    return output


def save_figure(combined: pd.DataFrame, meta: pd.DataFrame) -> None:
    selected_features = [
        "gene_SLAMF8",
        "gene_TREM2",
        "gene_CD9",
        "gene_GPNMB",
        "gene_LGALS3",
        "module_FARG95",
        "module_Iron_homeostasis",
        "module_Lipid_peroxidation",
        "module_Senescence_SASP",
        "module_Inflammation",
    ]
    plot_data = combined.loc[combined["feature"].isin(selected_features)].copy()
    order = (
        meta.loc[meta["feature"].isin(selected_features)]
        .sort_values("signed_stouffer_z")["feature"]
        .tolist()
    )
    sns.set_theme(style="whitegrid", context="notebook")
    fig, ax = plt.subplots(figsize=(10, 7))
    sns.pointplot(
        data=plot_data,
        y="feature",
        x="median_paired_difference",
        hue="dataset",
        order=order,
        dodge=0.55,
        join=False,
        ax=ax,
    )
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Median within-patient/sample difference: SLAMF8-high minus low")
    ax.set_ylabel("")
    ax.set_title("Cross-dataset SLAMF8-high macrophage state")
    ax.legend(frameon=False, title="Dataset")
    fig.tight_layout()
    fig.savefig(FIGURES / "008_cross_dataset_SLAMF8_high_state.png", dpi=260)
    plt.close(fig)


def main() -> None:
    combined, meta = integrate_state_tests()
    integrate_disease_tests()
    save_figure(combined, meta)


if __name__ == "__main__":
    main()
