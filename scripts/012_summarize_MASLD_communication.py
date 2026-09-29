#!/usr/bin/env python
import os
"""Filter and compare Healthy versus MASLD macrophage-HSC communication."""

from pathlib import Path
import re
import pandas as pd


PROJECT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
COMM = PROJECT / "results" / "communication"


def pathway(ligand: str, receptor: str = "") -> str:
    text = f"{ligand}_{receptor}".upper()
    for key, label in [("TGFB", "TGFb"), ("PDGF", "PDGF"), ("TNF", "TNF"), ("IL1", "IL1"), ("NOTCH", "NOTCH"), ("JAG", "NOTCH"), ("DLL", "NOTCH"), ("SPP1", "SPP1")]:
        if key in text:
            return label
    return "Other"


def liana_table(condition: str) -> pd.DataFrame:
    path = COMM / f"017_GSE298719_{condition}_LIANA_consensus.tsv"
    data = pd.read_csv(path, sep="\t")
    data = data[data.source.astype(str).str.startswith("Mac_") & data.target.astype(str).str.startswith("HSC_")].copy()
    ligand = "ligand.complex" if "ligand.complex" in data.columns else "ligand_complex"
    receptor = "receptor.complex" if "receptor.complex" in data.columns else "receptor_complex"
    data["focus_pathway"] = [pathway(a, b) for a, b in zip(data[ligand], data[receptor])]
    data["condition"] = condition
    return data


def cellchat_table(condition: str) -> pd.DataFrame:
    path = COMM / f"019_GSE298719_{condition}_CellChat_macrophage_to_HSC.tsv"
    data = pd.read_csv(path, sep="\t")
    data["condition"] = condition
    data["focus_pathway"] = [pathway(a, b) for a, b in zip(data.ligand, data.receptor)]
    return data


def main() -> None:
    healthy = liana_table("Healthy")
    healthy.to_csv(COMM / "018_GSE298719_Healthy_LIANA_macrophage_to_HSC.tsv", sep="\t", index=False)
    masld = liana_table("MASLD")
    both = pd.concat([healthy, masld], ignore_index=True)
    rank_col = "aggregate_rank" if "aggregate_rank" in both.columns else "mean_rank"
    both["active_consensus"] = both[rank_col] < 0.5
    summary = (both.groupby(["condition", "focus_pathway"], as_index=False)
               .agg(interactions=(rank_col, "size"), active_fraction=("active_consensus", "mean"), median_rank=(rank_col, "median"),
                    median_mean_rank=("mean_rank", "median"), median_natmi=("natmi.prod_weight", "median"),
                    median_cellphonedb=("cellphonedb.lr.mean", "median"), median_sca=("sca.LRscore", "median")))
    pivot = summary.pivot(index="focus_pathway", columns="condition")
    for metric in ["interactions", "median_rank", "median_mean_rank", "median_natmi", "median_cellphonedb", "median_sca"]:
        if (metric, "Healthy") in pivot and (metric, "MASLD") in pivot:
            summary.loc[:, f"MASLD_minus_Healthy_{metric}"] = summary[metric].where(summary.condition.eq("MASLD"), pd.NA)
    wide = summary.pivot(index="focus_pathway", columns="condition")
    rows = []
    for label in sorted(both.focus_pathway.unique()):
        row = {"focus_pathway": label}
        for metric in ["interactions", "active_fraction", "median_rank", "median_mean_rank", "median_natmi", "median_cellphonedb", "median_sca"]:
            for cond in ["Healthy", "MASLD"]:
                vals = summary.loc[(summary.focus_pathway == label) & (summary.condition == cond), metric]
                row[f"{cond}_{metric}"] = float(vals.iloc[0]) if len(vals) else float("nan")
        row["MASLD_minus_Healthy_interactions"] = row.get("MASLD_interactions", float("nan")) - row.get("Healthy_interactions", float("nan"))
        row["MASLD_minus_Healthy_active_fraction"] = row.get("MASLD_active_fraction", float("nan")) - row.get("Healthy_active_fraction", float("nan"))
        row["MASLD_minus_Healthy_median_rank"] = row.get("MASLD_median_rank", float("nan")) - row.get("Healthy_median_rank", float("nan"))
        rows.append(row)
    pd.DataFrame(rows).to_csv(COMM / "021_GSE298719_LIANA_Healthy_vs_MASLD_pathway_summary.tsv", sep="\t", index=False)

    cell = pd.concat([cellchat_table("Healthy"), cellchat_table("MASLD")], ignore_index=True)
    cell_summary = (cell.groupby(["condition", "focus_pathway"], as_index=False)
                    .agg(interactions=("prob", "size"), median_probability=("prob", "median"),
                         max_probability=("prob", "max"), significant_fraction=("pval", lambda x: float((x < 0.05).mean()))))
    cell_rows = []
    for label in sorted(cell.focus_pathway.unique()):
        row = {"focus_pathway": label}
        for metric in ["interactions", "median_probability", "max_probability", "significant_fraction"]:
            for cond in ["Healthy", "MASLD"]:
                vals = cell_summary.loc[(cell_summary.focus_pathway == label) & (cell_summary.condition == cond), metric]
                row[f"{cond}_{metric}"] = float(vals.iloc[0]) if len(vals) else float("nan")
        row["MASLD_minus_Healthy_interactions"] = row.get("MASLD_interactions", float("nan")) - row.get("Healthy_interactions", float("nan"))
        row["MASLD_minus_Healthy_median_probability"] = row.get("MASLD_median_probability", float("nan")) - row.get("Healthy_median_probability", float("nan"))
        cell_rows.append(row)
    pd.DataFrame(cell_rows).to_csv(COMM / "022_GSE298719_CellChat_Healthy_vs_MASLD_pathway_summary.tsv", sep="\t", index=False)

    top = (cell.sort_values(["condition", "prob"], ascending=[True, False])
           .groupby("condition", as_index=False, group_keys=False).head(30))
    top.to_csv(COMM / "023_GSE298719_CellChat_top_macrophage_HSC_interactions.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
