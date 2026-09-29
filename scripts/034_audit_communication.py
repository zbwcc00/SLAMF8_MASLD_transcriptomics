#!/usr/bin/env python
import os
"""Audit macrophage-to-HSC communication without treating cells as replicates."""

from pathlib import Path
import gzip
import re

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import mannwhitneyu, spearmanr


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
COMM = ROOT / "results" / "communication"
DATA = ROOT / "data_processed"
OUT = COMM / "audit"
CANDIDATES = {
    "SPP1": ["CD44", "ITGAV_ITGB1", "ITGAV_ITGB5", "ITGA9_ITGB1"],
    "PDGFC": ["PDGFRA", "PDGFRB"],
    "PDGFB": ["LRP1", "PDGFRA_PDGFRB", "ITGAV"],
    "TGFB1": ["ENG", "TGFBR1_TGFBR2", "ACVR1_TGFBR2", "ITGB5"],
    "CXCL16": ["CXCR6"],
}


def candidate_mask(value: str, candidates: list[str]) -> bool:
    return str(value) in candidates


def audit_inferred_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    cellchat_rows, liana_rows = [], []
    for condition in ["Healthy", "MASLD"]:
        cellchat = pd.read_csv(COMM / f"019_GSE298719_{condition}_CellChat_macrophage_to_HSC.tsv", sep="\t")
        liana = pd.read_csv(COMM / f"018_GSE298719_{condition}_LIANA_macrophage_to_HSC.tsv", sep="\t")
        liana["active_consensus"] = liana["aggregate_rank"] < 0.5
        for ligand, receptors in CANDIDATES.items():
            cc = cellchat[cellchat.ligand.eq(ligand)].copy()
            cc = cc[cc.receptor.isin(receptors)]
            cellchat_rows.append({
                "condition": condition, "ligand": ligand, "method": "CellChat",
                "n_pairs": len(cc), "n_significant": int((cc.pval < 0.05).sum()),
                "significant_fraction": float((cc.pval < 0.05).mean()) if len(cc) else np.nan,
                "median_score": cc.prob.median() if len(cc) else np.nan,
                "max_score": cc.prob.max() if len(cc) else np.nan,
                "sources": ";".join(sorted(cc.source.unique())) if len(cc) else "",
                "targets": ";".join(sorted(cc.target.unique())) if len(cc) else "",
            })
            li = liana[liana["ligand.complex"].eq(ligand)].copy()
            li = li[li["receptor.complex"].isin(receptors)]
            li_active = li[li.active_consensus]
            liana_rows.append({
                "condition": condition, "ligand": ligand, "method": "LIANA",
                "n_pairs": len(li_active), "n_inferred_rows": len(li),
                "n_active": len(li_active),
                "median_aggregate_rank": li_active.aggregate_rank.median() if len(li_active) else np.nan,
                "median_mean_rank": li_active.mean_rank.median() if len(li_active) else np.nan,
                "sources": ";".join(sorted(li_active.source.unique())) if len(li_active) else "",
                "targets": ";".join(sorted(li_active.target.unique())) if len(li_active) else "",
            })
    return pd.DataFrame(cellchat_rows), pd.DataFrame(liana_rows)


def load_bundle() -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = pd.read_csv(DATA / "GSE298719_communication_cell_metadata.tsv.gz", sep="\t")
    with gzip.open(DATA / "GSE298719_communication_features.tsv.gz", "rt") as handle:
        features = pd.read_csv(handle, sep="\t")
    with gzip.open(DATA / "GSE298719_communication_counts.mtx.gz", "rb") as handle:
        counts = mmread(handle).tocsr()
    features["gene"] = features["var_name"].astype(str).str.upper()
    gene_positions = {}
    for index, gene in enumerate(features.gene):
        gene_positions.setdefault(gene, index)
    return metadata, counts, gene_positions


def expression_by_sample() -> pd.DataFrame:
    metadata, counts, positions = load_bundle()
    genes = sorted(set(CANDIDATES) | {receptor for receptors in CANDIDATES.values() for receptor in receptors for receptor in receptor.split("_")} | {"SLAMF8"})
    genes = [gene for gene in genes if gene in positions]
    rows = []
    for (disease, sample_id), sample in metadata.groupby(["disease", "sample_id"], observed=True):
        row = {"condition": disease, "sample_id": sample_id}
        for compartment, prefix in [("Macrophage", "mac"), ("HSC", "hsc")]:
            indices = sample.index[sample.compartment.eq(compartment)].to_numpy()
            if len(indices) < 20:
                continue
            library = np.asarray(counts[:, indices].sum(axis=0)).ravel()
            for gene in genes:
                values = np.asarray(counts[positions[gene], indices].toarray()).ravel()
                row[f"{prefix}_{gene}"] = np.log2(1 + np.mean(values / np.maximum(library, 1) * 1e4))
        rows.append(row)
    return pd.DataFrame(rows)


def subtype_expression() -> pd.DataFrame:
    metadata, counts, positions = load_bundle()
    genes = [gene for gene in ["SLAMF8", "SPP1", "PDGFC", "PDGFB", "TGFB1", "CXCL16", "GPNMB"] if gene in positions]
    rows = []
    source_groups = ["Mac_GPNMB_TREM2_SAMac", "Mac_TREM2_CD9_SAMac", "Mac_FCN1_VCAN_monocyte"]
    for (disease, sample_id, group), sample in metadata.groupby(["disease", "sample_id", "cell_group"], observed=True):
        if group not in source_groups or len(sample) < 10:
            continue
        indices = sample.index.to_numpy()
        library = np.asarray(counts[:, indices].sum(axis=0)).ravel()
        row = {"condition": disease, "sample_id": sample_id, "cell_group": group, "n_cells": len(indices)}
        for gene in genes:
            values = np.asarray(counts[positions[gene], indices].toarray()).ravel()
            row[gene] = np.log2(1 + np.mean(values / np.maximum(library, 1) * 1e4))
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_expression(table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene in sorted(set(CANDIDATES) | {receptor for receptors in CANDIDATES.values() for receptor in receptors for receptor in receptor.split("_")} | {"SLAMF8"}):
        for prefix in ["mac", "hsc"]:
            col = f"{prefix}_{gene}"
            if col not in table:
                continue
            groups = [table.loc[table.condition.eq(condition), col].dropna().to_numpy() for condition in ["Healthy", "MASLD"]]
            if not all(len(values) >= 3 for values in groups):
                continue
            stat, pvalue = mannwhitneyu(groups[1], groups[0], alternative="two-sided")
            rows.append({"feature": col, "healthy_n": len(groups[0]), "masld_n": len(groups[1]),
                         "healthy_median": np.median(groups[0]), "masld_median": np.median(groups[1]),
                         "masld_minus_healthy": np.median(groups[1]) - np.median(groups[0]), "p_value": pvalue})
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cellchat, liana = audit_inferred_tables()
    cellchat.to_csv(OUT / "cellchat_candidate_audit.tsv", sep="\t", index=False)
    liana.to_csv(OUT / "liana_candidate_audit.tsv", sep="\t", index=False)
    expression = expression_by_sample()
    expression.to_csv(OUT / "patient_level_ligand_receptor_expression.tsv", sep="\t", index=False)
    summarize_expression(expression).to_csv(OUT / "patient_level_ligand_receptor_expression_tests.tsv", sep="\t", index=False)
    subtype = subtype_expression()
    subtype.to_csv(OUT / "patient_level_macrophage_subtype_expression.tsv", sep="\t", index=False)
    association_rows = []
    for group, data in subtype.groupby("cell_group", observed=True):
        masld = data[data.condition.eq("MASLD")]
        if len(masld) >= 4 and masld.SLAMF8.nunique() > 1:
            for gene in ["SPP1", "PDGFC", "PDGFB", "TGFB1", "CXCL16"]:
                if gene in masld and masld[gene].nunique() > 1:
                    rho, pvalue = spearmanr(masld.SLAMF8, masld[gene])
                    association_rows.append({"condition": "MASLD", "cell_group": group, "feature": gene,
                                             "n": len(masld), "spearman_rho": rho, "p_value": pvalue})
    pd.DataFrame(association_rows).to_csv(OUT / "patient_level_SLAMF8_ligand_associations.tsv", sep="\t", index=False)
    print(cellchat.to_string(index=False), flush=True)
    print(liana.to_string(index=False), flush=True)
    print(summarize_expression(expression).to_string(index=False), flush=True)
    print(pd.DataFrame(association_rows).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
