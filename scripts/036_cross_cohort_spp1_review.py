#!/usr/bin/env python
import os
"""Final donor-level cross-cohort review of the SLAMF8-SPP1-HSC hypothesis."""

from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import spearmanr, pearsonr
from statsmodels.stats.multitest import multipletests


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
DATA = ROOT / "data_processed"
OUT = ROOT / "results" / "communication" / "cross_cohort_spp1"
DATASETS = ["GSE344087", "GSE298719", "GSE212837"]
MACROPHAGE_SOURCE = {
    "GSE344087": ["TREM2_CD9_SAMac"],
    "GSE298719": ["GPNMB_TREM2_SAMac", "TREM2_CD9_SAMac"],
    "GSE212837": ["SPP1_GPNMB_macrophage", "TREM2_CD9_SPP1_SAMac"],
}
MACROPHAGE_EXCLUDE = r"cDC|pDC|Cycling|contaminant|ambiguous|low_quality|Neutrophil"
HSC_EXCLUDE = r"contaminant|low_quality|ambiguous"
MAC_GENES = ["SLAMF8", "SPP1", "GPNMB", "TREM2", "CD9", "PDGFC", "PDGFB", "TGFB1", "CXCL16"]
HSC_GENES = ["CD44", "ITGAV", "ITGB1", "ITGB5", "ITGA9", "PDGFRA", "PDGFRB", "ENG", "TGFBR1", "TGFBR2"]


def sample_col(obs: pd.DataFrame) -> str:
    return "patient_id" if "patient_id" in obs.columns else "sample_id"


def valid_mask(values: pd.Series, pattern: str) -> np.ndarray:
    return ~values.astype(str).str.contains(pattern, case=False, regex=True).to_numpy()


def pseudobulk(adata: ad.AnnData, sample_key: str, genes: list[str], label: str) -> pd.DataFrame:
    present = [gene for gene in genes if gene in adata.var_names]
    positions = {gene: adata.var_names.get_loc(gene) for gene in present}
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    rows = []
    for sample, indexes in adata.obs.groupby(sample_key, observed=True).indices.items():
        indexes = np.asarray(indexes)
        matrix = counts[indexes, :]
        total = np.asarray(matrix.sum(axis=1)).ravel().sum()
        if total <= 0:
            continue
        row = {"sample_id": str(sample), "n_cells": len(indexes), "label": label}
        for gene, position in positions.items():
            value = float(np.asarray(matrix[:, position].sum()).ravel()[0])
            row[gene] = np.log1p(1e6 * value / total)
        rows.append(row)
    return pd.DataFrame(rows)


def load_dataset(dataset: str) -> pd.DataFrame:
    mac = ad.read_h5ad(DATA / f"{dataset}_myeloid_reclustered.h5ad")
    hsc = ad.read_h5ad(DATA / f"{dataset}_hsc_reclustered.h5ad")
    mac_key, hsc_key = sample_col(mac.obs), sample_col(hsc.obs)
    mac.obs["_sample"] = mac.obs[mac_key].astype(str).values
    hsc.obs["_sample"] = hsc.obs[hsc_key].astype(str).values
    mac.obs["_source_state"] = mac.obs["subtype"].astype(str).isin(MACROPHAGE_SOURCE[dataset]).values
    mac.obs["_valid_mac"] = valid_mask(mac.obs["subtype"], MACROPHAGE_EXCLUDE)
    hsc.obs["_valid_hsc"] = valid_mask(hsc.obs["subtype"], HSC_EXCLUDE)
    # Use all non-contaminant myeloid cells for the broad state and the curated SAMac subset for the focused state.
    all_mac = mac[mac.obs._valid_mac].copy()
    source_mac = mac[mac.obs._valid_mac & mac.obs._source_state].copy()
    valid_hsc = hsc[hsc.obs._valid_hsc].copy()
    broad = pseudobulk(all_mac, "_sample", MAC_GENES, "All_macrophage")
    focused = pseudobulk(source_mac, "_sample", MAC_GENES, "SAMac_source")
    hsc_values = pseudobulk(valid_hsc, "_sample", HSC_GENES, "HSC")
    modules = valid_hsc.obs.groupby("_sample", observed=True)["module_HSC_activation"].mean().rename("HSC_activation")
    hsc_values = hsc_values.merge(modules, left_on="sample_id", right_index=True, how="left")
    meta_cols = ["_sample"] + [column for column in ["disease", "fibrosis"] if column in mac.obs]
    meta = mac.obs[meta_cols].drop_duplicates("_sample").rename(columns={"_sample": "sample_id"})
    values = broad.merge(focused, on="sample_id", how="outer", suffixes=("_mac", "_source"))
    values = values.merge(hsc_values, on="sample_id", how="inner", suffixes=("", "_hsc"))
    values = values.merge(meta, on="sample_id", how="left")
    values["dataset"] = dataset
    return values


def correlations(values: pd.DataFrame) -> pd.DataFrame:
    dataset = values.dataset.iloc[0]
    rows = []
    pairs = [
        ("SLAMF8_source", "SPP1_source", "SLAMF8_vs_SPP1_in_SAMac"),
        ("SLAMF8_mac", "HSC_activation", "SLAMF8_vs_HSC_activation"),
        ("SPP1_mac", "HSC_activation", "SPP1_vs_HSC_activation"),
        ("SPP1_source", "HSC_activation", "SPP1_source_vs_HSC_activation"),
    ]
    for receptor in HSC_GENES:
        pairs.append(("SPP1_source", receptor, f"SPP1_source_vs_HSC_{receptor}"))
    for x_name, y_name, feature in pairs:
        if x_name not in values or y_name not in values:
            continue
        for scope, subset in [("all_samples", values), ("disease_samples", values[values.get("disease", pd.Series(index=values.index, dtype=str)).astype(str).str.lower().isin(["masld", "nash", "fibrosis"])])]:
            pair = subset[[x_name, y_name]].dropna()
            if len(pair) < 5 or pair[x_name].nunique() < 2 or pair[y_name].nunique() < 2:
                continue
            rho, pvalue = spearmanr(pair[x_name], pair[y_name])
            rows.append({"dataset": dataset, "scope": scope, "feature": feature, "n": len(pair), "rho": rho, "p_value": pvalue})
    condition_name = "disease" if "disease" in values.columns else "fibrosis"
    for x_name, y_name, feature in pairs:
        if x_name not in values or y_name not in values or condition_name not in values:
            continue
        pair = values[[x_name, y_name, condition_name]].dropna()
        if len(pair) < 6 or pair[condition_name].nunique() < 2:
            continue
        design = pd.get_dummies(pair[condition_name].astype(str), drop_first=True, dtype=float)
        design.insert(0, "intercept", 1.0)
        x_rank = pd.Series(pair[x_name]).rank(method="average").to_numpy(dtype=float)
        y_rank = pd.Series(pair[y_name]).rank(method="average").to_numpy(dtype=float)
        design_values = design.to_numpy()
        x_resid = x_rank - design_values @ np.linalg.lstsq(design_values, x_rank, rcond=None)[0]
        y_resid = y_rank - design_values @ np.linalg.lstsq(design_values, y_rank, rcond=None)[0]
        rho, pvalue = pearsonr(x_resid, y_resid)
        rows.append({"dataset": dataset, "scope": "partial_condition", "feature": feature,
                     "n": len(pair), "rho": rho, "p_value": pvalue})
    return pd.DataFrame(rows)


def meta_correlations(table: pd.DataFrame, scope_name: str = "all_samples") -> pd.DataFrame:
    rows = []
    for feature, group in table.groupby("feature"):
        group = group[group.scope.eq(scope_name)].dropna(subset=["rho", "p_value"])
        if len(group) < 2:
            continue
        weights = np.maximum(group.n.to_numpy(dtype=float) - 3, 1)
        z = np.arctanh(np.clip(group.rho.to_numpy(dtype=float), -0.999999, 0.999999))
        pooled_z = np.sum(weights * z) / np.sum(weights)
        statistic = pooled_z * np.sqrt(np.sum(weights))
        from scipy.stats import norm
        pvalue = 2 * norm.sf(abs(statistic))
        rows.append({"feature": feature, "n_datasets": len(group), "n_total": int(group.n.sum()),
                     "positive_datasets": int((group.rho > 0).sum()), "pooled_rho": np.tanh(pooled_z),
                     "p_value": pvalue})
    output = pd.DataFrame(rows)
    if len(output):
        output["fdr"] = multipletests(output.p_value, method="fdr_bh")[1]
    return output


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    values = []
    tests = []
    for dataset in DATASETS:
        table = load_dataset(dataset)
        values.append(table)
        tests.append(correlations(table))
        table.to_csv(OUT / f"{dataset}_patient_level_SPP1_values.tsv", sep="\t", index=False)
    values = pd.concat(values, ignore_index=True)
    tests = pd.concat(tests, ignore_index=True)
    tests["fdr_within_scope"] = tests.groupby("scope").p_value.transform(lambda x: multipletests(x, method="fdr_bh")[1])
    tests.to_csv(OUT / "cross_cohort_SPP1_correlations.tsv", sep="\t", index=False)
    meta_correlations(tests).to_csv(OUT / "cross_cohort_SPP1_meta_correlations.tsv", sep="\t", index=False)
    meta_correlations(tests, "partial_condition").to_csv(OUT / "cross_cohort_SPP1_partial_condition_meta.tsv", sep="\t", index=False)
    values.to_csv(OUT / "all_cohorts_patient_level_SPP1_values.tsv", sep="\t", index=False)
    print(tests.to_string(index=False), flush=True)
    print(meta_correlations(tests).to_string(index=False), flush=True)
    print(meta_correlations(tests, "partial_condition").to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
