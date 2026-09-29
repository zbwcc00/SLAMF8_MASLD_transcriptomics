#!/usr/bin/env python
"""Continuous SLAMF8, threshold sensitivity, random-effects and LOCO audits."""
from __future__ import annotations
import os

from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import spearmanr, norm, wilcoxon
from statsmodels.stats.multitest import multipletests

ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
FARG = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1])) / "data" / "reference" / "FARG95_gene_set_provenance.tsv"
OUT = ROOT / "results" / "revision_p0" / "continuous_meta"
OUT.mkdir(parents=True, exist_ok=True)
DATASETS = ["GSE344087", "GSE298719", "GSE212837"]
ANCHORS = {"SLAMF8", "SPP1", "GPNMB", "TREM2", "CD9"}
GENESETS = {
    "FARG95_SLAMF8_excluded": None,
    "FARG95_anchor_excluded": None,
    "Iron_homeostasis": ["TFRC", "SLC11A2", "FTH1", "FTL", "SLC40A1", "NCOA4", "HMOX1", "STEAP3", "PCBP1", "PCBP2", "HAMP"],
    "Senescence_SASP": ["CDKN1A", "CDKN2A", "GLB1", "SERPINE1", "GDF15", "IL6", "IL1B", "CXCL8", "MMP3", "MMP9", "CCL2", "TGFB1", "IGFBP7"],
}


def read_raw(adata: ad.AnnData, genes: list[str]) -> pd.DataFrame:
    genes = [g for g in genes if g in adata.raw.var_names]
    x = adata.raw[:, genes].X
    if hasattr(x, "toarray"):
        x = x.toarray()
    return pd.DataFrame(np.asarray(x, dtype=float), index=adata.obs_names, columns=genes)


def valid_mask(obs: pd.DataFrame) -> np.ndarray:
    return ~obs["subtype"].astype(str).str.contains(r"cDC|pDC|_DC|Cycling|ambiguous|contaminant|Neutrophil|low_quality", case=False, regex=True).to_numpy()


def dataset_values(dataset: str, fargs: list[str]) -> pd.DataFrame:
    adata = ad.read_h5ad(ROOT / "data_processed" / f"{dataset}_myeloid_reclustered.h5ad")
    adata = adata[valid_mask(adata.obs)].copy()
    sample_col = "patient_id" if "patient_id" in adata.obs.columns else "gsm"
    gene_sets = {
        "FARG95_SLAMF8_excluded": [g for g in fargs if g != "SLAMF8"],
        "FARG95_anchor_excluded": [g for g in fargs if g not in ANCHORS],
        "Iron_homeostasis": GENESETS["Iron_homeostasis"],
        "Senescence_SASP": GENESETS["Senescence_SASP"],
    }
    genes = sorted(set(["SLAMF8"] + sum(gene_sets.values(), [])))
    expr = read_raw(adata, genes)
    obs = adata.obs.copy()
    obs["_sample"] = obs[sample_col].astype(str).to_numpy()
    obs["SLAMF8"] = expr["SLAMF8"].to_numpy()
    for name, members in gene_sets.items():
        members = [g for g in members if g in expr.columns]
        obs[name] = expr[members].mean(axis=1) if members else np.nan
    condition_col = {"GSE344087": "fibrosis", "GSE298719": "disease", "GSE212837": "disease"}[dataset]
    metadata = obs.groupby("_sample", observed=True).agg(
        SLAMF8=("SLAMF8", "mean"),
        **{name: (name, "mean") for name in gene_sets},
        n_cells=("SLAMF8", "size"),
    ).reset_index().rename(columns={"_sample": "sample"})
    conditions = obs[["_sample", condition_col]].drop_duplicates("_sample").rename(columns={"_sample": "sample"})
    return metadata.merge(conditions, on="sample", how="left").assign(dataset=dataset)


def continuous_tests(values: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for dataset, frame in values.groupby("dataset", sort=False):
        for feature in GENESETS:
            pair = frame[["SLAMF8", feature]].dropna()
            if len(pair) < 5 or pair["SLAMF8"].nunique() < 2 or pair[feature].nunique() < 2:
                continue
            rho, p = spearmanr(pair["SLAMF8"], pair[feature])
            rows.append({"dataset": dataset, "feature": feature, "n_samples": len(pair), "rho": rho, "p_value": p})
    tests = pd.DataFrame(rows)
    if len(tests):
        tests["fdr_within_dataset"] = tests.groupby("dataset")["p_value"].transform(lambda x: multipletests(x, method="fdr_bh")[1])
    return tests


def random_effects(group: pd.DataFrame) -> dict[str, float]:
    group = group.dropna(subset=["rho", "n_samples"])
    k = len(group)
    if k < 2:
        return {"n_datasets": k, "pooled_rho_random": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p_value": np.nan, "tau2": np.nan, "I2": np.nan}
    rho = np.clip(group.rho.to_numpy(float), -0.999999, 0.999999)
    z = np.arctanh(rho)
    var = 1.0 / np.maximum(group.n_samples.to_numpy(float) - 3.0, 1.0)
    w = 1.0 / var
    fixed = np.sum(w * z) / np.sum(w)
    q = np.sum(w * (z - fixed) ** 2)
    df = k - 1
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    wr = 1.0 / (var + tau2)
    pooled = np.sum(wr * z) / np.sum(wr)
    se = np.sqrt(1.0 / np.sum(wr))
    return {
        "n_datasets": k,
        "pooled_rho_random": float(np.tanh(pooled)),
        "ci_low": float(np.tanh(pooled - 1.96 * se)),
        "ci_high": float(np.tanh(pooled + 1.96 * se)),
        "p_value": float(2 * norm.sf(abs(pooled / se))),
        "tau2": float(tau2),
        "I2": float(max(0.0, (q - df) / q) if q > 0 else 0.0),
    }


def meta_and_loco(tests: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, loco = [], []
    for feature, group in tests.groupby("feature", sort=False):
        result = random_effects(group)
        result["feature"] = feature
        result["positive_datasets"] = int((group.rho > 0).sum())
        rows.append(result)
        for excluded in group.dataset:
            remain = group[group.dataset != excluded]
            sub = random_effects(remain)
            sub.update({"feature": feature, "excluded_dataset": excluded, "positive_datasets": int((remain.rho > 0).sum())})
            loco.append(sub)
    meta = pd.DataFrame(rows)
    if len(meta):
        meta["fdr"] = multipletests(meta.p_value.fillna(1), method="fdr_bh")[1]
    return meta.sort_values("p_value"), pd.DataFrame(loco)


def threshold_sensitivity(fargs: list[str]) -> pd.DataFrame:
    rows = []
    for dataset in DATASETS:
        adata = ad.read_h5ad(ROOT / "data_processed" / f"{dataset}_myeloid_reclustered.h5ad")
        adata = adata[valid_mask(adata.obs)].copy()
        sample_col = "patient_id" if "patient_id" in adata.obs.columns else "gsm"
        expr = read_raw(adata, sorted(set(["SLAMF8"] + [g for g in fargs if g != "SLAMF8"])))
        obs = adata.obs.copy()
        obs["_sample"] = obs[sample_col].astype(str).to_numpy()
        obs["SLAMF8"] = expr["SLAMF8"].to_numpy()
        members = [g for g in fargs if g not in ANCHORS and g in expr.columns]
        obs["FARG95_anchor_excluded"] = expr[members].mean(axis=1)
        min_high = 5 if dataset == "GSE212837" else 20
        for quantile in (0.70, 0.75, 0.80):
            deltas = []
            for sample, group in obs.groupby("_sample", observed=True):
                cut = group.SLAMF8.quantile(quantile)
                high, low = group[group.SLAMF8 > cut], group[group.SLAMF8 <= cut]
                if len(high) >= min_high and len(low) >= 20:
                    deltas.append(high.FARG95_anchor_excluded.mean() - low.FARG95_anchor_excluded.mean())
            if len(deltas) >= 5 and not np.allclose(deltas, 0):
                stat, p = wilcoxon(deltas)
            else:
                stat, p = np.nan, np.nan
            rows.append({"dataset": dataset, "quantile": quantile, "n_samples": len(deltas), "median_delta": np.median(deltas) if deltas else np.nan, "positive_fraction": np.mean(np.asarray(deltas) > 0) if deltas else np.nan, "wilcoxon_statistic": stat, "p_value": p})
    out = pd.DataFrame(rows)
    if len(out):
        out["fdr_all_thresholds"] = multipletests(out.p_value.fillna(1), method="fdr_bh")[1]
    return out


def main() -> None:
    farg_table = pd.read_csv(FARG, sep="\t")
    farg_col = "gene" if "gene" in farg_table.columns else farg_table.columns[0]
    fargs = farg_table[farg_col].dropna().astype(str).str.upper().drop_duplicates().tolist()
    values = pd.concat([dataset_values(ds, fargs) for ds in DATASETS], ignore_index=True)
    values.to_csv(OUT / "patient_level_continuous_scores.tsv", sep="\t", index=False)
    tests = continuous_tests(values)
    tests.to_csv(OUT / "continuous_slamf8_correlations.tsv", sep="\t", index=False)
    meta, loco = meta_and_loco(tests)
    meta.to_csv(OUT / "continuous_random_effects_meta.tsv", sep="\t", index=False)
    loco.to_csv(OUT / "continuous_leave_one_cohort_out.tsv", sep="\t", index=False)
    threshold_sensitivity(fargs).to_csv(OUT / "threshold_sensitivity.tsv", sep="\t", index=False)
    (OUT / "README.md").write_text(
        "# Continuous SLAMF8 and heterogeneity audit\n\n"
        "Continuous sample-level SLAMF8 associations use FARG95 scores excluding SLAMF8 and anchor genes. "
        "Random-effects summaries use Fisher-z effects with DerSimonian–Laird tau2; LOCO rows show whether conclusions depend on one cohort. "
        "Threshold sensitivity repeats the paired contrast at the 70th, 75th and 80th percentiles.\n",
        encoding="utf-8",
    )
    print(meta.to_string(index=False))
    print(threshold_sensitivity(fargs).to_string(index=False))


if __name__ == "__main__":
    main()
