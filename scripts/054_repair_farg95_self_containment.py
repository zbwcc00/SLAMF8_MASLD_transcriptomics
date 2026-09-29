#!/usr/bin/env python
"""Repair the SLAMF8/FARG95 self-containment issue without overwriting old results."""
from __future__ import annotations
import os

from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import chi2, norm, wilcoxon
from statsmodels.stats.multitest import multipletests

ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
FARG = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1])) / "data" / "reference" / "FARG95_gene_set_provenance.tsv"
OUT = ROOT / "results" / "revision_p0" / "farg95_self_containment"
OUT.mkdir(parents=True, exist_ok=True)
DATASETS = ["GSE344087", "GSE298719", "GSE212837"]
ANCHORS = {"SLAMF8", "SPP1", "GPNMB", "TREM2", "CD9"}
COMPARISONS = {
    "GSE344087": ("fibrosis", "fibrosis", "no fibrosis"),
    "GSE298719": ("disease", "MASLD", "Healthy"),
    "GSE212837": ("disease", "NASH", "Control"),
}


def bh(values: pd.Series) -> pd.Series:
    return pd.Series(multipletests(values.to_numpy(dtype=float), method="fdr_bh")[1], index=values.index)


def read_matrix(adata: ad.AnnData, genes: list[str]) -> pd.DataFrame:
    genes = [g for g in genes if g in adata.raw.var_names]
    x = adata.raw[:, genes].X
    if hasattr(x, "toarray"):
        x = x.toarray()
    return pd.DataFrame(np.asarray(x, dtype=float), index=adata.obs_names, columns=genes)


def valid_mask(obs: pd.DataFrame) -> np.ndarray:
    excluded = r"cDC|pDC|_DC|Cycling|ambiguous|contaminant|Neutrophil|low_quality"
    return ~obs["subtype"].astype(str).str.contains(excluded, case=False, regex=True).to_numpy()


def run_dataset(dataset: str, fargs: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    adata = ad.read_h5ad(ROOT / "data_processed" / f"{dataset}_myeloid_reclustered.h5ad")
    adata = adata[valid_mask(adata.obs)].copy()
    sample_col = "patient_id" if "patient_id" in adata.obs.columns else "gsm"
    genes = sorted(set(fargs) | {"SLAMF8"})
    expression = read_matrix(adata, genes)
    if "SLAMF8" not in expression:
        raise RuntimeError(f"{dataset}: SLAMF8 is absent from normalized raw layer")
    score_sets = {
        "FARG95_original": [g for g in fargs if g in expression.columns],
        "FARG95_SLAMF8_excluded": [g for g in fargs if g not in {"SLAMF8"} and g in expression.columns],
        "FARG95_anchor_excluded": [g for g in fargs if g not in ANCHORS and g in expression.columns],
    }
    scores = pd.DataFrame(index=expression.index)
    for name, members in score_sets.items():
        scores[name] = expression[members].mean(axis=1) if members else np.nan
    scores["SLAMF8_expression"] = expression["SLAMF8"]
    obs = adata.obs.copy()
    obs["_sample"] = obs[sample_col].astype(str).to_numpy()
    obs["SLAMF8_expression"] = scores["SLAMF8_expression"].to_numpy()
    for col in scores.columns:
        if col != "SLAMF8_expression":
            obs[col] = scores[col].to_numpy()
    cells = obs[["_sample", "SLAMF8_expression"] + list(score_sets)].copy()
    cells.insert(0, "dataset", dataset)
    cells.insert(1, "barcode", obs.index.astype(str))
    cells.to_csv(OUT / f"{dataset}_cell_scores.tsv.gz", sep="\t", index=False, compression="gzip")

    min_high = 5 if dataset == "GSE212837" else 20
    rows = []
    for sample, group in obs.groupby("_sample", observed=True):
        threshold = float(group["SLAMF8_expression"].quantile(0.75))
        high = group[group["SLAMF8_expression"] > threshold]
        low = group[group["SLAMF8_expression"] <= threshold]
        if len(high) < min_high or len(low) < 20:
            continue
        for feature in score_sets:
            rows.append({
                "dataset": dataset,
                "sample": sample,
                "feature": feature,
                "n_high": len(high),
                "n_low": len(low),
                "threshold": threshold,
                "mean_high": high[feature].mean(),
                "mean_low": low[feature].mean(),
                "paired_difference": high[feature].mean() - low[feature].mean(),
                "n_genes_present": len(score_sets[feature]),
            })
    differences = pd.DataFrame(rows)
    differences.to_csv(OUT / f"{dataset}_paired_differences.tsv.gz", sep="\t", index=False, compression="gzip")
    tests = []
    for feature, group in differences.groupby("feature", sort=False):
        delta = group["paired_difference"].dropna().to_numpy(dtype=float)
        if len(delta) < 5 or np.allclose(delta, 0):
            continue
        stat, p = wilcoxon(delta, alternative="two-sided", zero_method="wilcox")
        tests.append({
            "dataset": dataset,
            "feature": feature,
            "n_samples": len(delta),
            "median_paired_difference": float(np.median(delta)),
            "mean_paired_difference": float(np.mean(delta)),
            "positive_fraction": float(np.mean(delta > 0)),
            "wilcoxon_statistic": float(stat),
            "p_value": float(p),
        })
    tests = pd.DataFrame(tests)
    if len(tests):
        tests["fdr_within_dataset"] = bh(tests["p_value"])
    tests.to_csv(OUT / f"{dataset}_paired_tests.tsv", sep="\t", index=False)
    return differences, tests


def meta_analysis(tests: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for feature, group in tests.groupby("feature", sort=False):
        p = np.clip(group["p_value"].to_numpy(float), np.finfo(float).tiny, 1.0)
        signs = np.sign(group["median_paired_difference"].to_numpy(float))
        weights = np.sqrt(group["n_samples"].to_numpy(float))
        z = signs * norm.isf(p / 2.0)
        z_stouffer = float(np.sum(weights * z) / np.sqrt(np.sum(weights ** 2)))
        fisher_stat = float(-2 * np.log(p).sum())
        rows.append({
            "feature": feature,
            "n_datasets": len(group),
            "n_samples_total": int(group["n_samples"].sum()),
            "positive_datasets": int((group["median_paired_difference"] > 0).sum()),
            "median_dataset_effect": float(group["median_paired_difference"].median()),
            "signed_stouffer_z": z_stouffer,
            "signed_stouffer_p": float(2 * norm.sf(abs(z_stouffer))),
            "fisher_p": float(chi2.sf(fisher_stat, 2 * len(group))),
        })
    meta = pd.DataFrame(rows)
    if len(meta):
        meta["signed_stouffer_fdr"] = bh(meta["signed_stouffer_p"])
        meta["fisher_fdr"] = bh(meta["fisher_p"])
    return meta.sort_values("signed_stouffer_p")


def main() -> None:
    farg_table = pd.read_csv(FARG, sep="\t")
    farg_col = "gene" if "gene" in farg_table.columns else farg_table.columns[0]
    fargs = farg_table[farg_col].dropna().astype(str).str.upper().drop_duplicates().tolist()
    anchors_present = sorted(set(fargs).intersection(ANCHORS))
    all_tests = []
    for dataset in DATASETS:
        _, tests = run_dataset(dataset, fargs)
        all_tests.append(tests)
    tests = pd.concat(all_tests, ignore_index=True)
    tests.to_csv(OUT / "cross_dataset_paired_tests.tsv", sep="\t", index=False)
    meta = meta_analysis(tests)
    meta.to_csv(OUT / "cross_dataset_farg95_exclusion_meta.tsv", sep="\t", index=False)
    manifest = {
        "original_gene_count": len(fargs),
        "slamf8_in_original": "SLAMF8" in fargs,
        "prespecified_anchor_rule": sorted(ANCHORS),
        "anchor_genes_present_in_farg95": anchors_present,
        "score_definitions": {
            "FARG95_original": "all FARG95 genes present",
            "FARG95_SLAMF8_excluded": "FARG95 excluding SLAMF8",
            "FARG95_anchor_excluded": (
                "FARG95 after applying the prespecified SLAMF8/SPP1/GPNMB/TREM2/CD9 rule; "
                "only SLAMF8 is a member of FARG95"
            ),
        },
        "threshold": "within-sample SLAMF8 > 75th percentile; low <= 75th percentile",
        "seed": 20260920,
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    (OUT / "README.md").write_text(
        "# FARG95 self-containment repair\n\n"
        "This is a sensitivity re-analysis only. The original results are not overwritten. "
        "The primary contrast is unchanged; only FARG95 gene membership is altered. "
        "The prespecified anchor rule names SLAMF8, SPP1, GPNMB, TREM2 and CD9, but only SLAMF8 "
        "is a member of FARG95; the SLAMF8-excluded and anchor-excluded scores are therefore identical.\n",
        encoding="utf-8",
    )
    print(meta.to_string(index=False))


if __name__ == "__main__":
    main()
