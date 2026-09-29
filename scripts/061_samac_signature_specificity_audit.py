#!/usr/bin/env python
"""Audit whether the bulk SAMac-like/fibrosis association is gene-set specific."""
from __future__ import annotations
import os

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
SCORES = ROOT / "results" / "revision_p2_5" / "bulk_immune_and_models" / "bulk_cell_state_signature_scores.tsv"
OUT = ROOT / "results" / "revision_p0" / "signature_specificity"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260923
N_RANDOM = 10000

SAMAC = ["TREM2", "GPNMB", "CD9", "SPP1", "LGALS3", "LPL", "APOC1", "FABP5", "CTSB", "CTSD"]
SAMAC_NO_SPP1 = [g for g in SAMAC if g != "SPP1"]
HSC_ECM_EXTERNAL = ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "DCN", "LUM", "ACTA2", "TAGLN", "PDGFRB"]


def score_genes(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    found = [g for g in genes if g in expr.index]
    if len(found) < 3:
        return pd.Series(np.nan, index=expr.columns), found
    x = expr.loc[found]
    z = x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0), found


def hsc_sets() -> dict[str, list[str]]:
    p = ROOT / "results" / "revision_p2_5" / "bulk_immune_and_models" / "bulk_cell_state_signature_coverage.tsv"
    cov = pd.read_csv(p, sep="\t")
    rows = cov[(cov.dataset == "GSE162694") & (cov.signature == "Activated_HSC_fibrogenic")]
    if rows.empty:
        return {"Bulk_HSC_fibrogenic_signature": HSC_ECM_EXTERNAL}
    found = str(rows.iloc[0].genes_found).split(";")
    return {"Bulk_HSC_fibrogenic_signature": [g for g in HSC_ECM_EXTERNAL if g in found]}


def random_matched_sets(expr: pd.DataFrame, signature: list[str], rng: np.random.Generator) -> tuple[list[list[str]], int]:
    available = list(expr.index)
    signature = [g for g in signature if g in expr.index]
    if len(signature) < 3:
        return [], 0
    means = expr.mean(axis=1)
    bins = pd.qcut(means.rank(method="first"), q=min(10, max(2, len(expr) // 20)), labels=False, duplicates="drop")
    gene_bin = pd.Series(bins, index=expr.index)
    groups: dict[object, list[str]] = {}
    for gene in available:
        groups.setdefault(gene_bin[gene], []).append(gene)
    sets = []
    for _ in range(N_RANDOM):
        draw = []
        for gene in signature:
            pool = [g for g in groups[gene_bin[gene]] if g not in signature and g not in draw]
            if not pool:
                pool = [g for g in available if g not in signature and g not in draw]
            if not pool:
                break
            draw.append(str(rng.choice(pool)))
        if len(draw) == len(signature):
            sets.append(draw)
    return sets, len(signature)


def main() -> None:
    scores = pd.read_csv(SCORES, sep="\t")
    rng = np.random.default_rng(SEED)
    association_rows, loo_rows, random_rows, coverage_rows = [], [], [], []
    for dataset, samples in scores.groupby("dataset", sort=True):
        if "normal_histology" in samples:
            samples = samples.loc[~samples["normal_histology"].astype(bool)].copy()
        if dataset == "GSE162694":
            raw = pd.read_csv(ROOT / "data_processed" / "bulk_geo" / "full" / "GSE162694_raw_counts.csv.gz", index_col=0)
            import anndata as ad
            gene_table = ad.read_h5ad(ROOT / "data_processed" / "GSE212837_myeloid_reclustered.h5ad", backed="r").var
            mapper = gene_table.drop_duplicates("gene_id").set_index("gene_id")["gene"].astype(str).str.upper()
            raw.index = mapper.reindex(raw.index.astype(str).str.split(".").str[0]).to_numpy()
            raw = raw.loc[raw.index.notna()].groupby(level=0).sum()
            expr = np.log2(raw.div(raw.sum(axis=0), axis=1).mul(1e6) + 1)
            expr.columns = expr.columns.astype(str)
        elif dataset == "GSE130970":
            raw = pd.read_csv(ROOT / "data_processed" / "bulk_geo" / "full" / "GSE130970_all_sample_salmon_tximport_TPM_entrez_gene_ID.csv.gz", index_col=0)
            gene_map = pd.read_csv(ROOT / "data_processed" / "bulk_geo" / "full" / "Homo_sapiens.gene_info.gz", sep="\t", compression="gzip", dtype=str)
            gene_map.columns = gene_map.columns.str.removeprefix("#")
            mapper = gene_map.drop_duplicates("GeneID").set_index("GeneID")["Symbol"].str.upper()
            raw.index = mapper.reindex(raw.index.astype(str)).to_numpy()
            raw = raw.loc[raw.index.notna()].groupby(level=0).sum()
            expr = np.log2(raw + 1)
            expr.columns = expr.columns.astype(str)
        else:
            continue
        common = [s for s in samples["sample"].astype(str) if s in expr.columns]
        samples = samples.set_index(samples["sample"].astype(str)).loc[common]
        expr = expr.loc[:, common]
        sample_ids = samples.index
        print(dataset, "matched_biopsies", len(sample_ids), "expression_samples", expr.shape[1], flush=True)

        signatures = {"SAMac_like": SAMAC, "SAMac_like_no_SPP1": SAMAC_NO_SPP1, **hsc_sets()}
        for name, members in signatures.items():
            val, found = score_genes(expr, members)
            q = pd.DataFrame({"fibrosis": samples.loc[sample_ids, "fibrosis"].to_numpy(), "score": val.reindex(sample_ids).to_numpy()}).dropna()
            rho, pval = spearmanr(q.fibrosis, q.score)
            association_rows.append({"dataset": dataset, "signature": name, "n": len(q), "n_genes_present": len(found),
                                     "genes_present": ";".join(found), "rho": rho, "p_value": pval})
            coverage_rows.append({"dataset": dataset, "signature": name, "n_requested": len(members),
                                  "n_found": len(found), "genes_found": ";".join(found)})
            for gene in found:
                reduced = [g for g in found if g != gene]
                s_loo, _ = score_genes(expr, reduced)
                qloo = pd.DataFrame({"fibrosis": samples.loc[sample_ids, "fibrosis"].to_numpy(), "score": s_loo.reindex(sample_ids).to_numpy()}).dropna()
                lrho, lp = spearmanr(qloo.fibrosis, qloo.score)
                loo_rows.append({"dataset": dataset, "signature": name, "gene_removed": gene, "n": len(qloo), "rho": lrho, "p_value": lp})

        present_samac = [g for g in SAMAC if g in expr.index]
        matched, size = random_matched_sets(expr, present_samac, rng)
        null_rhos = []
        for genes in matched:
            s, _ = score_genes(expr, genes)
            q = pd.DataFrame({"fibrosis": samples.loc[sample_ids, "fibrosis"].to_numpy(), "score": s.reindex(sample_ids).to_numpy()}).dropna()
            null_rhos.append(float(spearmanr(q.fibrosis, q.score).statistic))
        samac, found = score_genes(expr, present_samac)
        qs = pd.DataFrame({"fibrosis": samples.loc[sample_ids, "fibrosis"].to_numpy(), "score": samac.reindex(sample_ids).to_numpy()}).dropna()
        observed = float(spearmanr(qs.fibrosis, qs.score).statistic)
        exceed = int(np.sum(np.asarray(null_rhos) >= observed)) if null_rhos else 0
        random_rows.append({"dataset": dataset, "observed_rho": observed, "n_signature_genes": size,
                            "n_random_sets": len(null_rhos), "random_rho_median": float(np.median(null_rhos)) if null_rhos else np.nan,
                            "random_rho_q95": float(np.quantile(null_rhos, .95)) if null_rhos else np.nan,
                            "random_rho_q99": float(np.quantile(null_rhos, .99)) if null_rhos else np.nan,
                            "empirical_one_sided_p": (exceed + 1) / (len(null_rhos) + 1) if null_rhos else np.nan, "seed": SEED,
                            "matching": "gene-wise mean-expression decile; no signature genes"})

    assoc = pd.DataFrame(association_rows)
    assoc["fdr_within_dataset"] = assoc.groupby("dataset").p_value.transform(lambda x: multipletests(x, method="fdr_bh")[1])
    loo = pd.DataFrame(loo_rows)
    loo["fdr_within_dataset_signature"] = loo.groupby(["dataset", "signature"]).p_value.transform(lambda x: multipletests(x, method="fdr_bh")[1])
    assoc.to_csv(OUT / "signature_associations.tsv", sep="\t", index=False)
    loo.to_csv(OUT / "leave_one_gene_out.tsv", sep="\t", index=False)
    pd.DataFrame(random_rows).to_csv(OUT / "matched_random_signature_null.tsv", sep="\t", index=False)
    pd.DataFrame(coverage_rows).to_csv(OUT / "gene_coverage.tsv", sep="\t", index=False)
    (OUT / "README.md").write_text(
        "# SAMac-like signature specificity audit\n\n"
        "The audit was limited to the two primary complete NAFLD bulk cohorts. It reports leave-one-gene-out associations for the pre-specified SAMac-like and externally defined activated-HSC/ECM panels, plus 10,000 random gene sets matched by per-gene mean-expression decile for the SAMac-like score. The matched random-set test is a specificity sensitivity analysis, not an independent validation cohort.\n\n"
        f"Random seed: {SEED}. Sample unit: biopsy.\n", encoding="utf-8")
    print(assoc.to_string(index=False))
    print(pd.DataFrame(random_rows).to_string(index=False))


if __name__ == "__main__":
    main()
