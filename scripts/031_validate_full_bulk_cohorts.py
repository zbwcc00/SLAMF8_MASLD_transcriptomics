#!/usr/bin/env python
import os
"""Independent complete-transcriptome MASLD fibrosis validation."""

from pathlib import Path
import gzip
import re
import requests
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm

ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
FULL = ROOT / "data_processed" / "bulk_geo" / "full"
OUT = ROOT / "results" / "bulk_validation"
FARG = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1])) / "data" / "reference" / "FARG95_gene_set_provenance.tsv"
SETS = {
    "FARG95": pd.read_csv(FARG, sep="\t")["gene"].str.upper().dropna().unique().tolist(),
    "Iron_homeostasis": ["TFRC", "SLC11A2", "FTH1", "FTL", "SLC40A1", "NCOA4", "HMOX1", "STEAP3", "PCBP1", "PCBP2", "HAMP"],
    "Senescence_SASP": ["CDKN1A", "CDKN2A", "GLB1", "SERPINE1", "GDF15", "IL6", "IL1B", "CXCL8", "MMP3", "MMP9", "CCL2", "TGFB1", "IGFBP7"],
    "HSC_activation": ["YAP1", "NUAK2", "COL1A1", "COL1A2", "COL3A1", "LOX", "ACTA2", "TAGLN", "TGFB1", "PDGFRA", "PDGFRB"],
}
TARGETS = ["SLAMF8", "GPNMB", "LGALS3", "SPP1", "PDGFB", "CXCL16", "TNF"]


def gene_maps():
    sc = ROOT / "data_processed" / "GSE212837_myeloid_reclustered.h5ad"
    import anndata as ad
    obj = ad.read_h5ad(sc, backed="r")
    genes = obj.var[["gene_id", "gene"]].copy()
    genes["gene_id"] = genes.gene_id.astype(str).str.replace(r"\..*", "", regex=True)
    genes["gene"] = genes.gene.astype(str).str.upper()
    gene_map = genes.drop_duplicates("gene_id").set_index("gene_id")["gene"]
    path = FULL / "Homo_sapiens.gene_info.gz"
    if not path.exists():
        url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
        with requests.get(url, stream=True, timeout=120) as response:
            response.raise_for_status()
            with path.open("wb") as handle:
                for chunk in response.iter_content(1024 * 1024):
                    if chunk:
                        handle.write(chunk)
    ncbi = pd.read_csv(path, sep="\t", compression="gzip", dtype=str, low_memory=False)
    ncbi.columns = ncbi.columns.str.removeprefix("#")
    entrez = ncbi.drop_duplicates("GeneID").set_index("GeneID")["Symbol"].str.upper()
    return gene_map, entrez


def load_metadata(accession):
    frame = pd.read_csv(OUT / f"{accession}_sample_metadata.tsv", sep="\t")
    columns = [name for name in frame if "characteristics" in name]
    rows = []
    for _, item in frame.iterrows():
        row = {"gsm": item["!Sample_geo_accession"], "title": item["!Sample_title"]}
        for field in columns:
            value = str(item[field])
            if ":" in value:
                key, val = value.split(":", 1)
                row[key.strip().lower()] = val.strip()
        rows.append(row)
    return pd.DataFrame(rows)


def load_expression(accession, gene_map, entrez):
    if accession == "GSE162694":
        frame = pd.read_csv(FULL / "GSE162694_raw_counts.csv.gz", index_col=0)
        ids = frame.index.astype(str).str.split(".").str[0]
        frame.index = gene_map.reindex(ids).to_numpy()
        frame = frame.loc[frame.index.notna()]
        frame = frame.groupby(level=0).sum()
        lib = frame.sum(axis=0)
        return np.log2(frame.div(lib, axis=1).mul(1e6) + 1)
    frame = pd.read_csv(FULL / "GSE130970_all_sample_salmon_tximport_TPM_entrez_gene_ID.csv.gz", index_col=0)
    frame.index = entrez.reindex(frame.index.astype(str)).to_numpy()
    frame = frame.loc[frame.index.notna()].groupby(level=0).sum()
    return np.log2(frame + 1)


def score(frame, genes):
    matched = [gene for gene in genes if gene in frame.index]
    if len(matched) < 3:
        return pd.Series(np.nan, index=frame.columns), len(matched)
    values = frame.loc[matched]
    gene_sd = values.std(axis=1).replace(0, np.nan)
    zscores = values.sub(values.mean(axis=1), axis=0).div(gene_sd, axis=0)
    return zscores.mean(axis=0), len(matched)


def run(accession, gene_map, entrez):
    matrix = load_expression(accession, gene_map, entrez)
    metadata = load_metadata(accession)
    if accession == "GSE162694":
        metadata["sample"] = metadata.title.str.extract(r"(548nash\d+)", flags=re.I, expand=False)
        metadata["fibrosis"] = pd.to_numeric(metadata["fibrosis stage"], errors="coerce")
        metadata["normal_histology"] = metadata["fibrosis stage"].eq("normal liver histology")
        metadata.loc[metadata.normal_histology, "fibrosis"] = 0
    else:
        metadata["sample"] = metadata.title
        metadata["fibrosis"] = pd.to_numeric(metadata["fibrosis stage"], errors="coerce")
        metadata["normal_histology"] = False
    metadata["age_numeric"] = pd.to_numeric(metadata.get("age", metadata.get("age at biopsy")), errors="coerce")
    metadata["sex_female"] = metadata["sex"].astype(str).str.strip().str.lower().isin(["female", "f"]).astype(int)
    metadata["nas_numeric"] = pd.to_numeric(metadata.get("nas score", metadata.get("nafld activity score")), errors="coerce")
    joined = metadata.set_index("sample").join(matrix.T, how="inner")
    if joined.empty:
        raise ValueError(f"No sample IDs matched {accession}")
    if joined.index.duplicated().any():
        raise ValueError(f"Duplicate specimen ID in {accession}")
    output = joined[["gsm", "title", "fibrosis", "normal_histology", "age_numeric", "sex_female", "nas_numeric"]].copy()
    coverage = []
    for gene in TARGETS:
        output[gene] = joined[gene] if gene in matrix.index else np.nan
        coverage.append({"dataset": accession, "feature": gene, "genes_present": int(gene in matrix.index), "genes_requested": 1})
    for name, genes in SETS.items():
        values, n = score(matrix.loc[:, output.index], genes)
        output[name] = values
        coverage.append({"dataset": accession, "feature": name, "genes_present": n, "genes_requested": len(genes)})
    output.to_csv(OUT / f"{accession}_full_bulk_scores.tsv", sep="\t", index_label="sample")
    tests = []
    for scope in (["all_biopsies", "NAFLD_only"] if accession == "GSE162694" else ["NAFLD_only"]):
        scoped = output if scope == "all_biopsies" else output[~output.normal_histology]
        for feature in TARGETS + list(SETS):
            subset = scoped[["fibrosis", feature]].dropna()
            if len(subset) < 15 or subset.fibrosis.nunique() < 2:
                continue
            rho, p = spearmanr(subset.fibrosis, subset[feature])
            cases = subset.loc[subset.fibrosis >= 2, feature]
            controls = subset.loc[subset.fibrosis < 2, feature]
            auc = float(np.mean(cases.values[:, None] > controls.values[None, :]) + 0.5 * np.mean(cases.values[:, None] == controls.values[None, :])) if len(cases) and len(controls) else np.nan
            tests.append({"dataset": accession, "scope": scope, "feature": feature, "n": len(subset), "fibrosis_ge2": len(cases), "rho": rho, "p_value": p, "auc_f_ge2": auc})
    adjusted = []
    for scope in (["all_biopsies", "NAFLD_only"] if accession == "GSE162694" else ["NAFLD_only"]):
        scoped = output if scope == "all_biopsies" else output[~output.normal_histology]
        for feature in ["SLAMF8", "FARG95", "Iron_homeostasis", "HSC_activation"]:
            for covariates in [["age_numeric", "sex_female"], ["age_numeric", "sex_female", "nas_numeric"]]:
                fields = ["fibrosis", feature] + covariates
                subset = scoped[fields].dropna()
                if len(subset) < 25 or subset.fibrosis.nunique() < 2:
                    continue
                std = subset[feature].std()
                if not np.isfinite(std) or std == 0:
                    continue
                design = subset[covariates].copy()
                design["feature_per_sd"] = (subset[feature] - subset[feature].mean()) / std
                model = sm.OLS(subset.fibrosis, sm.add_constant(design)).fit(cov_type="HC3")
                adjusted.append({"dataset": accession, "scope": scope, "feature": feature,
                                 "adjustment": "age_sex_NAS" if "nas_numeric" in covariates else "age_sex",
                                 "n": len(subset), "beta_fibrosis_stage_per_sd": model.params["feature_per_sd"],
                                 "p_value": model.pvalues["feature_per_sd"]})
    links = []
    cases = output.loc[~output.normal_histology]
    for feature in ["GPNMB", "FARG95", "Iron_homeostasis", "Senescence_SASP", "HSC_activation"]:
        subset = cases[["SLAMF8", feature]].dropna()
        rho, p = spearmanr(subset["SLAMF8"], subset[feature])
        links.append({"dataset": accession, "feature": feature, "n": len(subset),
                      "rho_with_SLAMF8": rho, "p_value": p})
    return pd.DataFrame(tests), pd.DataFrame(coverage), pd.DataFrame(adjusted), pd.DataFrame(links), len(output)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    gene_map, entrez = gene_maps()
    tested, covered, adjusted, linked = [], [], [], []
    for accession in ["GSE162694", "GSE130970"]:
        table, coverage, sensitivity, links, size = run(accession, gene_map, entrez)
        tested.append(table)
        covered.append(coverage)
        adjusted.append(sensitivity)
        linked.append(links)
        print(accession, size, "samples", flush=True)
    tests = pd.concat(tested, ignore_index=True)
    tests["fdr_within_dataset_scope"] = tests.groupby(["dataset", "scope"])["p_value"].transform(lambda values: multipletests(values, method="fdr_bh")[1])
    tests.to_csv(OUT / "full_bulk_fibrosis_validation.tsv", sep="\t", index=False)
    pd.concat(covered, ignore_index=True).to_csv(OUT / "full_bulk_gene_coverage.tsv", sep="\t", index=False)
    pd.concat(adjusted, ignore_index=True).to_csv(OUT / "full_bulk_covariate_sensitivity.tsv", sep="\t", index=False)
    links = pd.concat(linked, ignore_index=True)
    links["fdr_within_dataset"] = links.groupby("dataset")["p_value"].transform(lambda values: multipletests(values, method="fdr_bh")[1])
    links.to_csv(OUT / "full_bulk_SLAMF8_module_correlations.tsv", sep="\t", index=False)
    print(tests[tests.feature.isin(["SLAMF8", "FARG95", "Iron_homeostasis", "HSC_activation"])].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
