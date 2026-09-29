#!/usr/bin/env python
import os
"""Validate SLAMF8/FARG95/HSC signals in newly downloaded public cohorts."""

from pathlib import Path
import gzip
import re
import requests

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
DATA = ROOT / "data_processed" / "bulk_geo" / "expanded_public"
META = ROOT / "results" / "bulk_validation"
OUT = ROOT / "results" / "bulk_validation" / "expanded_public"
FARG = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1])) / "data" / "reference" / "FARG95_gene_set_provenance.tsv"
SETS = {
    "FARG95": pd.read_csv(FARG, sep="\t")["gene"].astype(str).str.upper().dropna().unique().tolist(),
    "Iron_homeostasis": ["TFRC", "SLC11A2", "FTH1", "FTL", "SLC40A1", "NCOA4", "HMOX1", "STEAP3", "PCBP1", "PCBP2", "HAMP"],
    "HSC_activation": ["YAP1", "NUAK2", "COL1A1", "COL1A2", "COL3A1", "LOX", "ACTA2", "TAGLN", "TGFB1", "PDGFRA", "PDGFRB"],
}
TARGETS = ["SLAMF8", "SPP1", "GPNMB", "LGALS3", "PDGFB", "PDGFC", "CXCL16", "TGFB1"]


def gene_map() -> pd.Series:
    path = ROOT / "data_processed" / "bulk_geo" / "full" / "Homo_sapiens.gene_info.gz"
    if not path.exists():
        url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
        with requests.get(url, stream=True, timeout=120) as response:
            response.raise_for_status()
            with path.open("wb") as handle:
                for chunk in response.iter_content(1024 * 1024):
                    if chunk:
                        handle.write(chunk)
    table = pd.read_csv(path, sep="\t", compression="gzip", dtype=str, low_memory=False)
    table.columns = table.columns.str.removeprefix("#")
    return table.drop_duplicates("GeneID").set_index("GeneID")["Symbol"].str.upper()


def score(matrix: pd.DataFrame, genes: list[str]) -> pd.Series:
    matched = [gene for gene in genes if gene in matrix.index]
    if len(matched) < 3:
        return pd.Series(np.nan, index=matrix.columns)
    values = matrix.loc[matched]
    z = values.sub(values.mean(axis=1), axis=0).div(values.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def parse_metadata(accession: str) -> pd.DataFrame:
    frame = pd.read_csv(META / f"{accession}_sample_metadata.tsv", sep="\t")
    rows = []
    for _, item in frame.iterrows():
        row = {"sample": str(item["!Sample_title"]), "description": str(item.get("!Sample_description", "")), "gsm": item["!Sample_geo_accession"]}
        for field in [name for name in frame if "characteristics" in name]:
            value = str(item[field])
            if ":" in value:
                key, val = value.split(":", 1)
                row[key.strip().lower()] = val.strip()
        rows.append(row)
    return pd.DataFrame(rows)


def load_gse167523() -> tuple[pd.DataFrame, pd.DataFrame]:
    matrix = pd.read_csv(DATA / "GSE167523_Raw_gene_counts_matrix.txt.gz", sep="\t", index_col=0)
    matrix.index = matrix.index.astype(str).str.upper()
    matrix = matrix.groupby(level=0).sum()
    lib = matrix.sum(axis=0).replace(0, np.nan)
    matrix = np.log2(matrix.div(lib, axis=1).mul(1e6) + 1)
    metadata = parse_metadata("GSE167523")
    metadata["sample"] = [f"NAFLD{i:02d}" for i in range(1, len(metadata) + 1)]
    metadata["group"] = metadata["disease subtype"].str.upper()
    return matrix, metadata


def load_gse193066(filename: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(DATA / filename, sep="\t", compression="gzip", skiprows=2)
    raw = raw.rename(columns={raw.columns[0]: "gene_name", raw.columns[1]: "gene_id"}).set_index("gene_name")
    raw = raw.drop(columns=["gene_id"])
    raw.index = raw.index.astype(str).str.upper()
    raw = raw[~raw.index.duplicated(keep="first")]
    matrix = np.log2(raw.astype(float) + 1)
    metadata = parse_metadata("GSE193066")
    metadata = metadata[metadata["sample"].isin(matrix.columns)].copy()
    metadata["fibrosis_numeric"] = pd.to_numeric(metadata["fibrosis stage"], errors="coerce")
    metadata["biopsy"] = metadata["biopsy"].str.extract(r"(\d+)", expand=False)
    return matrix, metadata


def load_gse126848(ensembl: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(DATA / "GSE126848_Gene_counts_raw.txt.gz", sep="\t", index_col=0)
    ids = raw.index.astype(str).str.replace(r"\..*", "", regex=True)
    raw.index = ensembl.reindex(ids).to_numpy()
    raw = raw.loc[raw.index.notna()].groupby(level=0).sum()
    lib = raw.sum(axis=0).replace(0, np.nan)
    matrix = np.log2(raw.div(lib, axis=1).mul(1e6) + 1)
    metadata = parse_metadata("GSE126848")
    matrix_columns = set(matrix.columns.astype(str))
    description_ids = metadata["description"].str.extract(r"(\d+)$", expand=False)
    metadata["sample"] = np.where(description_ids.isin(matrix_columns), description_ids, metadata["sample"])
    metadata["sample"] = metadata["sample"].astype(str)
    metadata = metadata[metadata["sample"].isin(matrix.columns)].copy()
    metadata["group"] = metadata["disease"].str.upper()
    return matrix, metadata


def analyze(accession: str, matrix: pd.DataFrame, metadata: pd.DataFrame, mode: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = metadata.loc[metadata["sample"].isin(matrix.columns)].set_index("sample").copy()
    output = metadata.copy()
    for gene in TARGETS:
        output[gene] = matrix.reindex(index=[gene], columns=output.index).iloc[0] if gene in matrix.index else np.nan
    for name, genes in SETS.items():
        output[name] = score(matrix.loc[:, output.index], genes)
    tests = []
    if mode == "fibrosis":
        subset = output.dropna(subset=["fibrosis_numeric"])
        for feature in TARGETS + list(SETS):
            pair = subset[["fibrosis_numeric", feature]].dropna()
            if len(pair) >= 20 and pair.fibrosis_numeric.nunique() >= 2:
                rho, p = spearmanr(pair.fibrosis_numeric, pair[feature])
                tests.append({"accession": accession, "analysis": "fibrosis_stage", "feature": feature, "n": len(pair), "rho": rho, "p_value": p})
    else:
        for contrast, case, control in [("NASH_vs_NAFL", "NASH", "NAFL"), ("NASH_vs_nonNASH", "NASH", None)]:
            if contrast == "NASH_vs_NAFL":
                subset = output[output.group.isin([case, control])]
            else:
                subset = output.copy()
                subset["_control"] = ~subset.group.eq(case)
            case_values = subset[subset.group.eq(case)]
            control_values = subset[subset.group.eq(control)] if control else subset[~subset.group.eq(case)]
            for feature in TARGETS + list(SETS):
                a = case_values[feature].dropna().to_numpy()
                b = control_values[feature].dropna().to_numpy()
                if len(a) >= 5 and len(b) >= 5:
                    stat, p = mannwhitneyu(a, b, alternative="two-sided")
                    auc = (np.greater.outer(a, b).mean() + 0.5 * np.equal.outer(a, b).mean())
                    tests.append({"accession": accession, "analysis": contrast, "feature": feature, "n_case": len(a), "n_control": len(b), "median_difference": np.median(a) - np.median(b), "auc": auc, "p_value": p})
    return output.reset_index(), pd.DataFrame(tests)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gene_map()
    import anndata as ad
    reference = ad.read_h5ad(ROOT / "data_processed" / "GSE212837_myeloid_reclustered.h5ad", backed="r")
    ensembl = reference.var.set_index(reference.var["gene_id"].astype(str).str.replace(r"\..*", "", regex=True))["gene"].astype(str).str.upper()
    all_values, all_tests = [], []
    for accession, loader, mode in [("GSE167523", load_gse167523, "group"), ("GSE126848", lambda: load_gse126848(ensembl), "group")]:
        matrix, metadata = loader()
        values, tests = analyze(accession, matrix, metadata, mode)
        values.to_csv(OUT / f"{accession}_scores.tsv", sep="\t", index=False)
        all_values.append(values)
        all_tests.append(tests)
    for filename in ["GSE193066_NAFLD.HUn106.gct.gz", "GSE193066_NAFLD.HUn164.gct.gz"]:
        matrix, metadata = load_gse193066(filename)
        values, tests = analyze("GSE193066_" + ("HUn106" if "HUn106" in filename else "HUn164"), matrix, metadata, "fibrosis")
        values.to_csv(OUT / f"{accession if False else filename.replace('.gct.gz','')}_scores.tsv", sep="\t", index=False)
        all_values.append(values)
        all_tests.append(tests)
    tests = pd.concat(all_tests, ignore_index=True)
    tests["fdr_within_analysis"] = tests.groupby(["accession", "analysis"])["p_value"].transform(lambda x: multipletests(x, method="fdr_bh")[1])
    tests.to_csv(OUT / "expanded_public_validation.tsv", sep="\t", index=False)
    pd.DataFrame({"accession": ["GSE167523", "GSE193066_HUn106", "GSE193066_HUn164", "GSE126848"], "matrix_type": ["raw counts", "RLE-normalized GCT", "RLE-normalized GCT", "raw counts"], "clinical_endpoint": ["NAFL vs NASH", "fibrosis stage", "fibrosis stage/paired biopsies", "healthy/obese/NAFLD/NASH"]}).to_csv(OUT / "expanded_public_cohort_manifest.tsv", sep="\t", index=False)
    print(tests[tests.feature.isin(["SLAMF8", "SPP1", "FARG95", "Iron_homeostasis", "HSC_activation"])].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
