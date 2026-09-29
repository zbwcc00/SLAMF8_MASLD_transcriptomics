#!/usr/bin/env python
"""External bulk validation of SLAMF8 and ferro-aging-like programs."""

from __future__ import annotations
import os

import csv
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests


PROJECT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
FARG = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1])) / "data" / "reference" / "FARG95_gene_set_provenance.tsv"
IN = PROJECT / "data_processed" / "bulk_geo"
PLAT = IN / "platforms"
META = PROJECT / "results" / "bulk_validation"
OUT = PROJECT / "results" / "bulk_validation"

GENE_SETS = {
    "FARG95": pd.read_csv(FARG, sep="\t")["gene"].astype(str).str.upper().drop_duplicates().tolist(),
    "Iron_homeostasis": ["TFRC", "SLC11A2", "FTH1", "FTL", "SLC40A1", "NCOA4", "HMOX1", "STEAP3", "PCBP1", "PCBP2", "HAMP"],
    "PUFA_ACSL4": ["ACSL4", "LPCAT3", "AGPAT3", "ELOVL5", "FADS1", "FADS2", "ALOX15", "PLA2G6"],
    "Lipid_peroxidation": ["ACSL4", "ALOX5", "ALOX12", "ALOX15", "POR", "NOX1", "NOX2", "CYBB", "PTGS2"],
    "Antioxidant_defense": ["GPX4", "SLC7A11", "GCLC", "GCLM", "NFE2L2", "FSP1", "AIFM2", "DHODH", "GCH1"],
    "Senescence_SASP": ["CDKN1A", "CDKN2A", "GLB1", "SERPINE1", "GDF15", "IL6", "IL1B", "CXCL8", "MMP3", "MMP9", "CCL2", "TGFB1", "IGFBP7"],
    "Inflammation": ["IL1B", "TNF", "NFKBIA", "CXCL8", "CCL2", "CCL3", "CCL4", "NLRP3", "STAT1", "IRF1"],
    "HSC_activation": ["YAP1", "NUAK2", "COL1A1", "COL1A2", "COL3A1", "LOX", "ACTA2", "TAGLN", "TGFB1", "PDGFRA", "PDGFRB"],
}


PLATFORMS = {"GSE48452": "GPL11532", "GSE63067": "GPL570"}


def probe_map(platform: str) -> pd.DataFrame:
    path = PLAT / f"{platform}.annot.gz"
    marker_line = None
    with gzip.open(path, "rt", errors="replace") as handle:
        for line_number, line in enumerate(handle):
            if line.startswith("!platform_table_begin"):
                marker_line = line_number
                break
    table = pd.read_csv(path, sep="\t", compression="gzip", skiprows=marker_line + 1, dtype=str)
    mapping = table[["ID", "Gene symbol"]].dropna()
    mapping["Gene symbol"] = mapping["Gene symbol"].str.split("///")
    mapping = mapping.explode("Gene symbol")
    mapping["ID"] = mapping["ID"].astype(str).str.upper()
    mapping["Gene symbol"] = mapping["Gene symbol"].astype(str).str.strip().str.upper()
    mapping = mapping[mapping["Gene symbol"].ne("")]
    return mapping.drop_duplicates()


def read_matrix(accession: str) -> pd.DataFrame:
    path = IN / f"{accession}_series_matrix.txt.gz"
    marker_line = None
    with gzip.open(path, "rt", errors="replace") as handle:
        for line_number, line in enumerate(handle):
            if line.startswith("!series_matrix_table_begin"):
                marker_line = line_number
                break
    if marker_line is None:
        raise ValueError(f"No matrix table in {accession}")
    matrix = pd.read_csv(path, sep="\t", compression="gzip", skiprows=marker_line + 1, index_col=0)
    matrix.index = matrix.index.astype(str).str.strip('"').str.upper()
    matrix.columns = matrix.columns.astype(str).str.strip('"')
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    matrix = matrix.loc[~matrix.index.str.startswith("!")]
    matrix = matrix.groupby(level=0).mean()
    if accession in PLATFORMS:
        mapping = probe_map(PLATFORMS[accession])
        mapping = mapping[mapping["ID"].isin(matrix.index)]
        mapped = matrix.loc[mapping["ID"]].copy()
        mapped.index = mapping["Gene symbol"].to_numpy()
        matrix = mapped.groupby(level=0).mean()
    return matrix


def parse_characteristics(row: pd.Series) -> dict[str, str]:
    values = {}
    for value in row.dropna().astype(str):
        match = re.match(r"\s*([^:]+):\s*(.*)\s*$", value)
        if match:
            values[match.group(1).strip().lower()] = match.group(2).strip()
    return values


def metadata(accession: str) -> pd.DataFrame:
    table = pd.read_csv(META / f"{accession}_sample_metadata.tsv", sep="\t")
    characteristic_columns = [x for x in table.columns if "characteristics" in x]
    records = []
    for _, row in table.iterrows():
        record = {"sample": str(row["!Sample_geo_accession"]).strip()}
        for column in characteristic_columns:
            record.update(parse_characteristics(row[[column]]))
        records.append(record)
    return pd.DataFrame(records).set_index("sample")


def normalize(matrix: pd.DataFrame) -> pd.DataFrame:
    values = matrix.copy()
    finite_max = np.nanmax(values.to_numpy())
    if finite_max > 50:
        values = np.log2(values.clip(lower=0) + 1.0)
    return values


def score(matrix: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, int]:
    present = [gene for gene in genes if gene in matrix.index]
    if not present:
        return pd.Series(np.nan, index=matrix.columns), 0
    values = matrix.loc[present]
    z = values.sub(values.mean(axis=1), axis=0).div(values.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0), len(present)


def adjust(values: list[float]) -> list[float]:
    if not values:
        return []
    return multipletests(np.nan_to_num(values, nan=1.0), method="fdr_bh")[1].tolist()


def analyze(accession: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = normalize(read_matrix(accession))
    meta = metadata(accession)
    common = [sample for sample in expr.columns if sample in meta.index]
    expr = expr.loc[:, common]
    meta = meta.loc[common].copy()
    result = meta.copy()
    result["SLAMF8"] = expr.loc["SLAMF8"] if "SLAMF8" in expr.index else np.nan
    coverage = []
    for name, genes in GENE_SETS.items():
        result[name], count = score(expr, genes)
        coverage.append({"dataset": accession, "program": name, "genes_present": count, "genes_requested": len(genes)})

    tests = []
    stage_col = None
    if accession == "GSE163211":
        stage_col = "nafld stage"
        order = ["Normal", "Steatosis", "NASH_F0", "NASH_F1_F4"]
        result["stage_numeric"] = result[stage_col].map({x: i for i, x in enumerate(order)})
    elif accession == "GSE150734":
        stage_col = "fibrosis stage"
        result["stage_numeric"] = pd.to_numeric(result[stage_col], errors="coerce")
    elif accession == "GSE48452":
        stage_col = "fibrosis"
        result["stage_numeric"] = pd.to_numeric(result[stage_col], errors="coerce")
    elif accession == "GSE63067":
        stage_col = "disease status"
        result["stage_numeric"] = result[stage_col].map({"healthy": 0, "steatosis": 1, "non-alcoholic steatohepatitis": 2})

    for target in ["SLAMF8", "FARG95", "Iron_homeostasis", "PUFA_ACSL4", "Lipid_peroxidation", "Senescence_SASP", "Inflammation", "HSC_activation"]:
        if target not in result or stage_col is None:
            continue
        valid = result[[target, "stage_numeric"]].dropna()
        if len(valid) >= 4 and valid["stage_numeric"].nunique() >= 2:
            rho, p = spearmanr(valid["stage_numeric"], valid[target])
            tests.append({"dataset": accession, "target": target, "test": "Spearman_stage", "n": len(valid), "effect": rho, "p_value": p})

    if accession == "GSE163211":
        for target in ["SLAMF8", "FARG95", "Iron_homeostasis", "HSC_activation"]:
            groups = [result.loc[result[stage_col] == group, target].dropna().to_numpy() for group in ["Normal", "Steatosis", "NASH_F0", "NASH_F1_F4"]]
            groups = [x for x in groups if len(x)]
            if len(groups) >= 2:
                stat, p = kruskal(*groups)
                tests.append({"dataset": accession, "target": target, "test": "Kruskal_stage", "n": int(sum(map(len, groups))), "effect": stat, "p_value": p})
        for target in ["SLAMF8", "FARG95", "Iron_homeostasis", "HSC_activation"]:
            a = result.loc[result[stage_col] == "Normal", target].dropna()
            b = result.loc[result[stage_col] == "NASH_F1_F4", target].dropna()
            if len(a) and len(b):
                stat, p = mannwhitneyu(a, b, alternative="two-sided")
                effect = (b.median() - a.median()) / np.sqrt((a.var() + b.var()) / 2) if (a.var() + b.var()) > 0 else np.nan
                tests.append({"dataset": accession, "target": target, "test": "NASH_F1_F4_vs_Normal", "n": len(a) + len(b), "effect": effect, "p_value": p})

    result.to_csv(OUT / f"{accession}_bulk_scores.tsv", sep="\t")
    return pd.DataFrame(tests), pd.DataFrame(coverage)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    all_tests, all_coverage = [], []
    for accession in ["GSE163211", "GSE150734", "GSE48452", "GSE63067"]:
        print(f"analyzing {accession}", flush=True)
        tests, coverage = analyze(accession)
        all_tests.append(tests)
        all_coverage.append(coverage)
    tests = pd.concat(all_tests, ignore_index=True)
    tests["fdr"] = np.nan
    for test_name, idx in tests.groupby("test").groups.items():
        tests.loc[idx, "fdr"] = adjust(tests.loc[idx, "p_value"].tolist())
    tests.to_csv(OUT / "bulk_validation_tests.tsv", sep="\t", index=False)
    pd.concat(all_coverage, ignore_index=True).to_csv(OUT / "bulk_validation_gene_coverage.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
