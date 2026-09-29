#!/usr/bin/env python
"""Validate the SLAMF8/FARG95/HSC axis in additional human public cohorts."""

from __future__ import annotations
import os

import csv
import gzip
import re
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
DATA = ROOT / "data_processed" / "bulk_geo" / "expanded_public"
META = ROOT / "results" / "bulk_validation"
OUT = ROOT / "results" / "bulk_validation" / "expanded_public"
FARG = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1])) / "data" / "reference" / "FARG95_gene_set_provenance.tsv"

TARGETS = ["SLAMF8", "SPP1", "GPNMB", "LGALS3", "PDGFB", "PDGFC", "CXCL16", "TGFB1"]
SETS = {
    "FARG95": pd.read_csv(FARG, sep="\t")["gene"].astype(str).str.upper().dropna().unique().tolist(),
    "Iron_homeostasis": ["TFRC", "SLC11A2", "FTH1", "FTL", "SLC40A1", "NCOA4", "HMOX1", "STEAP3", "PCBP1", "PCBP2", "HAMP"],
    "HSC_activation": ["YAP1", "NUAK2", "COL1A1", "COL1A2", "COL3A1", "LOX", "ACTA2", "TAGLN", "TGFB1", "PDGFRA", "PDGFRB"],
}


def score(matrix: pd.DataFrame, genes: list[str]) -> pd.Series:
    matched = [gene for gene in genes if gene in matrix.index]
    if len(matched) < 3:
        return pd.Series(np.nan, index=matrix.columns)
    values = matrix.loc[matched].astype(float)
    z = values.sub(values.mean(axis=1), axis=0).div(values.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def parse_geo_matrix_metadata(path: Path) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    rows: dict[str, list[str]] = {}
    characteristic_rows: list[list[str]] = []
    with gzip.open(path, "rt", errors="replace", newline="") as handle:
        for line in handle:
            if not line.startswith("!"):
                if line.startswith("\"ID_REF\""):
                    break
                continue
            fields = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            if fields[0].startswith("!Sample_characteristics"):
                characteristic_rows.append(fields[1:])
            else:
                rows[fields[0]] = fields[1:]
    accessions = rows.get("!Sample_geo_accession", [])
    metadata = pd.DataFrame({"gsm": accessions})
    for key, values in rows.items():
        if key.startswith("!Sample_"):
            metadata[key] = values
    for values in characteristic_rows:
        for idx, value in enumerate(values):
            if idx >= len(metadata):
                continue
            if ":" in value:
                key, val = value.split(":", 1)
                metadata.loc[idx, key.strip().lower().replace(" ", "_")] = val.strip()
    return metadata, rows


def parse_series_table(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", errors="replace") as handle:
        table = pd.read_csv(handle, sep="\t", comment="!", index_col=0)
    table.index = table.index.astype(str).str.strip('"')
    table.columns = table.columns.astype(str).str.strip('"')
    return table.apply(pd.to_numeric, errors="coerce")


def probe_to_symbol() -> dict[str, str]:
    sqlite_path = DATA / "hugene20sttranscriptcluster.db" / "inst" / "extdata" / "hugene20sttranscriptcluster.sqlite"
    gene_info = ROOT / "data_processed" / "bulk_geo" / "full" / "Homo_sapiens.gene_info.gz"
    symbols = pd.read_csv(gene_info, sep="\t", compression="gzip", dtype=str, low_memory=False)
    symbols.columns = symbols.columns.str.removeprefix("#")
    symbols = symbols.drop_duplicates("GeneID").set_index("GeneID")["Symbol"].str.upper()
    conn = sqlite3.connect(sqlite_path)
    mapping = pd.read_sql_query("select probe_id, gene_id from probes where gene_id is not null", conn)
    conn.close()
    mapping["symbol"] = mapping["gene_id"].map(symbols)
    mapping = mapping.dropna(subset=["symbol"]).drop_duplicates("probe_id")
    return mapping.set_index("probe_id")["symbol"].to_dict()


def collapse_symbols(probe_matrix: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    symbols = pd.Series(probe_matrix.index, index=probe_matrix.index).map(mapping)
    keep = symbols.notna()
    collapsed = probe_matrix.loc[keep].copy()
    collapsed.index = symbols[keep].values
    return collapsed.groupby(level=0).mean()


def add_scores(matrix: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    metadata = metadata.copy()
    metadata["sample"] = metadata["sample"].astype(str)
    metadata = metadata[metadata["sample"].isin(matrix.columns)].set_index("sample")
    output = metadata.copy()
    for gene in TARGETS:
        output[gene] = matrix.reindex(index=[gene], columns=output.index).iloc[0] if gene in matrix.index else np.nan
    for name, genes in SETS.items():
        output[name] = score(matrix.loc[:, output.index], genes)
    return output


def test_groups(accession: str, output: pd.DataFrame, group_col: str, case: str, control: str | None) -> pd.DataFrame:
    rows = []
    case_values = output[output[group_col].astype(str).str.upper().eq(case.upper())]
    control_values = output[output[group_col].astype(str).str.upper().eq(control.upper())] if control else output[~output[group_col].astype(str).str.upper().eq(case.upper())]
    for feature in TARGETS + list(SETS):
        a = case_values[feature].dropna().to_numpy()
        b = control_values[feature].dropna().to_numpy()
        if len(a) >= 5 and len(b) >= 5:
            _, p = mannwhitneyu(a, b, alternative="two-sided")
            auc = (np.greater.outer(a, b).mean() + 0.5 * np.equal.outer(a, b).mean())
            rows.append({"accession": accession, "analysis": f"{case}_vs_{control or 'non' + case}", "feature": feature, "n_case": len(a), "n_control": len(b), "median_difference": np.median(a) - np.median(b), "auc": auc, "p_value": p})
    return pd.DataFrame(rows)


def load_gse281797() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(DATA / "GSE281797_tpm.txt.gz", sep="\t", compression="gzip")
    symbols = raw["HGNC symbol"].astype(str).str.upper()
    matrix = raw.drop(columns=["Gene stable ID", "HGNC symbol", "Gene description", "NCBI gene (formerly Entrezgene) ID"]).copy()
    matrix.index = symbols
    matrix = matrix[~matrix.index.isin(["NAN", "NONE"])].groupby(level=0).mean()
    soft = DATA / "GSE281797_sample_metadata.tsv"
    metadata = pd.read_csv(soft, sep="\t")
    metadata = metadata.rename(columns={"diagnosis": "group"})
    metadata["sample"] = metadata["!Sample_title"].str.extract(r"(S\d+)$", expand=False)
    metadata["group"] = metadata["group"].str.upper().str.replace("MASH", "MASH", regex=False)
    return np.log2(matrix + 1), metadata[["sample", "group"]]


def load_gse83452() -> tuple[pd.DataFrame, pd.DataFrame]:
    path = DATA / "GSE83452_series_matrix.txt.gz"
    metadata, _ = parse_geo_matrix_metadata(path)
    probes = parse_series_table(path)
    matrix = collapse_symbols(probes, probe_to_symbol())
    metadata = metadata.rename(columns={"!Sample_geo_accession": "sample"})
    metadata["sample"] = metadata["sample"].astype(str)
    metadata["group"] = metadata["liver_status"].str.upper().str.replace("NO NASH", "NO_NASH", regex=False)
    metadata["timepoint"] = metadata["time"].str.lower().str.replace(" ", "_", regex=False)
    metadata["intervention"] = metadata["type_of_intervention"].str.upper()
    metadata["patient_key"] = metadata["sample_name"].str.extract(r"\((\d+)\)", expand=False)
    fallback = metadata["sample_name"].str.extract(r"patient\s+(\d+)", expand=False)
    metadata["patient_key"] = metadata["patient_key"].fillna(fallback)
    return matrix, metadata[["sample", "group", "timepoint", "intervention", "patient_key"]]


def paired_gse83452(output: pd.DataFrame) -> pd.DataFrame:
    rows = []
    baseline = output[output.timepoint.eq("baseline")].set_index("patient_key")
    follow = output[output.timepoint.eq("follow-up")].set_index("patient_key")
    common = sorted(set(baseline.index) & set(follow.index))
    for feature in TARGETS + list(SETS):
        pairs = pd.DataFrame({"baseline": baseline.loc[common, feature], "follow_up": follow.loc[common, feature]}).dropna()
        if len(pairs) >= 5:
            try:
                _, p = wilcoxon(pairs["follow_up"], pairs["baseline"])
            except ValueError:
                p = np.nan
            rows.append({"accession": "GSE83452", "analysis": "paired_followup_minus_baseline", "feature": feature, "n_pairs": len(pairs), "median_delta": np.median(pairs["follow_up"] - pairs["baseline"]), "p_value": p})
    return pd.DataFrame(rows)


def paired_gse193066() -> pd.DataFrame:
    filename = DATA / "GSE193066_NAFLD.HUn164.gct.gz"
    raw = pd.read_csv(filename, sep="\t", compression="gzip", skiprows=2)
    raw = raw.rename(columns={raw.columns[0]: "gene_name", raw.columns[1]: "gene_id"}).set_index("gene_name").drop(columns=["gene_id"])
    raw.index = raw.index.astype(str).str.upper()
    matrix = np.log2(raw.astype(float) + 1)
    meta = pd.read_csv(META / "GSE193066_sample_metadata.tsv", sep="\t")
    meta["sample"] = meta["!Sample_title"].astype(str)
    for col in meta.columns:
        if "characteristics" in col:
            values = meta[col].astype(str).str.split(":", n=1, expand=True)
            if values.shape[1] == 2:
                key = values.iloc[0, 0].strip().lower().replace(" ", "_")
                meta[key] = values.iloc[:, 1].str.strip().to_numpy()
    meta["pair"] = meta["sample"].str.replace(r"_(1|2)$", "", regex=True)
    meta["visit"] = np.where(meta["sample"].str.endswith("_2"), 2, 1)
    meta["fibrosis_numeric"] = pd.to_numeric(meta["fibrosis_stage"], errors="coerce")
    values = add_scores(matrix, meta[["sample", "pair", "visit", "fibrosis_numeric"]])
    first = values[values.visit.eq(1)].set_index("pair")
    second = values[values.visit.eq(2)].set_index("pair")
    common = sorted(set(first.index) & set(second.index))
    rows = []
    for feature in TARGETS + list(SETS):
        pair_values = pd.DataFrame({"first": first.loc[common, feature], "second": second.loc[common, feature], "fibrosis_first": first.loc[common, "fibrosis_numeric"], "fibrosis_second": second.loc[common, "fibrosis_numeric"]}).dropna(subset=["first", "second"])
        delta = pair_values["second"] - pair_values["first"]
        fibrosis_delta = pair_values["fibrosis_second"] - pair_values["fibrosis_first"]
        valid = pd.DataFrame({"delta": delta, "fibrosis_delta": fibrosis_delta}).dropna()
        rho, p = (spearmanr(valid["delta"], valid["fibrosis_delta"]) if len(valid) >= 10 and valid["fibrosis_delta"].nunique() > 1 else (np.nan, np.nan))
        rows.append({"accession": "GSE193066_HUn164", "analysis": "paired_visit2_minus_visit1", "feature": feature, "n_pairs": len(pair_values), "median_delta": np.median(delta) if len(delta) else np.nan, "delta_fibrosis_rho": rho, "delta_fibrosis_p": p})
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    all_tests = []

    matrix, metadata = load_gse281797()
    output = add_scores(matrix, metadata)
    output.reset_index().to_csv(OUT / "GSE281797_scores.tsv", sep="\t", index=False)
    all_tests.extend([test_groups("GSE281797", output, "group", "MASH", "MASL"), test_groups("GSE281797", output, "group", "MASH", "NO PATHOLOGY")])

    matrix, metadata = load_gse83452()
    output = add_scores(matrix, metadata)
    output.reset_index().to_csv(OUT / "GSE83452_scores.tsv", sep="\t", index=False)
    all_tests.append(test_groups("GSE83452", output[output.timepoint.eq("baseline")], "group", "NASH", "NO_NASH"))
    paired = paired_gse83452(output)
    if len(paired):
        paired["fdr_within_analysis"] = multipletests(paired["p_value"], method="fdr_bh")[1]
    paired.to_csv(OUT / "GSE83452_paired_followup.tsv", sep="\t", index=False)
    all_tests.append(paired)

    paired193 = paired_gse193066()
    valid_delta = paired193["delta_fibrosis_p"].notna()
    if valid_delta.any():
        paired193.loc[valid_delta, "delta_fibrosis_fdr"] = multipletests(paired193.loc[valid_delta, "delta_fibrosis_p"], method="fdr_bh")[1]
    paired193.to_csv(OUT / "GSE193066_HUn164_paired.tsv", sep="\t", index=False)
    all_tests.append(paired193)

    tests = pd.concat([frame for frame in all_tests if len(frame)], ignore_index=True, sort=False)
    if "p_value" in tests:
        mask = tests["p_value"].notna()
        tests.loc[mask, "fdr_within_analysis"] = tests.loc[mask].groupby(tests.loc[mask, "analysis"])["p_value"].transform(lambda x: multipletests(x, method="fdr_bh")[1])
    tests.to_csv(OUT / "additional_human_validation.tsv", sep="\t", index=False)
    print(tests[tests.feature.isin(["SLAMF8", "SPP1", "FARG95", "Iron_homeostasis", "HSC_activation"])].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
