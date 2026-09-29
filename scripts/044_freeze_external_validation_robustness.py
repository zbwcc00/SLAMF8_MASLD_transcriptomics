#!/usr/bin/env python
"""Freeze external validation inputs and run patient/sample-unit robustness audits."""

from __future__ import annotations
import os

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
RESULTS = ROOT / "results" / "bulk_validation" / "expanded_public"
FREEZE = ROOT / "results" / "final_freeze"
FIGURES = ROOT / "figures"
FEATURES = ["SLAMF8", "SPP1", "FARG95", "Iron_homeostasis", "HSC_activation"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def freeze_manifest() -> pd.DataFrame:
    paths = [
        RESULTS / "expanded_public_validation.tsv",
        RESULTS / "additional_human_validation.tsv",
        RESULTS / "GSE167523_scores.tsv",
        RESULTS / "GSE126848_scores.tsv",
        RESULTS / "GSE193066_NAFLD.HUn106_scores.tsv",
        RESULTS / "GSE193066_NAFLD.HUn164_scores.tsv",
        RESULTS / "GSE193066_HUn164_paired.tsv",
        RESULTS / "GSE83452_scores.tsv",
        RESULTS / "GSE83452_paired_followup.tsv",
        RESULTS / "GSE281797_scores.tsv",
        ROOT / "scripts" / "039_validate_expanded_public_cohorts.py",
        ROOT / "scripts" / "042_validate_additional_human_cohorts.py",
    ]
    rows = []
    for path in paths:
        rows.append({"path": str(path), "exists": path.exists(), "size_bytes": path.stat().st_size if path.exists() else np.nan, "sha256": sha256(path) if path.exists() else ""})
    return pd.DataFrame(rows)


def bootstrap_auc(data: pd.DataFrame, label: str, case: str, control: str, n_boot: int = 2000) -> pd.DataFrame:
    rng = np.random.default_rng(20260921)
    subset = data[data[label].isin([case, control])].copy()
    rows = []
    y = (subset[label] == case).astype(int).to_numpy()
    for feature in FEATURES:
        values = pd.to_numeric(subset[feature], errors="coerce").to_numpy()
        keep = np.isfinite(values) & np.isfinite(y)
        values, labels = values[keep], y[keep]
        observed = roc_auc_score(labels, values) if len(np.unique(labels)) == 2 else np.nan
        boot = []
        for _ in range(n_boot):
            indices = rng.integers(0, len(values), len(values))
            sample_y, sample_values = labels[indices], values[indices]
            if len(np.unique(sample_y)) == 2:
                boot.append(roc_auc_score(sample_y, sample_values))
        rows.append({"analysis": f"{case}_vs_{control}", "feature": feature, "n": len(values), "auc": observed, "auc_ci_low": np.quantile(boot, 0.025) if boot else np.nan, "auc_ci_high": np.quantile(boot, 0.975) if boot else np.nan})
    return pd.DataFrame(rows)


def direction_audit() -> pd.DataFrame:
    tables = []
    for path in [RESULTS / "expanded_public_validation.tsv", RESULTS / "additional_human_validation.tsv"]:
        table = pd.read_csv(path, sep="\t")
        table["source_file"] = path.name
        tables.append(table)
    data = pd.concat(tables, ignore_index=True, sort=False)
    rows = []
    for _, row in data[data.feature.isin(FEATURES)].iterrows():
        if pd.notna(row.get("median_difference")):
            effect, metric = row["median_difference"], "median_difference"
        elif pd.notna(row.get("rho")):
            effect, metric = row["rho"], "rho"
        elif pd.notna(row.get("median_delta")):
            effect, metric = row["median_delta"], "median_delta"
        elif pd.notna(row.get("delta_fibrosis_rho")):
            effect, metric = row["delta_fibrosis_rho"], "delta_fibrosis_rho"
        else:
            continue
        rows.append({"accession": row.get("accession", ""), "analysis": row.get("analysis", ""), "feature": row["feature"], "effect_metric": metric, "effect": effect, "direction": "positive" if effect > 0 else "negative" if effect < 0 else "zero", "p_value": row.get("p_value", row.get("delta_fibrosis_p", np.nan)), "source_file": row["source_file"]})
    return pd.DataFrame(rows)


def patient_unit_audit() -> pd.DataFrame:
    rows = []
    for filename, unit, expected in [
        ("GSE167523_scores.tsv", "independent biopsy", 98),
        ("GSE126848_scores.tsv", "scored samples", 31),
        ("GSE193066_NAFLD.HUn106_scores.tsv", "first biopsy patients", 106),
        ("GSE193066_HUn164_paired.tsv", "paired patients", 58),
        ("GSE83452_paired_followup.tsv", "paired patients", 60),
        ("GSE281797_scores.tsv", "independent biopsy", 94),
    ]:
        path = RESULTS / filename
        if "paired" in filename:
            observed = int(pd.read_csv(path, sep="\t")["n_pairs"].max())
        else:
            observed = len(pd.read_csv(path, sep="\t"))
        rows.append({"file": filename, "unit": unit, "observed_n": observed, "expected_n": expected, "status": "PASS" if observed == expected else "CHECK"})
    return pd.DataFrame(rows)


def main() -> None:
    FREEZE.mkdir(parents=True, exist_ok=True)
    manifest = freeze_manifest()
    manifest.to_csv(FREEZE / "external_validation_freeze_manifest.tsv", sep="\t", index=False)

    gse167523 = pd.read_csv(RESULTS / "GSE167523_scores.tsv", sep="\t")
    gse281797 = pd.read_csv(RESULTS / "GSE281797_scores.tsv", sep="\t")
    gse83452 = pd.read_csv(RESULTS / "GSE83452_scores.tsv", sep="\t")
    auc = pd.concat([
        bootstrap_auc(gse167523, "group", "NASH", "NAFL"),
        bootstrap_auc(gse281797, "group", "MASH", "MASL"),
        bootstrap_auc(gse83452[gse83452.timepoint.eq("baseline")], "group", "NASH", "NO_NASH"),
    ], ignore_index=True)
    auc.to_csv(FREEZE / "external_validation_bootstrap_auc.tsv", sep="\t", index=False)

    directions = direction_audit()
    directions.to_csv(FREEZE / "external_validation_direction_audit.tsv", sep="\t", index=False)
    patient_unit_audit().to_csv(FREEZE / "patient_sample_unit_audit.tsv", sep="\t", index=False)

    paired = pd.read_csv(RESULTS / "GSE193066_HUn164_paired.tsv", sep="\t")
    valid = paired["delta_fibrosis_p"].notna()
    paired.loc[valid, "delta_fibrosis_fdr"] = multipletests(paired.loc[valid, "delta_fibrosis_p"], method="fdr_bh")[1]
    paired.to_csv(RESULTS / "GSE193066_HUn164_paired.tsv", sep="\t", index=False)
    print("freeze and robustness audit complete", flush=True)


if __name__ == "__main__":
    main()
