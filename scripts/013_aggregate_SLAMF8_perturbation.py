#!/usr/bin/env python
import os
"""Aggregate cell-level SLAMF8 counterfactuals at patient and cohort level."""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp, combine_pvalues
from statsmodels.stats.multitest import multipletests


PROJECT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
OUT = PROJECT / "results" / "perturbation"


def main() -> None:
    files = sorted(OUT.glob("031_*_SLAMF8_network_perturbation_cell.tsv.gz"))
    patient_tables = []
    for path in files:
        data = pd.read_csv(path, sep="\t")
        delta_cols = [x for x in data.columns if x.startswith("delta_")]
        long = data.melt(id_vars=["dataset", "cell", "donor", "subtype"], value_vars=delta_cols, var_name="target", value_name="delta")
        long["target"] = long["target"].str.replace("^delta_", "", regex=True)
        patient = (long.groupby(["dataset", "donor", "target"], as_index=False)
                   .agg(delta=("delta", "mean"), cell_n=("delta", "size")))
        patient_tables.append(patient)
    patients = pd.concat(patient_tables, ignore_index=True)
    patients.to_csv(OUT / "034_SLAMF8_network_perturbation_patient_level.tsv", sep="\t", index=False)

    rows = []
    for (dataset, target), group in patients.groupby(["dataset", "target"]):
        values = group.delta.to_numpy(dtype=float)
        p = float(ttest_1samp(values, 0.0, nan_policy="omit").pvalue) if len(values) >= 2 else float("nan")
        rows.append({"dataset": dataset, "target": target, "patient_n": len(values), "mean_delta": float(np.mean(values)), "median_delta": float(np.median(values)), "sd_delta": float(np.std(values, ddof=1)) if len(values) > 1 else float("nan"), "one_sample_p": p, "direction": "down" if np.mean(values) < 0 else "up" if np.mean(values) > 0 else "flat"})
    cohort = pd.DataFrame(rows)
    cohort.to_csv(OUT / "035_SLAMF8_network_perturbation_dataset_summary.tsv", sep="\t", index=False)

    meta_rows = []
    for target, group in cohort.groupby("target"):
        pvals = group.one_sample_p.dropna().to_numpy()
        signs = group.direction.tolist()
        meta_rows.append({
            "target": target,
            "datasets": len(group),
            "mean_of_dataset_delta": float(group.mean_delta.mean()),
            "datasets_down": int(sum(x == "down" for x in signs)),
            "datasets_up": int(sum(x == "up" for x in signs)),
            "direction_consistent": bool(len(set(signs)) == 1),
            "fisher_p": float(combine_pvalues(pvals, method="fisher")[1]) if len(pvals) else float("nan"),
        })
    meta = pd.DataFrame(meta_rows)
    meta["fisher_fdr"] = multipletests(meta.fisher_p.fillna(1.0).to_numpy(), method="fdr_bh")[1]
    meta.sort_values(["fisher_fdr", "mean_of_dataset_delta"], key=lambda x: x.abs() if x.name == "mean_of_dataset_delta" else x, inplace=True)
    meta.to_csv(OUT / "036_cross_dataset_SLAMF8_network_perturbation_meta.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
