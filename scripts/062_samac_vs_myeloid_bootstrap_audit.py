#!/usr/bin/env python3
import os
"""Compare paired fibrosis associations of two fixed bulk signatures."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
SCORES = ROOT / "results" / "revision_p2_5" / "bulk_immune_and_models" / "bulk_cell_state_signature_scores.tsv"
OUTPUT = ROOT / "results" / "revision_p2_5" / "samac_vs_myeloid_bootstrap_audit.tsv"
SEED = 20260923
BOOTSTRAPS = 5000


def correlation(frame: pd.DataFrame, feature: str) -> float:
    return float(spearmanr(frame["fibrosis"], frame[feature]).statistic)


def main() -> None:
    scores = pd.read_csv(SCORES, sep="\t")
    rows = []
    for dataset, cohort in scores.groupby("dataset", sort=True):
        cohort = cohort.loc[~cohort["normal_histology"].astype(bool)].copy()
        for signature in ("SAMac_like", "SAMac_like_no_SPP1"):
            columns = ["sample", "fibrosis", signature, "Macrophage_myeloid"]
            paired = cohort[columns].dropna()
            if paired["sample"].duplicated().any():
                raise ValueError(f"Duplicate biopsy identifiers in {dataset}")
            observed_signature = correlation(paired, signature)
            observed_myeloid = correlation(paired, "Macrophage_myeloid")
            differences = []
            rng = np.random.default_rng(SEED)
            for _ in range(BOOTSTRAPS):
                indices = rng.integers(0, len(paired), len(paired))
                sampled = paired.iloc[indices]
                if sampled["fibrosis"].nunique() < 2:
                    continue
                difference = correlation(sampled, signature) - correlation(sampled, "Macrophage_myeloid")
                if np.isfinite(difference):
                    differences.append(difference)
            lower, upper = np.quantile(differences, [0.025, 0.975])
            rows.append(
                {
                    "dataset": dataset,
                    "signature": signature,
                    "comparator": "Macrophage_myeloid",
                    "n_biopsies": len(paired),
                    "rho_signature": observed_signature,
                    "rho_comparator": observed_myeloid,
                    "delta_rho": observed_signature - observed_myeloid,
                    "bootstrap_ci_lower": lower,
                    "bootstrap_ci_upper": upper,
                    "n_valid_bootstraps": len(differences),
                    "seed": SEED,
                    "resampling_unit": "biopsy",
                }
            )
    output = pd.DataFrame(rows)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUTPUT, sep="\t", index=False)
    print(output.to_string(index=False))


if __name__ == "__main__":
    main()
