#!/usr/bin/env python
"""Re-audit virtual SLAMF8 perturbation with donor-level uncertainty.

The legacy Fisher P values are intentionally not reused.  Donor-level
sign-flip permutation and bootstrap intervals are descriptive sensitivity
statistics for a counterfactual prediction, not causal evidence.
"""
from __future__ import annotations
import os

from pathlib import Path
import itertools
import hashlib
import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
RES = ROOT / "results" / "perturbation"
OUT = ROOT / "results" / "revision_p0" / "virtual_perturbation"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260920


def sign_flip_p(values: np.ndarray, seed: int = SEED) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 2:
        return np.nan
    observed = abs(values.mean())
    if len(values) <= 20:
        all_signs = np.array(list(itertools.product([-1.0, 1.0], repeat=len(values))), dtype=float)
        perm = np.abs((all_signs * values).mean(axis=1))
    else:
        rng = np.random.default_rng(seed)
        perm = np.abs((rng.choice([-1.0, 1.0], size=(20000, len(values))) * values).mean(axis=1))
    return float((np.sum(perm >= observed) + 1) / (len(perm) + 1))


def bootstrap(values: np.ndarray, seed: int) -> tuple[float, float, float, float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 2:
        return (np.nan,) * 5
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(values), size=(20000, len(values)))
    means = values[indices].mean(axis=1)
    medians = np.median(values[indices], axis=1)
    return (float(values.mean()), float(np.quantile(means, .025)), float(np.quantile(means, .975)), float(np.quantile(medians, .025)), float(np.quantile(medians, .975)))


def main() -> None:
    patient = pd.read_csv(RES / "034_SLAMF8_network_perturbation_patient_level.tsv", sep="\t")
    summaries = []
    for (dataset, target), group in patient.groupby(["dataset", "target"], sort=False):
        values = group.delta.to_numpy(dtype=float)
        token = hashlib.sha256(f"{dataset}|{target}".encode()).hexdigest()
        mean, ci_low, ci_high, median_ci_low, median_ci_high = bootstrap(values, SEED + int(token[:8], 16) % 100000)
        summaries.append({
            "dataset": dataset,
            "target": target,
            "donor_n": len(values),
            "mean_delta": mean,
            "median_delta": float(np.median(values)) if len(values) else np.nan,
            "sd_delta": float(np.std(values, ddof=1)) if len(values) > 1 else np.nan,
            "standardized_mean_delta": mean / np.std(values, ddof=1) if len(values) > 1 and np.std(values, ddof=1) > 0 else np.nan,
            "fraction_down": float(np.mean(values < 0)) if len(values) else np.nan,
            "fraction_up": float(np.mean(values > 0)) if len(values) else np.nan,
            "sign_flip_p": sign_flip_p(values, SEED),
            "mean_bootstrap_ci_low": ci_low,
            "mean_bootstrap_ci_high": ci_high,
            "median_bootstrap_ci_low": median_ci_low,
            "median_bootstrap_ci_high": median_ci_high,
        })
    audit = pd.DataFrame(summaries)
    cv_frames = []
    for path in sorted(RES.glob("032_*_SLAMF8_network_perturbation_summary.tsv")):
        d = pd.read_csv(path, sep="\t")
        if "group_cv_r2_mean" in d:
            cv_frames.append(d[["dataset", "target", "group_cv_r2_mean", "n_donors"]])
    if cv_frames:
        cv = pd.concat(cv_frames, ignore_index=True).rename(columns={"n_donors": "model_donor_n"})
        audit = audit.merge(cv, on=["dataset", "target"], how="left")
    audit.to_csv(OUT / "donor_level_counterfactual_audit.tsv", sep="\t", index=False)
    rows = []
    for target, group in audit.groupby("target", sort=False):
        direction = np.sign(group.mean_delta.to_numpy(float))
        rows.append({
            "target": target,
            "n_datasets": len(group),
            "datasets_down": int(np.sum(direction < 0)),
            "datasets_up": int(np.sum(direction > 0)),
            "direction_consistent": bool(len(set(direction.tolist())) == 1),
            "median_standardized_mean_delta": float(group.standardized_mean_delta.median()),
            "median_group_cv_r2": float(group.group_cv_r2_mean.median()),
            "interpretation": "descriptive counterfactual audit; no causal P value",
        })
    cross = pd.DataFrame(rows)
    cross.to_csv(OUT / "cross_dataset_counterfactual_audit.tsv", sep="\t", index=False)
    (OUT / "README.md").write_text(
        "# Virtual perturbation re-audit\n\n"
        "The legacy Fisher P values are not used. Results are reported as donor-level "
        "predicted deltas, sign-flip permutation P values, bootstrap intervals, standardized "
        "effect sizes, direction consistency and donor-aware group-CV R2. These are model-based "
        "counterfactual predictions and are not experimental SLAMF8 knockout evidence.\n",
        encoding="utf-8",
    )
    print(cross.to_string(index=False))


if __name__ == "__main__":
    main()
