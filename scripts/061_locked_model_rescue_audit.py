#!/usr/bin/env python3
"""Transparent rescue audit for the locked external transcriptomic pilot.

This script does not select predictors or tune hyperparameters against the
external cohort. It quantifies model comparison, calibration-independent
discrimination, training-derived operating points, and explicitly exploratory
incremental value over available age/sex/NAS covariates.
"""
from __future__ import annotations
import os

from pathlib import Path
import json

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    confusion_matrix,
    roc_auc_score,
    roc_curve,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
SOURCE = ROOT / "results" / "revision_p2_5" / "bulk_immune_and_models"
OUT = ROOT / "results" / "revision_p2_6" / "model_rescue_audit"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260928
N_BOOT = 10000


def fixed_logistic() -> Pipeline:
    return Pipeline(
        [
            ("scale", StandardScaler()),
            (
                "logistic",
                LogisticRegression(
                    C=1.0,
                    penalty="l2",
                    solver="lbfgs",
                    max_iter=500,
                    random_state=SEED,
                ),
            ),
        ]
    )


def metric_row(name: str, y: np.ndarray, probability: np.ndarray) -> dict:
    prevalence = float(np.mean(y))
    brier = float(brier_score_loss(y, probability))
    null_brier = prevalence * (1.0 - prevalence)
    return {
        "model": name,
        "n": int(len(y)),
        "events": int(np.sum(y)),
        "prevalence": prevalence,
        "auc": float(roc_auc_score(y, probability)),
        "pr_auc": float(average_precision_score(y, probability)),
        "pr_auc_to_prevalence_ratio": float(average_precision_score(y, probability) / prevalence),
        "brier": brier,
        "null_brier": null_brier,
        "brier_skill_score": float(1.0 - brier / null_brier),
    }


def bootstrap_comparison(y: np.ndarray, probabilities: dict[str, np.ndarray]) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    names = list(probabilities)
    rows = []
    for left_index, left in enumerate(names):
        for right in names[left_index + 1 :]:
            observed_auc = roc_auc_score(y, probabilities[left]) - roc_auc_score(y, probabilities[right])
            observed_brier = brier_score_loss(y, probabilities[left]) - brier_score_loss(y, probabilities[right])
            auc_differences = []
            brier_differences = []
            for _ in range(N_BOOT):
                index = rng.integers(0, len(y), len(y))
                if np.unique(y[index]).size < 2:
                    continue
                auc_differences.append(
                    roc_auc_score(y[index], probabilities[left][index])
                    - roc_auc_score(y[index], probabilities[right][index])
                )
                brier_differences.append(
                    brier_score_loss(y[index], probabilities[left][index])
                    - brier_score_loss(y[index], probabilities[right][index])
                )
            rows.append(
                {
                    "left_model": left,
                    "right_model": right,
                    "auc_difference_left_minus_right": observed_auc,
                    "auc_difference_ci_low": float(np.quantile(auc_differences, 0.025)),
                    "auc_difference_ci_high": float(np.quantile(auc_differences, 0.975)),
                    "brier_difference_left_minus_right": observed_brier,
                    "brier_difference_ci_low": float(np.quantile(brier_differences, 0.025)),
                    "brier_difference_ci_high": float(np.quantile(brier_differences, 0.975)),
                }
            )
    return pd.DataFrame(rows)


def training_threshold_operating_point(
    model_name: str,
    train_y: np.ndarray,
    train_probability: np.ndarray,
    test_y: np.ndarray,
    test_probability: np.ndarray,
) -> dict:
    false_positive_rate, true_positive_rate, thresholds = roc_curve(train_y, train_probability)
    youden_index = int(np.argmax(true_positive_rate - false_positive_rate))
    threshold = float(thresholds[youden_index])
    predicted = (test_probability >= threshold).astype(int)
    true_negative, false_positive, false_negative, true_positive = confusion_matrix(
        test_y, predicted, labels=[0, 1]
    ).ravel()
    return {
        "model": model_name,
        "threshold_source": "training_youden",
        "threshold": threshold,
        "external_sensitivity": true_positive / (true_positive + false_negative),
        "external_specificity": true_negative / (true_negative + false_positive),
        "external_ppv": true_positive / (true_positive + false_positive) if true_positive + false_positive else np.nan,
        "external_npv": true_negative / (true_negative + false_negative) if true_negative + false_negative else np.nan,
        "external_accuracy": (true_positive + true_negative) / len(test_y),
        "external_true_positive": int(true_positive),
        "external_false_positive": int(false_positive),
        "external_true_negative": int(true_negative),
        "external_false_negative": int(false_negative),
    }


def main() -> None:
    scores = pd.read_csv(SOURCE / "bulk_cell_state_signature_scores.tsv", sep="\t")
    locked_predictions = pd.read_csv(
        SOURCE / "locked_external_transcriptomic_model_predictions.tsv", sep="\t"
    )
    scores = scores.loc[~scores["normal_histology"].astype(bool)].copy()
    scores["F_ge_2"] = (scores["fibrosis"] >= 2).astype(int)

    external_locked = locked_predictions.loc[
        locked_predictions["split"].eq("external_test")
    ].copy()
    external_wide = external_locked.pivot(
        index=["sample", "F_ge_2"], columns="model", values="probability"
    ).reset_index()
    external_y = external_wide["F_ge_2"].to_numpy(dtype=int)

    metrics = []
    locked_probabilities = {}
    for model_name in ["SLAMF8_only", "Locked_4_feature_transcriptomic"]:
        probability = external_wide[model_name].to_numpy(dtype=float)
        locked_probabilities[model_name] = probability
        metrics.append(metric_row(model_name, external_y, probability))

    train_locked = locked_predictions.loc[
        locked_predictions["split"].eq("training_apparent")
    ].copy()
    operating_points = []
    for model_name in ["SLAMF8_only", "Locked_4_feature_transcriptomic"]:
        train_model = train_locked.loc[train_locked["model"].eq(model_name)]
        test_model = external_locked.loc[external_locked["model"].eq(model_name)]
        operating_points.append(
            training_threshold_operating_point(
                model_name,
                train_model["F_ge_2"].to_numpy(dtype=int),
                train_model["probability"].to_numpy(dtype=float),
                test_model["F_ge_2"].to_numpy(dtype=int),
                test_model["probability"].to_numpy(dtype=float),
            )
        )

    incremental_specs = {
        "Age_sex": ["age_numeric", "sex_female"],
        "Age_sex_NAS": ["age_numeric", "sex_female", "nas_numeric"],
        "Age_sex_NAS_plus_SLAMF8": ["age_numeric", "sex_female", "nas_numeric", "SLAMF8"],
        "Age_sex_NAS_plus_4_state_features": [
            "age_numeric",
            "sex_female",
            "nas_numeric",
            "SLAMF8",
            "SAMac_like",
            "SPP1",
            "Activated_HSC_fibrogenic",
        ],
    }
    incremental_predictions = {}
    incremental_metrics = []
    coefficient_rows = []
    for model_name, fields in incremental_specs.items():
        train = scores.loc[scores["dataset"].eq("GSE162694")].dropna(
            subset=fields + ["F_ge_2"]
        )
        test = scores.loc[scores["dataset"].eq("GSE130970")].dropna(
            subset=fields + ["F_ge_2"]
        )
        model = fixed_logistic()
        model.fit(train[fields], train["F_ge_2"])
        probability = model.predict_proba(test[fields])[:, 1]
        incremental_predictions[model_name] = probability
        incremental_metrics.append(metric_row(model_name, test["F_ge_2"].to_numpy(), probability))
        coefficients = model.named_steps["logistic"].coef_[0]
        coefficient_rows.extend(
            {"model": model_name, "feature": field, "standardized_log_odds_coefficient": coefficient}
            for field, coefficient in zip(fields, coefficients)
        )

    incremental_y = scores.loc[scores["dataset"].eq("GSE130970")].dropna(
        subset=incremental_specs["Age_sex_NAS_plus_4_state_features"] + ["F_ge_2"]
    )["F_ge_2"].to_numpy(dtype=int)

    pd.DataFrame(metrics).to_csv(OUT / "locked_model_external_value.tsv", sep="\t", index=False)
    bootstrap_comparison(external_y, locked_probabilities).to_csv(
        OUT / "locked_model_paired_bootstrap_comparison.tsv", sep="\t", index=False
    )
    pd.DataFrame(operating_points).to_csv(
        OUT / "training_threshold_external_operating_points.tsv", sep="\t", index=False
    )
    pd.DataFrame(incremental_metrics).to_csv(
        OUT / "exploratory_covariate_incremental_models.tsv", sep="\t", index=False
    )
    bootstrap_comparison(incremental_y, incremental_predictions).to_csv(
        OUT / "exploratory_incremental_paired_bootstrap.tsv", sep="\t", index=False
    )
    pd.DataFrame(coefficient_rows).to_csv(
        OUT / "exploratory_incremental_model_coefficients.tsv", sep="\t", index=False
    )

    feature_columns = ["SLAMF8", "SAMac_like", "SPP1", "Activated_HSC_fibrogenic"]
    train_features = scores.loc[scores["dataset"].eq("GSE162694"), feature_columns]
    test_features = scores.loc[scores["dataset"].eq("GSE130970"), feature_columns]
    train_features.corr(method="spearman").to_csv(
        OUT / "training_feature_spearman_correlation.tsv", sep="\t"
    )
    test_features.corr(method="spearman").to_csv(
        OUT / "external_feature_spearman_correlation.tsv", sep="\t"
    )

    manifest = {
        "purpose": "post-result rescue audit without external-cohort tuning or predictor selection",
        "primary_locked_models_unchanged": True,
        "incremental_models_status": "exploratory post hoc sensitivity analysis",
        "external_cohort": "GSE130970",
        "endpoint": "fibrosis stage F>=2",
        "bootstrap_replicates": N_BOOT,
        "seed": SEED,
    }
    (OUT / "README.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print("Locked model external value")
    print(pd.DataFrame(metrics).to_string(index=False))
    print("\nTraining-derived operating points")
    print(pd.DataFrame(operating_points).to_string(index=False))
    print("\nExploratory covariate incremental models")
    print(pd.DataFrame(incremental_metrics).to_string(index=False))


if __name__ == "__main__":
    main()
