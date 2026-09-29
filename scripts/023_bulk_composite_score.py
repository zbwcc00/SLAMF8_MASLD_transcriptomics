#!/usr/bin/env python
import os
"""Pre-specified SLAMF8/ferro-aging composite score evaluation."""

from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve


PROJECT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
ROOT = PROJECT / "results" / "bulk_validation"


def z(series: pd.Series) -> pd.Series:
    return (series - series.mean()) / series.std(ddof=0)


def evaluate(data: pd.DataFrame, dataset: str, outcome_name: str, outcome: pd.Series) -> pd.DataFrame:
    features = ["SLAMF8", "FARG95", "Iron_homeostasis", "HSC_activation"]
    rows = []
    valid = data[features].notna().all(axis=1) & outcome.notna()
    data = data.loc[valid].copy()
    outcome = outcome.loc[valid].astype(int)
    scaled = pd.DataFrame({feature: z(data[feature]) for feature in features}, index=data.index)
    composites = {
        "SLAMF8_only": scaled["SLAMF8"],
        "FARG95_only": scaled["FARG95"],
        "Iron_only": scaled["Iron_homeostasis"],
        "SLAMF8_FARG95_Iron": scaled[["SLAMF8", "FARG95", "Iron_homeostasis"]].mean(axis=1),
        "SLAMF8_FARG95_Iron_HSC": scaled.mean(axis=1),
    }
    for name, score in composites.items():
        auc = roc_auc_score(outcome, score) if outcome.nunique() == 2 else np.nan
        fpr, tpr, thresholds = roc_curve(outcome, score) if outcome.nunique() == 2 else ([], [], [])
        j = int(np.argmax(tpr - fpr)) if len(tpr) else 0
        rows.append({"dataset": dataset, "outcome": outcome_name, "n": len(outcome), "positive_n": int(outcome.sum()), "score": name, "auc": auc, "youden_sensitivity": float(tpr[j]) if len(tpr) else np.nan, "youden_specificity": float(1 - fpr[j]) if len(fpr) else np.nan, "direction_positive": int(score[outcome == 1].mean() > score[outcome == 0].mean())})
    score_table = scaled.copy()
    for name, score in composites.items():
        score_table[name] = score
    score_table.insert(0, "outcome", outcome)
    score_table.to_csv(ROOT / f"{dataset}_{outcome_name}_composite_scores.tsv", sep="\t")
    return pd.DataFrame(rows)


def main() -> None:
    rows = []
    gse48452 = pd.read_csv(ROOT / "GSE48452_bulk_scores.tsv", sep="\t")
    fibrosis = pd.to_numeric(gse48452["fibrosis"], errors="coerce")
    rows.append(evaluate(gse48452, "GSE48452", "fibrosis_ge1", fibrosis.ge(1).where(fibrosis.notna())))
    rows.append(evaluate(gse48452, "GSE48452", "NASH_vs_other", (gse48452["group"].str.lower() == "nash").astype(int)))
    gse63067 = pd.read_csv(ROOT / "GSE63067_bulk_scores.tsv", sep="\t")
    rows.append(evaluate(gse63067, "GSE63067", "NASH_vs_other", (gse63067["disease status"].str.lower() == "non-alcoholic steatohepatitis").astype(int)))
    result = pd.concat(rows, ignore_index=True)
    result.to_csv(ROOT / "bulk_composite_score_performance.tsv", sep="\t", index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
