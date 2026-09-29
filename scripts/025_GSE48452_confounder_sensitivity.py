#!/usr/bin/env python
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PROJECT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
ROOT = PROJECT / "results" / "bulk_validation"


def main() -> None:
    data = pd.read_csv(ROOT / "GSE48452_bulk_scores.tsv", sep="\t")
    data["fibrosis_numeric"] = pd.to_numeric(data["fibrosis"], errors="coerce")
    rows = []
    filters = {
        "all_samples": data,
        "exclude_after_surgery": data[data["bariatric surgery"].fillna("").ne("after surgery")],
        "no_surgery_label_only": data[data["bariatric surgery"].isna()],
    }
    for name, subset in filters.items():
        valid = subset.dropna(subset=["fibrosis_numeric"])
        for target in ["SLAMF8", "FARG95", "Iron_homeostasis", "Senescence_SASP", "HSC_activation"]:
            rho, p = spearmanr(valid["fibrosis_numeric"], valid[target])
            rows.append({"filter": name, "target": target, "n": len(valid), "rho": rho, "p_value": p})
    result = pd.DataFrame(rows)
    result.to_csv(ROOT / "GSE48452_confounder_sensitivity.tsv", sep="\t", index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
