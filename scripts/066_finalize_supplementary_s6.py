#!/usr/bin/env python3
import os
"""Finalize the exploratory NicheNet ligand-prioritization supplement."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
FIG = ROOT / "figures"
OUT = FIG / "revision_p0"
SOURCE = ROOT / "source_data"
OUT.mkdir(parents=True, exist_ok=True)
SOURCE.mkdir(parents=True, exist_ok=True)
DATASETS = ["GSE344087", "GSE298719", "GSE212837"]


def main() -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.6, 5.0), dpi=300, constrained_layout=True)
    for index, (axis, dataset) in enumerate(zip(axes, DATASETS)):
        path = FIG / f"009_{dataset}_NicheNet_ligand_activity.png"
        axis.imshow(mpimg.imread(path))
        axis.set_axis_off()
        axis.set_title(f"{dataset}: exploratory ligand activity", loc="left", fontsize=9.2, pad=5)
        axis.text(0.0, 1.03, chr(65 + index), transform=axis.transAxes,
                  fontsize=9, fontweight="bold", va="bottom")
    fig.suptitle("NicheNet ligand prioritization across liver cohorts",
                 fontsize=11, fontweight="bold", y=1.08)
    fig.text(0.01, 0.01,
             "Expression-based ligand–target ranking; exploratory and not causal evidence.",
             fontsize=7.2, color="#566573")
    stem = OUT / "Supplementary_Figure6_NicheNet_ligand_prioritization"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    rows = []
    for dataset in DATASETS:
        for prefix in ("027", "028"):
            matches = sorted((ROOT / "results" / "communication").glob(
                f"{prefix}_{dataset}_NicheNet_*.tsv"))
            for path in matches:
                table = pd.read_csv(path, sep="\t")
                table.insert(0, "dataset", dataset)
                table.insert(1, "source_file", path.name)
                rows.append(table)
    if rows:
        pd.concat(rows, ignore_index=True, sort=False).to_csv(
            SOURCE / "Source_Data_Supplementary_Figure6_NicheNet_ligand_prioritization.tsv.gz",
            sep="\t", index=False, compression="gzip")


if __name__ == "__main__":
    main()
