#!/usr/bin/env python
import os
"""Assemble reviewed analysis panels without re-running any biological tests."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt


ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
FIG = ROOT / "figures"
ORDER = ["GSE344087", "GSE298719", "GSE212837"]
PREFIX = {dataset: f"{index:03d}" for index, dataset in enumerate(ORDER, start=1)}
PANELS = {
    "Figure1_myeloid_HSC_atlas": [
        [(f"{PREFIX[ds]}A_{ds}_myeloid_UMAP_curated.png", f"{ds}: myeloid") for ds in ORDER],
        [(f"{PREFIX[ds]}B_{ds}_hsc_UMAP_curated.png", f"{ds}: HSC") for ds in ORDER],
    ],
    "Figure2_SLAMF8_macrophage_state": [
        [("008_cross_dataset_SLAMF8_high_state.png", "Cross-cohort state"),
         ("010_cross_dataset_macrophage_SLAMF8_HSC_activation.png", "Donor-matched HSC association"), None],
        [(f"007_{ds}_SLAMF8_high_paired_effects.png", f"{ds}: within-donor states") for ds in ORDER],
    ],
    "Figure3_macrophage_trajectories": [
        [(f"004A_{ds}_monocle3_subtypes.png", f"{ds}: inferred trajectory") for ds in ORDER],
        [(f"004B_{ds}_monocle3_pseudotime.png", f"{ds}: pseudotime") for ds in ORDER],
    ],
    "Figure4_macrophage_HSC_communication_audit": [
        [("016_cross_cohort_SPP1_review.png", "Condition-adjusted cross-cohort SPP1 review")],
    ],
    "FigureS2_GSE298719_communication_audit": [
        [("015_communication_audit.png", "GSE298719 patient-level communication audit")],
    ],
    "FigureS1_NicheNet_ligand_prioritization": [
        [(f"009_{ds}_NicheNet_ligand_activity.png", f"{ds}: exploratory ligand activity") for ds in ORDER],
    ],
    "Figure5_independent_bulk_validation": [
        [("017_expanded_public_validation.png", "Expanded public-cohort validation and longitudinal robustness")],
    ],
}


def main():
    for name, rows in PANELS.items():
        columns = max(len(row) for row in rows)
        figure_width = 15.6 if columns == 1 else 6.1 * columns
        fig, axes = plt.subplots(len(rows), columns, figsize=(figure_width, 4.7 * len(rows)),
                                 squeeze=False, constrained_layout=True)
        for row_index, row in enumerate(rows):
            for column_index in range(columns):
                ax = axes[row_index, column_index]
                if column_index >= len(row) or row[column_index] is None:
                    ax.set_axis_off()
                    continue
                filename, title = row[column_index]
                path = FIG / filename
                if not path.is_file():
                    raise FileNotFoundError(path)
                ax.imshow(mpimg.imread(path))
                ax.set_title(title, fontsize=13, loc="left")
                ax.set_axis_off()
        fig.savefig(FIG / f"{name}.pdf", dpi=200, bbox_inches="tight")
        fig.savefig(FIG / f"{name}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(name, flush=True)


if __name__ == "__main__":
    main()
