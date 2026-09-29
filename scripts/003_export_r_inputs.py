from __future__ import annotations
import os

import gzip
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse


PROJECT_ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
PROCESSED = PROJECT_ROOT / "data_processed"
RESULTS = PROJECT_ROOT / "results"


def write_bundle(
    counts: sparse.spmatrix,
    obs: pd.DataFrame,
    var: pd.DataFrame,
    prefix: Path,
    embeddings: dict[str, np.ndarray] | None = None,
) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    matrix = counts.T.tocoo()
    with gzip.open(prefix.with_name(prefix.name + "_counts.mtx.gz"), "wb") as handle:
        io.mmwrite(handle, matrix, field="integer")
    features = var.copy()
    features.insert(0, "var_name", features.index.astype(str))
    features.to_csv(prefix.with_name(prefix.name + "_features.tsv.gz"), sep="\t", index=False, compression="gzip")
    obs_to_write = obs.drop(columns=["barcode"], errors="ignore")
    obs_to_write.to_csv(prefix.with_name(prefix.name + "_cell_metadata.tsv.gz"), sep="\t", index=True, index_label="barcode", compression="gzip")
    if embeddings:
        for name, values in embeddings.items():
            columns = [f"{name}_{index + 1}" for index in range(values.shape[1])]
            frame = pd.DataFrame(values, index=obs.index, columns=columns)
            frame.to_csv(prefix.with_name(prefix.name + f"_{name}.tsv.gz"), sep="\t", index=True, index_label="barcode", compression="gzip")


def trajectory_mask(obs: pd.DataFrame) -> np.ndarray:
    excluded = ("cDC", "pDC", "_DC", "Cycling", "ambiguous", "contaminant", "Neutrophil", "low_quality")
    labels = obs["subtype"].astype(str)
    return ~labels.str.contains("|".join(excluded), case=False, regex=True).to_numpy()


def communication_mask(obs: pd.DataFrame) -> np.ndarray:
    excluded = ("cDC", "pDC", "_DC", "Cycling", "ambiguous", "contaminant", "Neutrophil", "low_quality")
    labels = obs["subtype"].astype(str)
    return ~labels.str.contains("|".join(excluded), case=False, regex=True).to_numpy()


def hsc_communication_mask(obs: pd.DataFrame) -> np.ndarray:
    labels = obs["subtype"].astype(str)
    return ~labels.str.contains("contaminant|low_quality|ambiguous", case=False, regex=True).to_numpy()


def export_dataset(dataset: str) -> None:
    myeloid = ad.read_h5ad(PROCESSED / f"{dataset}_myeloid_reclustered.h5ad")
    hsc = ad.read_h5ad(PROCESSED / f"{dataset}_hsc_reclustered.h5ad")

    trajectory = myeloid[trajectory_mask(myeloid.obs)].copy()
    trajectory_counts = trajectory.layers["counts"]
    embeddings = {}
    for key in ("X_pca_harmony", "X_pca", "X_umap"):
        if key in trajectory.obsm:
            embeddings[key.removeprefix("X_")] = trajectory.obsm[key]
    write_bundle(
        trajectory_counts,
        trajectory.obs,
        trajectory.var,
        PROCESSED / f"{dataset}_trajectory",
        embeddings,
    )

    macrophages = myeloid[communication_mask(myeloid.obs)].copy()
    hsc = hsc[hsc_communication_mask(hsc.obs)].copy()
    common_genes = macrophages.var_names.intersection(hsc.var_names)
    macrophages = macrophages[:, common_genes]
    hsc = hsc[:, common_genes]
    macrophage_counts = macrophages.layers["counts"]
    hsc_counts = hsc.layers["counts"]
    combined_counts = sparse.vstack([macrophage_counts, hsc_counts], format="csr")
    macrophage_obs = macrophages.obs.copy()
    hsc_obs = hsc.obs.copy()
    macrophage_obs["compartment"] = "Macrophage"
    hsc_obs["compartment"] = "HSC"
    macrophage_obs["cell_group"] = "Mac_" + macrophage_obs["subtype"].astype(str)
    hsc_obs["cell_group"] = "HSC_" + hsc_obs["subtype"].astype(str)
    combined_obs = pd.concat([macrophage_obs, hsc_obs], axis=0)
    write_bundle(
        combined_counts,
        combined_obs,
        macrophages.var.copy(),
        PROCESSED / f"{dataset}_communication",
    )

    summary = pd.DataFrame(
        {
            "dataset": [dataset],
            "trajectory_cells": [trajectory.n_obs],
            "communication_macrophages": [macrophages.n_obs],
            "communication_hsc": [hsc.n_obs],
            "communication_genes": [len(common_genes)],
        }
    )
    summary.to_csv(RESULTS / f"006_{dataset}_R_input_summary.tsv", sep="\t", index=False)


def main() -> None:
    datasets = sys.argv[1:] or ["GSE344087", "GSE298719", "GSE212837"]
    for dataset in datasets:
        export_dataset(dataset)


if __name__ == "__main__":
    main()
