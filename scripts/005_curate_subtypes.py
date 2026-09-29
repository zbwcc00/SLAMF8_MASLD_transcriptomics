from __future__ import annotations
import os

import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import scanpy as sc


PROJECT_ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
PROCESSED = PROJECT_ROOT / "data_processed"
RESULTS = PROJECT_ROOT / "results"
FIGURES = PROJECT_ROOT / "figures"

CURATED = {
    ("GSE344087", "myeloid"): {
        "0": "FCN1_VCAN_monocyte",
        "1": "T_NK_contaminant",
        "2": "cDC2",
        "3": "TREM2_CD9_SAMac",
        "4": "FCGR3A_monocyte",
        "5": "T_NK_contaminant",
        "6": "pDC",
        "7": "Kupffer",
        "8": "T_NK_contaminant",
        "9": "cDC1",
        "10": "T_NK_contaminant",
        "11": "cDC2_low_quality",
        "12": "Kupffer",
        "13": "Cycling_myeloid",
        "14": "LAMP3_DC",
        "15": "Myeloid_ambiguous",
        "16": "Endothelial_contaminant",
    },
    ("GSE298719", "myeloid"): {
        "0": "Kupffer",
        "1": "Resident_macrophage",
        "2": "FCN1_VCAN_monocyte",
        "3": "cDC2",
        "4": "STAB1_FOLR2_macrophage",
        "5": "GPNMB_TREM2_SAMac",
        "6": "Resident_macrophage",
        "7": "Resident_macrophage",
        "8": "MARCO_HMOX1_Kupffer",
        "9": "Hepatocyte_contaminant",
        "10": "TREM2_CD9_SAMac",
        "11": "Endothelial_contaminant",
        "12": "Cycling_myeloid",
        "13": "cDC1",
        "14": "Neutrophil",
        "15": "T_NK_contaminant",
        "16": "Fibroblast_contaminant",
        "17": "Inflammatory_macrophage",
    },
    ("GSE344087", "hsc"): {
        "0": "Myofibroblast",
        "1": "Activated_HSC",
        "2": "Quiescent_HSC",
        "3": "HSC_low_quality",
        "4": "Activated_HSC_stress",
        "5": "Quiescent_HSC",
        "6": "Quiescent_HSC",
    },
    ("GSE298719", "hsc"): {
        "0": "Portal_fibroblast",
        "1": "Quiescent_HSC",
        "2": "Hepatocyte_contaminant",
        "3": "Portal_fibroblast",
        "4": "Quiescent_HSC",
        "5": "Myofibroblast",
        "6": "Activated_HSC",
        "7": "Macrophage_contaminant",
        "8": "Myofibroblast",
        "9": "Endothelial_contaminant",
        "10": "Portal_fibroblast",
    },
    ("GSE212837", "myeloid"): {
        "0": "Resident_macrophage",
        "1": "Resident_macrophage",
        "2": "Kupffer",
        "3": "Kupffer",
        "4": "Activated_Kupffer",
        "5": "TREM2_CD9_SPP1_SAMac",
        "6": "FCN1_VCAN_monocyte",
        "7": "Mixed_T_macrophage_low_confidence",
        "8": "SPP1_GPNMB_macrophage",
        "9": "SPP1_FCN1_transition",
        "10": "Endothelial_like_macrophage_low_confidence",
        "11": "Cycling_myeloid",
        "12": "cDC1",
        "13": "HSC_like_macrophage_low_confidence",
    },
    ("GSE212837", "hsc"): {
        "0": "Quiescent_HSC",
        "1": "Activated_HSC_acute_phase",
        "2": "Hepatocyte_contaminant",
        "3": "Hepatocyte_mito_contaminant",
        "4": "Myofibroblast",
        "5": "Activated_HSC",
        "6": "Quiescent_HSC",
        "7": "Activated_HSC_acute_phase",
        "8": "Hepatocyte_contaminant",
        "9": "T_cell_contaminant",
        "10": "Endothelial_contaminant",
        "11": "Macrophage_contaminant",
    },
}


def curate(dataset: str, lineage: str) -> None:
    mapping = CURATED.get((dataset, lineage))
    if mapping is None:
        print(f"No curated map for {dataset} {lineage}; retaining automated labels")
        return
    path = PROCESSED / f"{dataset}_{lineage}_reclustered.h5ad"
    adata = ad.read_h5ad(path)
    clusters = set(adata.obs["leiden"].astype(str).unique())
    if clusters != set(mapping):
        raise ValueError(f"Cluster map mismatch for {dataset} {lineage}: observed={clusters}, mapped={set(mapping)}")
    adata.obs["subtype_auto"] = adata.obs["subtype"].astype(str)
    adata.obs["subtype"] = adata.obs["leiden"].astype(str).map(mapping).astype("category")
    adata.uns["subtype_annotation"] = {
        "method": "curated marker-panel review",
        "mapping": mapping,
    }
    adata.write_h5ad(path, compression="gzip")

    prefix = "003" if lineage == "myeloid" else "004"
    score_path = RESULTS / f"{prefix}_{dataset}_{lineage}_cluster_signature_scores.tsv"
    scores = __import__("pandas").read_csv(score_path, sep="\t")
    scores["assigned_subtype_auto"] = scores["assigned_subtype"]
    scores["assigned_subtype"] = scores["leiden"].astype(str).map(mapping)
    scores.to_csv(score_path, sep="\t", index=False)

    count_columns = [column for column in ["gsm", "sample_id", "patient_id", "disease", "fibrosis", "subtype"] if column in adata.obs]
    subtype_counts = adata.obs.groupby(count_columns, observed=True).size().rename("n_cells").reset_index()
    subtype_counts.to_csv(RESULTS / f"{prefix}_{dataset}_{lineage}_sample_subtype_counts.tsv", sep="\t", index=False)

    colors = ["subtype", "disease"] + (["fibrosis"] if "fibrosis" in adata.obs else [])
    sc.pl.umap(adata, color=colors, ncols=len(colors), show=False, frameon=False)
    plt.gcf().set_size_inches(6 * len(colors), 5)
    plt.tight_layout()
    suffix = "A" if lineage == "myeloid" else "B"
    figure_number = {"GSE344087": "001", "GSE298719": "002", "GSE212837": "003"}[dataset]
    plt.savefig(FIGURES / f"{figure_number}{suffix}_{dataset}_{lineage}_UMAP_curated.png", dpi=240)
    plt.close("all")


def main() -> None:
    requested = sys.argv[1:] or ["GSE344087", "GSE298719"]
    for dataset in requested:
        for lineage in ("myeloid", "hsc"):
            curate(dataset, lineage)
    with open(RESULTS / "012_curated_annotation_map.json", "w", encoding="utf-8") as handle:
        json.dump({f"{dataset}_{lineage}": mapping for (dataset, lineage), mapping in CURATED.items()}, handle, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    main()
