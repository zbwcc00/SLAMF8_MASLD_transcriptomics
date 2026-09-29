from __future__ import annotations
import os

import gzip
import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import io, sparse


SEED = 20260920
PROJECT_ROOT = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))
SOURCE_ROOT = Path(os.environ.get("SLAMF8_MASLD_SOURCE_ROOT", Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1]))))
FARG_PATH = Path(os.environ.get("SLAMF8_MASLD_ROOT", Path(__file__).resolve().parents[1])) / "data" / "reference" / "FARG95_gene_set_provenance.tsv"
PROCESSED = PROJECT_ROOT / "data_processed"
RESULTS = PROJECT_ROOT / "results"
FIGURES = PROJECT_ROOT / "figures"

for directory in (PROCESSED, RESULTS, FIGURES):
    directory.mkdir(parents=True, exist_ok=True)

sc.settings.seed = SEED
sc.settings.verbosity = 2

MYELOID_LABELS = {
    "GSE344087": {"Monocyte", "Macrophage", "Dendritic_cDC", "Dendritic_pDC", "Immune_ambiguous"},
    "GSE298719": {"Macrophage", "Dendritic_cDC", "Neutrophil"},
    "GSE212837": {"Macrophage", "Dendritic_cDC"},
}
HSC_LABELS = {"Stellate_HSC", "Stellate_fibroblast"}

MYELOID_SIGNATURES = {
    "FCN1_VCAN_monocyte": ["FCN1", "VCAN", "S100A8", "S100A9", "LILRB1", "CTSS", "CCR2"],
    "Kupffer": ["C1QA", "C1QB", "C1QC", "MARCO", "VSIG4", "TIMD4", "FOLR2", "CD5L"],
    "TREM2_CD9_SAMac": ["TREM2", "CD9", "GPNMB", "SPP1", "LGALS3", "CTSL", "FABP5"],
    "OLR1_IL1B_SAMac": ["OLR1", "IL1B", "SPP1", "TREM2", "FCN1", "EREG"],
    "Inflammatory_macrophage": ["IL1B", "TNF", "NFKBIA", "CXCL8", "CCL3", "CCL4", "NLRP3"],
    "cDC2": ["FCER1A", "CD1C", "CLEC10A", "CST3", "HLA-DQA1"],
    "cDC1": ["CLEC9A", "BATF3", "XCR1", "CADM1"],
    "pDC": ["GZMB", "JCHAIN", "TCF4", "IL3RA", "CLEC4C"],
    "Cycling_myeloid": ["MKI67", "TOP2A", "TYMS", "UBE2C"],
}

HSC_SIGNATURES = {
    "Quiescent_HSC": ["RBP1", "LRAT", "CYGB", "RELN", "PDGFRB", "DES", "COLEC11"],
    "Activated_HSC": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "LOX", "PDGFRA"],
    "Myofibroblast": ["ACTA2", "TAGLN", "MYL9", "CNN1", "TPM2", "CALD1"],
    "Portal_fibroblast": ["COL15A1", "PI16", "DPT", "MMP2", "COL14A1", "FBLN1"],
}


def load_gene_sets() -> dict[str, list[str]]:
    farg_table = pd.read_csv(FARG_PATH, sep="\t")
    gene_column = "gene" if "gene" in farg_table.columns else farg_table.columns[0]
    fargs = farg_table[gene_column].dropna().astype(str).str.upper().drop_duplicates().tolist()
    return {
        "FARG95": fargs,
        "Iron_homeostasis": ["TFRC", "SLC11A2", "FTH1", "FTL", "SLC40A1", "NCOA4", "HMOX1", "STEAP3", "PCBP1", "PCBP2", "HAMP"],
        "PUFA_ACSL4": ["ACSL4", "LPCAT3", "AGPAT3", "ELOVL5", "FADS1", "FADS2", "ALOX15", "PLA2G6"],
        "Lipid_peroxidation": ["ACSL4", "ALOX5", "ALOX12", "ALOX15", "POR", "NOX1", "NOX2", "CYBB", "PTGS2"],
        "Antioxidant_defense": ["GPX4", "SLC7A11", "GCLC", "GCLM", "NFE2L2", "FSP1", "AIFM2", "DHODH", "GCH1"],
        "Senescence_SASP": ["CDKN1A", "CDKN2A", "GLB1", "SERPINE1", "GDF15", "IL6", "IL1B", "CXCL8", "MMP3", "MMP9", "CCL2", "TGFB1", "IGFBP7"],
        "Inflammation": ["IL1B", "TNF", "NFKBIA", "CXCL8", "CCL2", "CCL3", "CCL4", "NLRP3", "STAT1", "IRF1"],
        "HSC_activation": ["YAP1", "NUAK2", "COL1A1", "COL1A2", "COL3A1", "LOX", "ACTA2", "TAGLN", "TGFB1", "PDGFRA", "PDGFRB"],
    }


GENE_SETS = load_gene_sets()


def make_unique(values: pd.Series | pd.Index) -> pd.Index:
    return ad.utils.make_index_unique(pd.Index(values.astype(str)))


def read_triplet(gsm: str, raw_dir: Path) -> tuple[sparse.csr_matrix, pd.DataFrame, np.ndarray]:
    matrix_path = next(raw_dir.glob(f"{gsm}_*_matrix.mtx.gz"))
    feature_path = next(raw_dir.glob(f"{gsm}_*_features.tsv.gz"))
    barcode_path = next(raw_dir.glob(f"{gsm}_*_barcodes.tsv.gz"))
    with gzip.open(matrix_path, "rb") as handle:
        matrix = io.mmread(handle).tocsr().T.tocsr()
    features = pd.read_csv(feature_path, sep="\t", header=None, compression="gzip")
    barcodes = pd.read_csv(barcode_path, sep="\t", header=None, compression="gzip")[0].astype(str).to_numpy()
    var = pd.DataFrame({"gene_id": features[0].astype(str), "gene": features[1].astype(str)})
    var.index = make_unique(var["gene"])
    return matrix, var, barcodes


def load_gse212837() -> ad.AnnData:
    metadata = pd.read_csv(
        SOURCE_ROOT / "results" / "032_GSE212837_cell_metadata.tsv.gz", sep="\t"
    )
    wanted_labels = MYELOID_LABELS["GSE212837"] | HSC_LABELS
    metadata = metadata.loc[
        metadata["passes_qc"].astype(bool) & metadata["cell_type"].isin(wanted_labels)
    ].copy()
    raw_dir = SOURCE_ROOT / "data_raw" / "GSE212837"
    matrices: list[sparse.csr_matrix] = []
    observations: list[pd.DataFrame] = []
    reference_var: pd.DataFrame | None = None
    for gsm, group in metadata.groupby("gsm", sort=True):
        matrix, var, barcodes = read_triplet(gsm, raw_dir)
        if reference_var is None:
            reference_var = var
        elif not np.array_equal(reference_var["gene_id"].to_numpy(), var["gene_id"].to_numpy()):
            raise ValueError(f"Gene order differs in {gsm}")
        barcode_to_position = pd.Series(np.arange(len(barcodes)), index=barcodes)
        positions = barcode_to_position.reindex(group["original_barcode"]).to_numpy()
        if np.isnan(positions).any():
            raise ValueError(f"Missing selected barcodes in {gsm}")
        matrices.append(matrix[positions.astype(int), :])
        observations.append(group.copy())
    if reference_var is None:
        raise RuntimeError("No GSE212837 cells selected")
    obs = pd.concat(observations, ignore_index=True).set_index("barcode", drop=False)
    return ad.AnnData(X=sparse.vstack(matrices, format="csr"), obs=obs, var=reference_var)


def load_source(dataset: str) -> ad.AnnData:
    if dataset == "GSE212837":
        return load_gse212837()
    source_path = SOURCE_ROOT / "data_processed" / f"{dataset}_combined_counts_qc_annotations.h5ad"
    backed = ad.read_h5ad(source_path, backed="r")
    wanted = MYELOID_LABELS[dataset] | HSC_LABELS
    mask = backed.obs["passes_qc"].astype(bool) & backed.obs["cell_type"].isin(wanted)
    subset = backed[mask.to_numpy(), :].to_memory()
    backed.file.close()
    if "gene" in subset.var.columns:
        subset.var_names = make_unique(subset.var["gene"])
    return subset


def attach_doublets(adata: ad.AnnData, dataset: str) -> ad.AnnData:
    path = RESULTS / f"001_{dataset}_scDblFinder_cell_calls.tsv.gz"
    calls = pd.read_csv(path, sep="\t").set_index("barcode")
    if calls.index.has_duplicates:
        raise ValueError(f"Duplicate doublet calls in {dataset}")
    joined = calls.reindex(adata.obs_names)
    if joined["scDblFinder_class"].isna().any():
        missing = int(joined["scDblFinder_class"].isna().sum())
        raise ValueError(f"{dataset}: {missing} selected cells lack scDblFinder calls")
    adata.obs["scDblFinder_score"] = joined["scDblFinder_score"].to_numpy()
    adata.obs["scDblFinder_class"] = joined["scDblFinder_class"].astype(str).to_numpy()
    return adata[adata.obs["scDblFinder_class"].eq("singlet")].copy()


def add_signature_scores(adata: ad.AnnData, signatures: dict[str, list[str]]) -> pd.DataFrame:
    unique_genes = sorted({gene for genes in signatures.values() for gene in genes})
    present = [gene for gene in unique_genes if gene in adata.var_names]
    matrix = adata[:, present].X
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix, dtype=np.float32)
    means = matrix.mean(axis=0)
    standard_deviations = matrix.std(axis=0)
    standard_deviations[standard_deviations == 0] = 1
    standardized = (matrix - means) / standard_deviations
    gene_to_column = {gene: index for index, gene in enumerate(present)}
    coverage_rows = []
    for name, genes in signatures.items():
        available = [gene for gene in genes if gene in gene_to_column]
        coverage_rows.append({"signature": name, "n_defined": len(genes), "n_present": len(available)})
        if available:
            positions = [gene_to_column[gene] for gene in available]
            adata.obs[f"score_{name}"] = standardized[:, positions].mean(axis=1)
        else:
            adata.obs[f"score_{name}"] = np.nan
    return pd.DataFrame(coverage_rows)


def annotate_clusters(adata: ad.AnnData, signatures: dict[str, list[str]], lineage: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    coverage = add_signature_scores(adata, signatures)
    score_columns = [f"score_{name}" for name in signatures]
    cluster_scores = adata.obs.groupby("leiden", observed=True)[score_columns].mean()
    cluster_scores.columns = [column.removeprefix("score_") for column in cluster_scores.columns]
    assignments = {}
    for cluster, row in cluster_scores.iterrows():
        ordered = row.sort_values(ascending=False)
        best = str(ordered.index[0])
        margin = float(ordered.iloc[0] - ordered.iloc[1]) if len(ordered) > 1 else np.inf
        if not np.isfinite(ordered.iloc[0]) or (ordered.iloc[0] < 0 and margin < 0.10):
            label = f"{lineage}_ambiguous"
        elif margin < 0.05:
            label = f"{best}_mixed"
        else:
            label = best
        assignments[str(cluster)] = label
    adata.obs["subtype"] = adata.obs["leiden"].astype(str).map(assignments).astype("category")
    score_table = cluster_scores.reset_index()
    score_table["assigned_subtype"] = score_table["leiden"].astype(str).map(assignments)
    return score_table, coverage


def add_biology_modules(adata: ad.AnnData) -> pd.DataFrame:
    rows = []
    for name, genes in GENE_SETS.items():
        available = [gene for gene in genes if gene in adata.var_names]
        rows.append({"module": name, "n_defined": len(genes), "n_present": len(available)})
        if not available:
            adata.obs[f"module_{name}"] = np.nan
            continue
        values = adata[:, available].X
        if sparse.issparse(values):
            values = np.asarray(values.mean(axis=1)).ravel()
        else:
            values = np.asarray(values).mean(axis=1)
        adata.obs[f"module_{name}"] = values
    return pd.DataFrame(rows)


def run_reclustering(adata: ad.AnnData, dataset: str, lineage: str) -> ad.AnnData:
    adata.var_names_make_unique()
    adata.var.index.name = "var_name"
    sc.pp.filter_genes(adata, min_cells=3)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    add_biology_modules(adata)
    batch_key = "gsm"
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=min(3000, adata.n_vars), batch_key=batch_key)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=min(50, adata.n_obs - 1, int(adata.var["highly_variable"].sum()) - 1), use_highly_variable=True, random_state=SEED)
    representation = "X_pca"
    if adata.obs[batch_key].nunique() > 1:
        try:
            sc.external.pp.harmony_integrate(adata, key=batch_key, basis="X_pca", adjusted_basis="X_pca_harmony", random_state=SEED)
            representation = "X_pca_harmony"
        except Exception as error:
            adata.uns["harmony_error"] = str(error)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=min(30, adata.obsm[representation].shape[1]), use_rep=representation, random_state=SEED)
    sc.tl.umap(adata, min_dist=0.35, random_state=SEED)
    resolution = 0.65 if lineage == "myeloid" else 0.45
    sc.tl.leiden(adata, resolution=resolution, key_added="leiden", random_state=SEED)
    signatures = MYELOID_SIGNATURES if lineage == "myeloid" else HSC_SIGNATURES
    cluster_scores, signature_coverage = annotate_clusters(adata, signatures, lineage)
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", pts=True, use_raw=True)
    markers = sc.get.rank_genes_groups_df(adata, group=None)

    prefix = "003" if lineage == "myeloid" else "004"
    cluster_scores.to_csv(RESULTS / f"{prefix}_{dataset}_{lineage}_cluster_signature_scores.tsv", sep="\t", index=False)
    signature_coverage.to_csv(RESULTS / f"{prefix}_{dataset}_{lineage}_signature_coverage.tsv", sep="\t", index=False)
    markers.to_csv(RESULTS / f"{prefix}_{dataset}_{lineage}_cluster_markers.tsv.gz", sep="\t", index=False, compression="gzip")
    count_columns = [column for column in ["gsm", "sample_id", "patient_id", "disease", "fibrosis", "subtype"] if column in adata.obs]
    subtype_counts = adata.obs.groupby(count_columns, observed=True).size().rename("n_cells").reset_index()
    subtype_counts.to_csv(RESULTS / f"{prefix}_{dataset}_{lineage}_sample_subtype_counts.tsv", sep="\t", index=False)

    output_path = PROCESSED / f"{dataset}_{lineage}_reclustered.h5ad"
    adata.write_h5ad(output_path, compression="gzip")
    save_umap(adata, dataset, lineage)
    return adata


def save_umap(adata: ad.AnnData, dataset: str, lineage: str) -> None:
    sns.set_theme(style="white", context="notebook")
    coordinates = adata.obsm["X_umap"]
    columns = ["subtype", "disease"]
    if "fibrosis" in adata.obs:
        columns.append("fibrosis")
    color_values = adata.raw[:, "SLAMF8"].X if lineage == "myeloid" and "SLAMF8" in adata.raw.var_names else None
    n_panels = len(columns) + (1 if color_values is not None else 0)
    fig, axes = plt.subplots(1, n_panels, figsize=(6 * n_panels, 5))
    if n_panels == 1:
        axes = [axes]
    for axis, column in zip(axes, columns):
        categories = adata.obs[column].astype(str)
        for category in sorted(categories.unique()):
            mask = categories.eq(category).to_numpy()
            axis.scatter(coordinates[mask, 0], coordinates[mask, 1], s=4, alpha=0.55, label=category, rasterized=True)
        axis.set_title(f"{dataset} {lineage}: {column}")
        axis.set_xticks([])
        axis.set_yticks([])
        axis.legend(markerscale=3, frameon=False, fontsize=8)
    if color_values is not None:
        if sparse.issparse(color_values):
            color_values = color_values.toarray().ravel()
        else:
            color_values = np.asarray(color_values).ravel()
        plot = axes[-1].scatter(coordinates[:, 0], coordinates[:, 1], c=color_values, cmap="magma", s=4, rasterized=True)
        axes[-1].set_title(f"{dataset} {lineage}: SLAMF8")
        axes[-1].set_xticks([])
        axes[-1].set_yticks([])
        fig.colorbar(plot, ax=axes[-1], fraction=0.046)
    fig.tight_layout()
    figure_number = {"GSE344087": "001", "GSE298719": "002", "GSE212837": "003"}[dataset]
    suffix = "A" if lineage == "myeloid" else "B"
    fig.savefig(FIGURES / f"{figure_number}{suffix}_{dataset}_{lineage}_UMAP.png", dpi=240)
    plt.close(fig)


def main() -> None:
    requested = sys.argv[1:]
    datasets = requested if requested else list(MYELOID_LABELS)
    unknown = sorted(set(datasets) - set(MYELOID_LABELS))
    if unknown:
        raise ValueError(f"Unknown datasets: {unknown}")
    manifest_path = RESULTS / "005_reclustering_manifest.json"
    if manifest_path.exists():
        with open(manifest_path, encoding="utf-8") as handle:
            manifest = json.load(handle)
    else:
        manifest = {"seed": SEED, "datasets": {}}
    for dataset in datasets:
        print(f"Loading {dataset}", flush=True)
        relevant = attach_doublets(load_source(dataset), dataset)
        myeloid = relevant[relevant.obs["cell_type"].isin(MYELOID_LABELS[dataset])].copy()
        hsc = relevant[relevant.obs["cell_type"].isin(HSC_LABELS)].copy()
        print(f"{dataset}: {myeloid.n_obs} myeloid singlets, {hsc.n_obs} HSC/fibroblast singlets", flush=True)
        myeloid = run_reclustering(myeloid, dataset, "myeloid")
        hsc = run_reclustering(hsc, dataset, "hsc")
        manifest["datasets"][dataset] = {
            "myeloid_singlets": int(myeloid.n_obs),
            "hsc_singlets": int(hsc.n_obs),
            "myeloid_subtypes": myeloid.obs["subtype"].value_counts().to_dict(),
            "hsc_subtypes": hsc.obs["subtype"].value_counts().to_dict(),
            "harmony_myeloid": "X_pca_harmony" in myeloid.obsm,
            "harmony_hsc": "X_pca_harmony" in hsc.obsm,
        }
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    main()
