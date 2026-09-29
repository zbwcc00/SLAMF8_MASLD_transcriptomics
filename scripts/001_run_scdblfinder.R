#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: 001_run_scdblfinder.R <GSE344087|GSE298719|GSE212837> [limit]")
}

dataset <- args[[1L]]
limit <- if (length(args) >= 2L) as.integer(args[[2L]]) else NA_integer_
allowed <- c("GSE344087", "GSE298719", "GSE212837")
if (!dataset %in% allowed) {
  stop("Unknown dataset: ", dataset)
}

seed <- 20260920L
project_root <- normalizePath(Sys.getenv("SLAMF8_MASLD_ROOT", "."), winslash = "/", mustWork = FALSE)
source_root <- normalizePath(Sys.getenv("SLAMF8_MASLD_SOURCE_ROOT", project_root), winslash = "/", mustWork = FALSE)
raw_dir <- file.path(source_root, "data_raw", dataset)
source_results <- file.path(source_root, "results")
output_dir <- file.path(project_root, "results", "doublet", dataset)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata_files <- c(
  GSE344087 = "005_GSE344087_cell_metadata.tsv.gz",
  GSE298719 = "018_GSE298719_cell_metadata.tsv.gz",
  GSE212837 = "032_GSE212837_cell_metadata.tsv.gz"
)
cell_metadata <- read.delim(
  file.path(source_results, metadata_files[[dataset]]),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
cell_metadata <- cell_metadata[cell_metadata$passes_qc %in% c(TRUE, "True", "TRUE"), , drop = FALSE]
gsms <- sort(unique(cell_metadata$gsm))
if (!is.na(limit)) {
  gsms <- head(gsms, limit)
}

read_triplet <- function(gsm) {
  matrix_path <- list.files(
    raw_dir,
    pattern = paste0("^", gsm, ".*_matrix\\.mtx\\.gz$"),
    full.names = TRUE
  )
  feature_path <- list.files(
    raw_dir,
    pattern = paste0("^", gsm, ".*_features\\.tsv\\.gz$"),
    full.names = TRUE
  )
  barcode_path <- list.files(
    raw_dir,
    pattern = paste0("^", gsm, ".*_barcodes\\.tsv\\.gz$"),
    full.names = TRUE
  )
  if (length(matrix_path) != 1L || length(feature_path) != 1L || length(barcode_path) != 1L) {
    stop("Incomplete or ambiguous 10x triplet for ", gsm)
  }
  counts <- readMM(gzfile(matrix_path))
  features <- read.delim(gzfile(feature_path), header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(gzfile(barcode_path), header = FALSE, stringsAsFactors = FALSE)[[1L]]
  rownames(counts) <- make.unique(as.character(features[[2L]]))
  colnames(counts) <- as.character(barcodes)
  counts
}

read_library <- function(gsm) {
  if (dataset == "GSE344087") {
    suppressPackageStartupMessages(library(Seurat))
    h5_path <- list.files(
      raw_dir,
      pattern = paste0("^", gsm, ".*filtered_feature_bc_matrix\\.h5$"),
      full.names = TRUE
    )
    if (length(h5_path) != 1L) {
      stop("Missing or ambiguous 10x H5 for ", gsm)
    }
    counts <- Read10X_h5(h5_path, use.names = TRUE, unique.features = TRUE)
    if (is.list(counts)) {
      if (!"Gene Expression" %in% names(counts)) {
        stop("Gene Expression assay absent in ", gsm)
      }
      counts <- counts[["Gene Expression"]]
    }
    return(as(counts, "dgCMatrix"))
  }
  as(read_triplet(gsm), "dgCMatrix")
}

run_one <- function(gsm) {
  output_path <- file.path(output_dir, paste0(gsm, "_scDblFinder.rds"))
  if (file.exists(output_path)) {
    message("Reusing ", output_path)
    return(readRDS(output_path))
  }
  message("Reading ", dataset, " / ", gsm)
  counts <- read_library(gsm)
  wanted <- cell_metadata$original_barcode[cell_metadata$gsm == gsm]
  keep <- match(wanted, colnames(counts), nomatch = 0L)
  missing <- sum(keep == 0L)
  if (missing > 0L) {
    stop(gsm, ": ", missing, " QC-passing barcodes absent from raw matrix")
  }
  counts <- counts[, keep, drop = FALSE]
  if (ncol(counts) < 50L) {
    stop(gsm, ": fewer than 50 QC-passing cells")
  }
  set.seed(seed + match(gsm, sort(unique(cell_metadata$gsm))))
  sce <- SingleCellExperiment(list(counts = counts))
  sce <- scDblFinder(sce, clusters = TRUE, verbose = TRUE)
  result <- data.frame(
    dataset = dataset,
    gsm = gsm,
    original_barcode = colnames(sce),
    barcode = paste0(gsm, "_", colnames(sce)),
    scDblFinder_score = colData(sce)$scDblFinder.score,
    scDblFinder_class = as.character(colData(sce)$scDblFinder.class),
    stringsAsFactors = FALSE
  )
  saveRDS(result, output_path)
  message(
    "Finished ", gsm, ": ", nrow(result), " cells; ",
    sum(result$scDblFinder_class == "doublet"), " predicted doublets"
  )
  result
}

results <- lapply(gsms, run_one)
combined <- do.call(rbind, results)
combined_path <- file.path(project_root, "results", paste0("001_", dataset, "_scDblFinder_cell_calls.tsv.gz"))
write.table(
  combined,
  gzfile(combined_path),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

summary_table <- aggregate(
  barcode ~ dataset + gsm + scDblFinder_class,
  data = combined,
  FUN = length
)
names(summary_table)[names(summary_table) == "barcode"] <- "n_cells"
totals <- aggregate(n_cells ~ dataset + gsm, data = summary_table, FUN = sum)
names(totals)[names(totals) == "n_cells"] <- "n_total"
summary_table <- merge(summary_table, totals, by = c("dataset", "gsm"))
summary_table$fraction <- summary_table$n_cells / summary_table$n_total
write.table(
  summary_table,
  file.path(project_root, "results", paste0("002_", dataset, "_scDblFinder_summary.tsv")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Wrote ", combined_path)
