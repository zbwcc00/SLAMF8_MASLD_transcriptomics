#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(scuttle)
  library(liana)
  library(CellChat)
})

args <- commandArgs(trailingOnly = TRUE)
datasets <- if (length(args) >= 1L) args[[1L]] else c("GSE344087", "GSE298719", "GSE212837")
requested_condition <- if (length(args) >= 2L) args[[2L]] else NULL
seed <- 20260920L
project_root <- normalizePath(Sys.getenv("SLAMF8_MASLD_ROOT", "."), winslash = "/", mustWork = FALSE)
processed_dir <- file.path(project_root, "data_processed")
results_dir <- file.path(project_root, "results")
figures_dir <- file.path(project_root, "figures")
dir.create(file.path(results_dir, "communication"), recursive = TRUE, showWarnings = FALSE)

condition_columns <- c(
  GSE344087 = "fibrosis",
  GSE298719 = "disease",
  GSE212837 = "disease"
)

read_bundle <- function(dataset) {
  prefix <- file.path(processed_dir, paste0(dataset, "_communication"))
  counts <- as(readMM(gzfile(paste0(prefix, "_counts.mtx.gz"))), "dgCMatrix")
  features <- read.delim(gzfile(paste0(prefix, "_features.tsv.gz")), check.names = FALSE)
  metadata <- read.delim(gzfile(paste0(prefix, "_cell_metadata.tsv.gz")), check.names = FALSE)
  rownames(counts) <- make.unique(as.character(features$var_name))
  colnames(counts) <- metadata$barcode
  rownames(metadata) <- metadata$barcode
  list(counts = counts, metadata = metadata)
}

downsample_groups <- function(counts, metadata, max_cells = 1000L) {
  set.seed(seed)
  split_cells <- split(seq_len(nrow(metadata)), metadata$cell_group)
  selected <- unlist(lapply(split_cells, function(index) sample(index, min(length(index), max_cells))), use.names = FALSE)
  selected <- sort(selected)
  list(counts = counts[, selected, drop = FALSE], metadata = metadata[selected, , drop = FALSE])
}

write_tsv <- function(data, path) {
  write.table(data, path, sep = "\t", row.names = FALSE, quote = FALSE)
}

focus_pathway <- function(ligand, receptor) {
  pair <- toupper(paste(ligand, receptor, sep = "_"))
  if (grepl("TGFB", pair)) return("TGFb")
  if (grepl("PDGF", pair)) return("PDGF")
  if (grepl("TNF", pair)) return("TNF")
  if (grepl("IL1", pair)) return("IL1")
  if (grepl("JAG|DLL|NOTCH", pair)) return("NOTCH")
  if (grepl("SPP1", pair)) return("SPP1")
  "Other"
}

run_liana <- function(dataset, condition, counts, metadata, methods) {
  sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = metadata
  )
  sce <- logNormCounts(sce)
  result <- liana_wrap(
    sce,
    method = methods,
    resource = "Consensus",
    idents_col = "cell_group",
    min_cells = 10,
    return_all = TRUE,
    assay = "logcounts",
    verbose = TRUE
  )
  aggregated <- liana_aggregate(
    result,
    aggregate_how = "magnitude",
    get_ranks = TRUE,
    get_agrank = TRUE,
    verbose = TRUE
  )
  saveRDS(result, file.path(processed_dir, paste0(dataset, "_", condition, "_LIANA_raw.rds")))
  saveRDS(aggregated, file.path(processed_dir, paste0(dataset, "_", condition, "_LIANA_consensus.rds")))
  write_tsv(aggregated, file.path(results_dir, "communication", paste0("017_", dataset, "_", condition, "_LIANA_consensus.tsv")))

  macrophage_hsc <- aggregated[
    grepl("^Mac_", aggregated$source) & grepl("^HSC_", aggregated$target),
    , drop = FALSE
  ]
  ligand_col <- intersect(c("ligand_complex", "ligand.complex"), names(macrophage_hsc))[[1L]]
  receptor_col <- intersect(c("receptor_complex", "receptor.complex"), names(macrophage_hsc))[[1L]]
  macrophage_hsc$focus_pathway <- mapply(focus_pathway, macrophage_hsc[[ligand_col]], macrophage_hsc[[receptor_col]])
  write_tsv(macrophage_hsc, file.path(results_dir, "communication", paste0("018_", dataset, "_", condition, "_LIANA_macrophage_to_HSC.tsv")))
  macrophage_hsc
}

run_cellchat <- function(dataset, condition, counts, metadata) {
  library_sizes <- Matrix::colSums(counts)
  normalized <- t(t(counts) / pmax(library_sizes, 1)) * 1e4
  normalized@x <- log1p(normalized@x)
  cellchat <- createCellChat(object = normalized, meta = metadata, group.by = "cell_group")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  future::plan("sequential")
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE, population.size = FALSE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  saveRDS(cellchat, file.path(processed_dir, paste0(dataset, "_", condition, "_CellChat.rds")))
  sources <- grep("^Mac_", levels(cellchat@idents), value = TRUE)
  targets <- grep("^HSC_", levels(cellchat@idents), value = TRUE)
  interactions <- subsetCommunication(cellchat, sources.use = sources, targets.use = targets)
  write_tsv(interactions, file.path(results_dir, "communication", paste0("019_", dataset, "_", condition, "_CellChat_macrophage_to_HSC.tsv")))
  interactions
}

for (dataset in datasets) {
  message("Communication analysis: ", dataset)
  bundle <- read_bundle(dataset)
  condition_col <- condition_columns[[dataset]]
  conditions <- unique(as.character(bundle$metadata[[condition_col]]))
  if (!is.null(requested_condition)) {
    conditions <- conditions[tolower(conditions) == tolower(requested_condition)]
    if (!length(conditions)) stop("Requested condition not found: ", requested_condition)
  }
  for (condition in conditions) {
    message("  Condition: ", condition)
    keep <- which(as.character(bundle$metadata[[condition_col]]) == condition)
    max_cells <- if (!is.null(requested_condition)) 700L else 1000L
    liana_methods <- if (!is.null(requested_condition)) c("natmi", "cellphonedb", "sca") else c("natmi", "connectome", "cellphonedb", "sca")
    local <- downsample_groups(bundle$counts[, keep, drop = FALSE], bundle$metadata[keep, , drop = FALSE], max_cells = max_cells)
    safe_condition <- gsub("[^A-Za-z0-9]+", "_", condition)
    error_path <- file.path(results_dir, "communication", paste0("020_", dataset, "_", safe_condition, "_errors.txt"))
    errors <- character()
    tryCatch(
      run_cellchat(dataset, safe_condition, local$counts, local$metadata),
      error = function(error) errors <<- c(errors, paste("CellChat:", conditionMessage(error)))
    )
    tryCatch(
      run_liana(dataset, safe_condition, local$counts, local$metadata, methods = liana_methods),
      error = function(error) errors <<- c(errors, paste("LIANA:", conditionMessage(error)))
    )
    if (length(errors)) {
      writeLines(errors, error_path)
      warning(paste(errors, collapse = " | "))
    } else if (file.exists(error_path)) {
      unlink(error_path)
    }
  }
}
