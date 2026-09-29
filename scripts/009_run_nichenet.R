#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(nichenetr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
datasets <- if (length(args)) args else c("GSE344087", "GSE298719", "GSE212837")
project_root <- normalizePath(Sys.getenv("SLAMF8_MASLD_ROOT", "."), winslash = "/", mustWork = FALSE)
processed_dir <- file.path(project_root, "data_processed")
results_dir <- file.path(project_root, "results", "communication")
figures_dir <- file.path(project_root, "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

comparisons <- list(
  GSE344087 = list(column = "fibrosis", case = "fibrosis", control = "no fibrosis"),
  GSE298719 = list(column = "disease", case = "MASLD", control = "Healthy"),
  GSE212837 = list(column = "disease", case = "NASH", control = "Control")
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

sample_column <- function(metadata) {
  if ("patient_id" %in% names(metadata) && any(nzchar(as.character(metadata$patient_id)))) "patient_id" else "gsm"
}

aggregate_counts <- function(counts, metadata, unit_col) {
  groups <- factor(metadata[[unit_col]], levels = unique(metadata[[unit_col]]))
  design <- sparseMatrix(
    i = seq_along(groups),
    j = as.integer(groups),
    x = 1,
    dims = c(length(groups), nlevels(groups))
  )
  output <- counts %*% design
  colnames(output) <- levels(groups)
  output
}

run_edger <- function(counts, sample_metadata, condition_col, case_label, control_label) {
  condition <- factor(sample_metadata[[condition_col]], levels = c(control_label, case_label))
  keep_samples <- !is.na(condition)
  counts <- counts[, keep_samples, drop = FALSE]
  condition <- droplevels(condition[keep_samples])
  if (nlevels(condition) != 2L || min(table(condition)) < 2L) {
    stop("Insufficient patient/sample replication for edgeR")
  }
  design <- model.matrix(~ condition)
  dge <- DGEList(counts = counts, group = condition)
  keep_genes <- filterByExpr(dge, design = design)
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design, robust = TRUE)
  fit <- glmQLFit(dge, design, robust = TRUE)
  test <- glmQLFTest(fit, coef = ncol(design))
  table <- topTags(test, n = Inf, sort.by = "none")$table
  table$gene <- rownames(table)
  table
}

expressed_genes <- function(counts, minimum_fraction = 0.10) {
  rownames(counts)[Matrix::rowMeans(counts > 0) >= minimum_fraction]
}

resource_dir <- file.path(processed_dir, "nichenet_resources")
ligand_target_matrix <- readRDS(file.path(resource_dir, "ligand_target_matrix.rds"))
lr_network <- readRDS(file.path(resource_dir, "lr_network.rds"))

for (dataset in datasets) {
  message("NicheNet: ", dataset)
  bundle <- read_bundle(dataset)
  metadata <- bundle$metadata
  counts <- bundle$counts
  comparison <- comparisons[[dataset]]
  unit_col <- sample_column(metadata)

  macrophage_cells <- which(metadata$compartment == "Macrophage")
  hsc_cells <- which(metadata$compartment == "HSC")
  macrophage_metadata <- metadata[macrophage_cells, , drop = FALSE]
  hsc_metadata <- metadata[hsc_cells, , drop = FALSE]
  macrophage_counts <- counts[, macrophage_cells, drop = FALSE]
  hsc_counts <- counts[, hsc_cells, drop = FALSE]

  macrophage_pb <- aggregate_counts(macrophage_counts, macrophage_metadata, unit_col)
  hsc_pb <- aggregate_counts(hsc_counts, hsc_metadata, unit_col)
  macrophage_sample_metadata <- macrophage_metadata[match(colnames(macrophage_pb), macrophage_metadata[[unit_col]]), , drop = FALSE]
  hsc_sample_metadata <- hsc_metadata[match(colnames(hsc_pb), hsc_metadata[[unit_col]]), , drop = FALSE]

  macrophage_de <- run_edger(macrophage_pb, macrophage_sample_metadata, comparison$column, comparison$case, comparison$control)
  hsc_de <- run_edger(hsc_pb, hsc_sample_metadata, comparison$column, comparison$case, comparison$control)
  write.table(macrophage_de, file.path(results_dir, paste0("025_", dataset, "_macrophage_patient_pseudobulk_DE.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(hsc_de, file.path(results_dir, paste0("026_", dataset, "_HSC_patient_pseudobulk_DE.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

  expressed_macrophage <- expressed_genes(macrophage_counts)
  expressed_hsc <- expressed_genes(hsc_counts)
  potential_ligands <- intersect(
    unique(lr_network$from[lr_network$to %in% expressed_hsc]),
    expressed_macrophage
  )
  background <- intersect(expressed_hsc, colnames(ligand_target_matrix))
  geneset <- hsc_de$gene[hsc_de$logFC > 0.25 & hsc_de$PValue < 0.10]
  geneset <- intersect(geneset, background)
  if (length(geneset) < 20L) {
    geneset <- hsc_de$gene[hsc_de$logFC > 0][order(hsc_de$PValue[hsc_de$logFC > 0])]
    geneset <- head(intersect(geneset, background), 150L)
  }
  if (length(geneset) < 10L || length(potential_ligands) < 5L) {
    stop(dataset, ": insufficient genes or ligands for NicheNet")
  }

  activities <- predict_ligand_activities(
    geneset = geneset,
    background_expressed_genes = background,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  )
  activities <- activities[order(activities$aupr_corrected, decreasing = TRUE), , drop = FALSE]
  activities$macrophage_logFC <- macrophage_de$logFC[match(activities$test_ligand, macrophage_de$gene)]
  activities$macrophage_FDR <- macrophage_de$FDR[match(activities$test_ligand, macrophage_de$gene)]
  write.table(activities, file.path(results_dir, paste0("027_", dataset, "_NicheNet_ligand_activities.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

  top_ligands <- head(activities$test_ligand, 30L)
  links <- do.call(
    rbind,
    lapply(
      top_ligands,
      get_weighted_ligand_target_links,
      geneset = geneset,
      ligand_target_matrix = ligand_target_matrix,
      n = 250
    )
  )
  hsc_targets <- c("YAP1", "NUAK2", "COL1A1", "COL1A2", "COL3A1", "LOX", "ACTA2", "TAGLN", "TGFB1", "PDGFRA", "PDGFRB")
  links$HSC_activation_target <- links$target %in% hsc_targets
  write.table(links, file.path(results_dir, paste0("028_", dataset, "_NicheNet_ligand_target_links.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

  plot_data <- head(activities, 20L)
  plot_data$test_ligand <- factor(plot_data$test_ligand, levels = rev(plot_data$test_ligand))
  plot <- ggplot(plot_data, aes(x = aupr_corrected, y = test_ligand, fill = macrophage_logFC)) +
    geom_col() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey70") +
    labs(title = paste(dataset, "macrophage ligands explaining HSC response"), x = "NicheNet corrected AUPR", y = NULL, fill = "Macrophage\nlogFC") +
    theme_bw(base_size = 11)
  ggsave(file.path(figures_dir, paste0("009_", dataset, "_NicheNet_ligand_activity.png")), plot, width = 8, height = 7, dpi = 240)
}
