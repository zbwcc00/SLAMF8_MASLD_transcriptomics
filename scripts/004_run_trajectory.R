#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(monocle3)
  library(slingshot)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
datasets <- if (length(args)) args else c("GSE344087", "GSE298719", "GSE212837")
seed <- 20260920L
project_root <- normalizePath(Sys.getenv("SLAMF8_MASLD_ROOT", "."), winslash = "/", mustWork = FALSE)
processed_dir <- file.path(project_root, "data_processed")
results_dir <- file.path(project_root, "results")
figures_dir <- file.path(project_root, "figures")

write_tsv_gz <- function(data, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection))
  write.table(data, connection, sep = "\t", row.names = FALSE, quote = FALSE)
}

read_bundle <- function(dataset) {
  prefix <- file.path(processed_dir, paste0(dataset, "_trajectory"))
  counts <- as(readMM(gzfile(paste0(prefix, "_counts.mtx.gz"))), "dgCMatrix")
  features <- read.delim(gzfile(paste0(prefix, "_features.tsv.gz")), check.names = FALSE)
  metadata <- read.delim(gzfile(paste0(prefix, "_cell_metadata.tsv.gz")), check.names = FALSE)
  rownames(counts) <- make.unique(as.character(features$var_name))
  colnames(counts) <- metadata$barcode
  rownames(metadata) <- metadata$barcode
  list(counts = counts, features = features, metadata = metadata, prefix = prefix)
}

choose_label <- function(labels, pattern, fallback = NULL) {
  hits <- labels[grepl(pattern, labels, ignore.case = TRUE)]
  if (length(hits)) return(hits[[1L]])
  fallback
}

choose_start_label <- function(labels) {
  exact_hits <- labels[grepl("FCN1|VCAN", labels, ignore.case = TRUE)]
  if (length(exact_hits)) return(exact_hits[[1L]])
  choose_label(labels, "monocyte", fallback = labels[[1L]])
}

root_principal_node <- function(cds, root_cells) {
  closest <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest <- as.matrix(closest[colnames(cds), , drop = FALSE])
  root_positions <- match(root_cells, rownames(closest), nomatch = 0L)
  root_positions <- root_positions[root_positions > 0L]
  if (!length(root_positions)) stop("No root cells mapped to principal graph")
  vertex <- names(which.max(table(closest[root_positions, 1L])))
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.integer(vertex)]
}

run_monocle <- function(dataset, bundle) {
  counts <- bundle$counts
  metadata <- bundle$metadata
  gene_metadata <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
  keep_genes <- Matrix::rowSums(counts > 0) >= 10L
  cds <- new_cell_data_set(
    counts[keep_genes, , drop = FALSE],
    cell_metadata = metadata,
    gene_metadata = gene_metadata[keep_genes, , drop = FALSE]
  )
  set.seed(seed)
  cds <- preprocess_cds(cds, num_dim = 50, method = "PCA")
  if (length(unique(colData(cds)$gsm)) > 1L) {
    cds <- align_cds(cds, alignment_group = "gsm")
  }
  cds <- reduce_dimension(cds, reduction_method = "UMAP", umap.fast_sgd = FALSE, cores = 1)
  cds <- cluster_cells(cds, reduction_method = "UMAP", resolution = 1e-3)
  cds <- learn_graph(cds, use_partition = TRUE, close_loop = FALSE)
  labels <- unique(as.character(colData(cds)$subtype))
  start_label <- choose_start_label(labels)
  root_cells <- rownames(colData(cds))[as.character(colData(cds)$subtype) == start_label]
  root_node <- root_principal_node(cds, root_cells)
  cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes = root_node)
  pseudotime <- monocle3::pseudotime(cds)
  output <- data.frame(
    barcode = names(pseudotime),
    monocle3_pseudotime = as.numeric(pseudotime),
    stringsAsFactors = FALSE
  )
  write_tsv_gz(output, file.path(results_dir, paste0("007_", dataset, "_monocle3_pseudotime.tsv.gz")))
  saveRDS(cds, file.path(processed_dir, paste0(dataset, "_monocle3_cds.rds")))

  subtype_plot <- plot_cells(
    cds,
    color_cells_by = "subtype",
    label_groups_by_cluster = FALSE,
    label_branch_points = FALSE,
    label_leaves = FALSE,
    cell_size = 0.35
  ) + ggtitle(paste(dataset, "Monocle3 trajectory by subtype"))
  ggsave(file.path(figures_dir, paste0("004A_", dataset, "_monocle3_subtypes.png")), subtype_plot, width = 9, height = 7, dpi = 240)
  time_plot <- plot_cells(
    cds,
    color_cells_by = "pseudotime",
    label_groups_by_cluster = FALSE,
    label_branch_points = TRUE,
    label_leaves = TRUE,
    cell_size = 0.35
  ) + ggtitle(paste(dataset, "Monocle3 pseudotime"))
  ggsave(file.path(figures_dir, paste0("004B_", dataset, "_monocle3_pseudotime.png")), time_plot, width = 9, height = 7, dpi = 240)
  list(cds = cds, output = output, start_label = start_label)
}

run_slingshot <- function(dataset, bundle) {
  metadata <- bundle$metadata
  pca_path <- paste0(bundle$prefix, "_pca_harmony.tsv.gz")
  if (!file.exists(pca_path)) pca_path <- paste0(bundle$prefix, "_pca.tsv.gz")
  embedding <- read.delim(gzfile(pca_path), check.names = FALSE)
  rownames(embedding) <- embedding$barcode
  embedding$barcode <- NULL
  embedding <- as.matrix(embedding[metadata$barcode, , drop = FALSE])
  labels <- unique(as.character(metadata$subtype))
  start_label <- choose_start_label(labels)
  end_label <- choose_label(labels, "TREM2|OLR1|SAMac", fallback = NULL)
  sce <- SingleCellExperiment(
    assays = list(counts = bundle$counts),
    colData = metadata
  )
  reducedDim(sce, "PCA") <- embedding
  set.seed(seed)
  if (is.null(end_label) || identical(start_label, end_label)) {
    sce <- slingshot(sce, clusterLabels = "subtype", reducedDim = "PCA", start.clus = start_label)
  } else {
    sce <- slingshot(
      sce,
      clusterLabels = "subtype",
      reducedDim = "PCA",
      start.clus = start_label,
      end.clus = end_label
    )
  }
  times <- slingPseudotime(sce, na = FALSE)
  weights <- slingCurveWeights(sce)
  output <- data.frame(barcode = colnames(sce), times, check.names = FALSE)
  write_tsv_gz(output, file.path(results_dir, paste0("008_", dataset, "_slingshot_pseudotime.tsv.gz")))
  saveRDS(sce, file.path(processed_dir, paste0(dataset, "_slingshot_sce.rds")))

  png(file.path(figures_dir, paste0("005_", dataset, "_slingshot_curves.png")), width = 2200, height = 1800, res = 240)
  palette <- grDevices::hcl.colors(length(labels), "Dynamic")
  colors <- palette[match(as.character(metadata$subtype), labels)]
  plot(embedding[, 1L:2L], col = colors, pch = 16, cex = 0.35, xlab = "PC1", ylab = "PC2", main = paste(dataset, "Slingshot lineages"))
  lines(SlingshotDataSet(sce), lwd = 3, col = "black")
  legend("topright", legend = labels, col = palette, pch = 16, cex = 0.7, bty = "n")
  dev.off()
  list(sce = sce, times = times, weights = weights, start_label = start_label, end_label = end_label)
}

run_tradeseq <- function(dataset, bundle, slingshot_result) {
  if (!requireNamespace("tradeSeq", quietly = TRUE)) {
    warning("tradeSeq unavailable; skipping dynamic-gene model for ", dataset)
    return(NULL)
  }
  counts <- bundle$counts
  times <- slingshot_result$times
  weights <- slingshot_result$weights
  usable <- which(rowSums(weights, na.rm = TRUE) > 0 & rowSums(!is.na(times)) > 0)
  if (length(usable) > 6000L) {
    set.seed(seed)
    groups <- split(usable, bundle$metadata$subtype[usable])
    allocations <- pmax(1L, round(vapply(groups, length, integer(1L)) / length(usable) * 6000L))
    selected_cells <- unlist(Map(function(index, number) sample(index, min(length(index), number)), groups, allocations), use.names = FALSE)
  } else {
    selected_cells <- usable
  }
  means <- Matrix::rowMeans(counts[, selected_cells, drop = FALSE])
  second <- Matrix::rowMeans(counts[, selected_cells, drop = FALSE]^2)
  dispersion <- pmax(second - means^2, 0) / (means + 0.01)
  candidates <- c("SLAMF8", "FCN1", "VCAN", "TREM2", "CD9", "GPNMB", "SPP1", "OLR1", "IL1B", "LGALS3", "ACSL4", "HMOX1", "FTH1", "FTL", "PTGS2", "CDKN1A")
  variable_genes <- names(sort(dispersion, decreasing = TRUE))
  variable_genes <- variable_genes[make.names(variable_genes) == variable_genes]
  selected_genes <- unique(c(variable_genes[seq_len(min(2000L, length(variable_genes)))], intersect(candidates, rownames(counts))))
  set.seed(seed)
  model <- tradeSeq::fitGAM(
    counts = counts[selected_genes, selected_cells, drop = FALSE],
    pseudotime = times[selected_cells, , drop = FALSE],
    cellWeights = weights[selected_cells, , drop = FALSE],
    nknots = 6,
    verbose = TRUE,
    parallel = FALSE
  )
  association <- tradeSeq::associationTest(model)
  association$gene <- rownames(association)
  association$fdr <- p.adjust(association$pvalue, method = "BH")
  start_end <- tradeSeq::startVsEndTest(model)
  start_end$gene <- rownames(start_end)
  start_end$fdr <- p.adjust(start_end$pvalue, method = "BH")
  write_tsv_gz(association, file.path(results_dir, paste0("009_", dataset, "_tradeSeq_association.tsv.gz")))
  write_tsv_gz(start_end, file.path(results_dir, paste0("010_", dataset, "_tradeSeq_start_vs_end.tsv.gz")))
  saveRDS(model, file.path(processed_dir, paste0(dataset, "_tradeSeq_fitGAM.rds")))
  invisible(model)
}

for (dataset in datasets) {
  message("Trajectory analysis: ", dataset)
  bundle <- read_bundle(dataset)
  monocle_result <- run_monocle(dataset, bundle)
  slingshot_result <- run_slingshot(dataset, bundle)
  run_tradeseq(dataset, bundle, slingshot_result)
  summary <- data.frame(
    dataset = dataset,
    n_cells = ncol(bundle$counts),
    n_genes = nrow(bundle$counts),
    monocle_start = monocle_result$start_label,
    slingshot_start = slingshot_result$start_label,
    slingshot_end = ifelse(is.null(slingshot_result$end_label), NA, slingshot_result$end_label),
    stringsAsFactors = FALSE
  )
  write.table(summary, file.path(results_dir, paste0("011_", dataset, "_trajectory_summary.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
}
