#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(monocle3)
  library(slingshot)
  library(igraph)
})

project_root <- normalizePath(Sys.getenv("SLAMF8_MASLD_ROOT", "."), winslash = "/", mustWork = FALSE)
processed_dir <- file.path(project_root, "data_processed")
results_dir <- file.path(project_root, "results")
figures_dir <- file.path(project_root, "figures")
datasets <- c("GSE344087", "GSE298719", "GSE212837")

broad_color <- function(label) {
  vapply(label, function(item) {
    value <- tolower(as.character(item))
    if (grepl("contaminant|low_quality|ambiguous", value)) return("#BDBDBD")
    if (grepl("monocyte|fcn1|fcgr3a", value)) return("#E69F00")
    if (grepl("kupffer|marco", value)) return("#0072B2")
    if (grepl("trem2|spp1|gpnmb|macrophage|samac", value)) return("#D55E00")
    if (grepl("dendritic|cdc|pdc|lamp3", value)) return("#009E73")
    if (grepl("cycling", value)) return("#CC79A7")
    "#666666"
  }, character(1))
}

for (dataset in datasets) {
  cds <- readRDS(file.path(processed_dir, paste0(dataset, "_monocle3_cds.rds")))
  umap <- reducedDims(cds)$UMAP
  graph <- principal_graph(cds)[["UMAP"]]
  nodes <- principal_graph_aux(cds)[["UMAP"]]$dp_mst
  edge_table <- as.data.frame(as_edgelist(graph, names = TRUE), stringsAsFactors = FALSE)
  colnames(edge_table) <- c("from", "to")
  edge_table$x1 <- nodes[1, match(edge_table$from, colnames(nodes))]
  edge_table$y1 <- nodes[2, match(edge_table$from, colnames(nodes))]
  edge_table$x2 <- nodes[1, match(edge_table$to, colnames(nodes))]
  edge_table$y2 <- nodes[2, match(edge_table$to, colnames(nodes))]
  write.table(edge_table[, c("from", "to", "x1", "y1", "x2", "y2")],
              file.path(results_dir, paste0("037_", dataset, "_monocle3_graph_edges.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  sce <- readRDS(file.path(processed_dir, paste0(dataset, "_slingshot_sce.rds")))
  curves <- slingCurves(SlingshotDataSet(sce))
  curve_rows <- list()
  for (index in seq_along(curves)) {
    curve <- curves[[index]]
    coords <- as.data.frame(curve$s[, 1:2, drop = FALSE])
    colnames(coords) <- c("x", "y")
    coords$lineage <- paste0("Lineage", index)
    curve_rows[[index]] <- coords
  }
  curves_table <- do.call(rbind, curve_rows)
  write.table(curves_table, file.path(results_dir, paste0("038_", dataset, "_slingshot_curves.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

png(file.path(figures_dir, "Supplementary_Figure3_Slingshot_lineages.png"), width = 3600, height = 1250, res = 300)
par(mfrow = c(1, 3), mar = c(4.0, 4.0, 3.0, 1.0), oma = c(0, 0, 2, 0))
for (dataset in datasets) {
  sce <- readRDS(file.path(processed_dir, paste0(dataset, "_slingshot_sce.rds")))
  embedding <- reducedDim(sce, "PCA")
  labels <- as.character(colData(sce)$subtype)
  curves <- slingCurves(SlingshotDataSet(sce))
  plot(embedding[, 1:2], col = broad_color(labels), pch = 16, cex = 0.28,
       xlab = "PCA 1", ylab = "PCA 2", main = dataset, axes = TRUE)
  for (curve in curves) lines(curve$s[, 1:2], col = "#17202A", lwd = 1.8)
  legend("topright", legend = c("Monocyte", "Kupffer", "Macrophage/SAMac", "Dendritic", "Contaminant"),
         col = c("#E69F00", "#0072B2", "#D55E00", "#009E73", "#BDBDBD"), pch = 16, cex = 0.6,
         bty = "n", title = paste0(length(curves), " lineages"))
}
mtext("Slingshot lineage curves in PCA space", outer = TRUE, cex = 1.35, font = 2)
dev.off()

pdf(file.path(figures_dir, "Supplementary_Figure3_Slingshot_lineages.pdf"), width = 12, height = 4.2, useDingbats = FALSE)
par(mfrow = c(1, 3), mar = c(4.0, 4.0, 3.0, 1.0), oma = c(0, 0, 2, 0))
for (dataset in datasets) {
  sce <- readRDS(file.path(processed_dir, paste0(dataset, "_slingshot_sce.rds")))
  embedding <- reducedDim(sce, "PCA")
  labels <- as.character(colData(sce)$subtype)
  curves <- slingCurves(SlingshotDataSet(sce))
  plot(embedding[, 1:2], col = broad_color(labels), pch = 16, cex = 0.28,
       xlab = "PCA 1", ylab = "PCA 2", main = dataset, axes = TRUE)
  for (curve in curves) lines(curve$s[, 1:2], col = "#17202A", lwd = 1.8)
  legend("topright", legend = c("Monocyte", "Kupffer", "Macrophage/SAMac", "Dendritic", "Contaminant"),
         col = c("#E69F00", "#0072B2", "#D55E00", "#009E73", "#BDBDBD"), pch = 16, cex = 0.6,
         bty = "n", title = paste0(length(curves), " lineages"))
}
mtext("Slingshot lineage curves in PCA space", outer = TRUE, cex = 1.35, font = 2)
dev.off()
