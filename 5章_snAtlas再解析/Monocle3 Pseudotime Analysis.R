#!/usr/bin/env Rscript

# ============================================================
# Monocle3 Pseudotime Analysis Pipeline (Seurat v5 Compatible)
#
# Description:
#   Performs trajectory inference and pseudotime analysis using Monocle3.
#   - Compatible with Seurat v5 object structure.
#   - Restores raw "pseudo-counts" from log-normalized data (inverse log transformation).
#   - Converts Seurat object to Monocle3 CellDataSet (CDS).
#   - Learns trajectory graph and allows interactive root selection.
#   - Visualizes pseudotime trajectories and gene expression trends.
#
# Input:
#   - Seurat RDS file (Log-normalized).
#
# Output:
#   - UMAP plots (Annotation, Pseudotime, Gene Expression).
#   - Gene expression trend plots along pseudotime.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(monocle3)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
})

# ==============================================================================
# Configuration
# ==============================================================================
# Input/Output
SEURAT_RDS    <- "path/to/your/scRNA_data.rds"
OUT_DIR       <- dirname(SEURAT_RDS)

# Target Gene for Trend Analysis
TARGET_GENE   <- "AT3G23430"

# Metadata Column for Cell Types
CELLTYPE_COL  <- "final_annotation"

# ==============================================================================
# 1. Load Data
# ==============================================================================
if (!file.exists(SEURAT_RDS)) stop("Seurat file not found: ", SEURAT_RDS)

seurat_obj <- readRDS(SEURAT_RDS)
seurat_obj <- UpdateSeuratObject(seurat_obj)

# ==============================================================================
# 2. Restore Pseudo-Counts (Crucial for Monocle3)
# ==============================================================================
message("Checking data layers and restoring pseudo-counts...")

DefaultAssay(seurat_obj) <- "RNA"

# Retrieve normalized data (Compatible with Seurat v5/v4)
norm_data <- tryCatch({
  if (is(seurat_obj[["RNA"]], "Assay5")) {
    LayerData(seurat_obj, assay = "RNA", layer = "data")
  } else {
    GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  }
}, error = function(e) {
  stop("Failed to retrieve 'data' layer/slot from RNA assay.")
})

# Inverse Log-Normalization: exp(x) - 1
# This restores approximate raw counts required by Monocle3
message("Restoring pseudo-counts (exp(x) - 1)...")
pseudo_counts <- exp(norm_data) - 1
pseudo_counts <- round(pseudo_counts) # Round to nearest integer

# Ensure sparse matrix format
if (!inherits(pseudo_counts, "dgCMatrix")) {
  pseudo_counts <- as(pseudo_counts, "dgCMatrix")
}

# Write back to 'counts' layer (Seurat v5 compatible)
if (is(seurat_obj[["RNA"]], "Assay5")) {
  LayerData(seurat_obj, assay = "RNA", layer = "counts") <- pseudo_counts
} else {
  seurat_obj <- SetAssayData(seurat_obj, assay = "RNA", slot = "counts", new.data = pseudo_counts)
}

message("Pseudo-count restoration complete.")
message("Max count value: ", max(pseudo_counts), " (Should be >> 10)")

# ==============================================================================
# 3. Convert to Monocle3 CDS
# ==============================================================================
message("Converting Seurat object to Monocle3 CDS...")
cds <- as.cell_data_set(seurat_obj)

# Recalculate size factors using restored counts
cds <- estimate_size_factors(cds)

# ==============================================================================
# 4. Trajectory Inference
# ==============================================================================
# Assign gene short names for plotting
rowData(cds)$gene_short_name <- rownames(cds)

# Cluster cells and learn graph
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

# ==============================================================================
# 5. Interactive Root Selection
# ==============================================================================
message("Plotting UMAP for root selection...")

p_root <- plot_cells(
  cds,
  color_cells_by = CELLTYPE_COL,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 3
) + ggtitle("Select Root Nodes")

print(p_root)

message("\n>>> ACTION REQUIRED: Select root nodes in the pop-up window and click 'Done'. <<<\n")
cds <- order_cells(cds)

# ==============================================================================
# 6. Visualization & Export
# ==============================================================================
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- 6.1 UMAP Colored by Annotation ---
if (!CELLTYPE_COL %in% colnames(colData(cds))) {
  warning("Column '", CELLTYPE_COL, "' not found in colData. Skipping annotation plot.")
} else {
  p_anno <- plot_cells(
    cds,
    color_cells_by = CELLTYPE_COL,
    label_cell_groups = TRUE,
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    group_label_size = 4,
    cell_size = 1
  ) +
    ggtitle("UMAP: Cell Type Annotation") +
    theme_minimal()
  
  print(p_anno)
  
  out_anno <- file.path(OUT_DIR, "umap_final_annotation_labeled.png")
  ggsave(out_anno, p_anno, width = 10, height = 8)
  message("Saved: ", normalizePath(out_anno))
}

# --- 6.2 UMAP Colored by Pseudotime ---
p_pt <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 3,
  cell_size = 1
) +
  scale_color_viridis_c(option = "plasma", na.value = "grey80", name = "Pseudotime") +
  ggtitle("UMAP: Pseudotime Trajectory") +
  theme_minimal()

print(p_pt)

out_pt <- file.path(OUT_DIR, "umap_pseudotime_corrected.png")
ggsave(out_pt, p_pt, width = 10, height = 8)
message("Saved: ", normalizePath(out_pt))

# --- 6.3 Gene Expression Trend along Pseudotime ---
if (TARGET_GENE %in% rownames(cds)) {
  
  # Plot Trend (Scatter + Spline)
  p_trend <- plot_genes_in_pseudotime(
    cds[TARGET_GENE, ],
    color_cells_by = "pseudotime",
    min_expr = 0,
    cell_size = 1.0,
    trend_formula = "~ splines::ns(pseudotime, df=3)"
  ) +
    scale_color_viridis_c(option = "plasma", name = "Pseudotime") +
    ggtitle(paste0("Expression Trend: ", TARGET_GENE)) +
    theme_minimal()
  
  print(p_trend)
  
  out_trend <- file.path(OUT_DIR, "gene_pseudotime_corrected.png")
  ggsave(out_trend, p_trend, width = 10, height = 6)
  message("Saved: ", normalizePath(out_trend))
  
  # Plot Expression on UMAP
  p_umap_gene <- plot_cells(
    cds,
    genes = TARGET_GENE,
    show_trajectory_graph = FALSE,
    cell_size = 1
  ) +
    scale_color_viridis_c(option = "plasma", name = "LogExpr") +
    ggtitle(paste0("Expression on UMAP: ", TARGET_GENE)) +
    theme_minimal()
  
  print(p_umap_gene)
  
  out_umap_gene <- file.path(OUT_DIR, "gene_umap_corrected.png")
  ggsave(out_umap_gene, p_umap_gene, width = 10, height = 8)
  message("Saved: ", normalizePath(out_umap_gene))
  
} else {
  warning("Target gene '", TARGET_GENE, "' not found in the dataset.")
}

message("Pipeline execution complete.")