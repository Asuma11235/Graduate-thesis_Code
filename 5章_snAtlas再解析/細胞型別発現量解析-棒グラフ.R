#!/usr/bin/env Rscript

# ============================================================
# Absolute Gene Expression Bar Plotter
#
# Description:
#   Calculates and visualizes absolute gene expression levels from a Seurat RDS dataset.
#   - Computes Average Expression (LogNormalized) per group.
#   - Generates a Bar Plot comparing mean expression across groups for specified genes.
#   - Supports optional log10 transformation for the Y-axis.
#
# Output:
#   - PDF file containing the bar plot.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(reshape2)
})

# ==============================================================================
# Configuration
# ==============================================================================
# Path to the input Seurat object (.RDS file)
RDS_PATH     <- "path/to/your/input_data.rds"

# Metadata column names
CELLTYPE_COL <- "integrated_annotation" 
STAGE_COL    <- NULL  # Set to NULL if stage information is not used

# ==============================================================================
# A) Function: Construct Gene Metrics
# ==============================================================================
build_gene_metrics <- function(x, celltype_col, stage_col = NULL) {
  stopifnot(inherits(x, "Seurat"))
  
  # Ensure RNA assay is active and normalized data exists
  DefaultAssay(x) <- "RNA"
  if (!"data" %in% slotNames(x[["RNA"]])) {
    message("No 'data' slot found. NormalizeData(LogNormalize, scale.factor=1e4).")
    x <- NormalizeData(x, normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = FALSE)
  }
  
  # Set Idents and define grouping column
  stopifnot(celltype_col %in% colnames(x@meta.data))
  Idents(x) <- factor(x@meta.data[[celltype_col]])
  
  if (!is.null(stage_col)) {
    stopifnot(stage_col %in% colnames(x@meta.data))
    x$.__group__ <- interaction(x@meta.data[[stage_col]],
                                as.character(Idents(x)),
                                drop = TRUE, sep = "|")
    grouping_desc <- paste0(stage_col, " x ", celltype_col)
  } else {
    x$.__group__ <- factor(as.character(Idents(x)))
    grouping_desc <- celltype_col
  }
  groups <- levels(x$.__group__)
  
  # Calculate Average Expression (slot = 'data')
  avg_obj <- AverageExpression(
    x, group.by = ".__group__", assays = DefaultAssay(x),
    slot = "data", verbose = FALSE
  )
  
  # Handle list output from AverageExpression
  if (is.list(avg_obj)) {
    if (DefaultAssay(x) %in% names(avg_obj)) {
      avg_mat <- avg_obj[[DefaultAssay(x)]]
    } else {
      avg_mat <- avg_obj[[1]]
    }
  } else {
    avg_mat <- avg_obj
  }
  
  # Convert to matrix
  if (inherits(avg_mat, "data.frame") || inherits(avg_mat, "tbl_df")) avg_mat <- as.matrix(avg_mat)
  if (inherits(avg_mat, "Matrix") || inherits(avg_mat, "dgCMatrix")) avg_mat <- as.matrix(avg_mat)
  if (!is.matrix(avg_mat)) stop("Cannot convert AverageExpression output to matrix: ", paste(class(avg_mat), collapse=", "))
  
  if (is.null(rownames(avg_mat))) rownames(avg_mat) <- rownames(x)
  if (is.null(colnames(avg_mat))) colnames(avg_mat) <- groups
  
  # Calculate Detection Rate (counts > 0). Use data > 0 if counts slot is empty.
  cnt <- GetAssayData(x, assay = DefaultAssay(x), slot = "counts")
  if (Matrix::nnzero(cnt) == 0) {
    message("counts slot is empty. Using data > 0 as detection proxy.")
    dat <- GetAssayData(x, assay = DefaultAssay(x), slot = "data")
    pct_mat <- sapply(groups, function(g) {
      cells_g <- colnames(x)[x$.__group__ == g]
      Matrix::rowMeans(dat[, cells_g, drop = FALSE] > 0)
    })
  } else {
    pct_mat <- sapply(groups, function(g) {
      cells_g <- colnames(x)[x$.__group__ == g]
      Matrix::rowMeans(cnt[, cells_g, drop = FALSE] > 0)
    })
  }
  
  # Format detection rate matrix
  if (is.null(dim(pct_mat))) {
    pct_mat <- matrix(pct_mat, ncol = 1, dimnames = list(rownames(x), groups))
  } else {
    rownames(pct_mat) <- rownames(x)
    colnames(pct_mat) <- groups
  }
  
  list(
    x = x,
    grouping_desc = grouping_desc,
    groups = groups,
    avg_mat = avg_mat,   # Row = Gene, Col = Group
    pct_mat = pct_mat    # Row = Gene, Col = Group (0-1)
  )
}

# ==============================================================================
# B) Function: Plot and Save Absolute Expression Bar Chart
#    - X-axis: Group
#    - Y-axis: Mean Expression (from 'data' slot)
#    - Fill: Gene (AGI)
# ==============================================================================
plot_bar_abs_expression <- function(agi_list,
                                    gm,
                                    outdir      = "gene_metrics_out",
                                    file_stub   = "AGI_bar_abs",
                                    y_transform = c("none", "log10"),
                                    pseudo_count = 0) {
  
  stopifnot(is.list(gm), !is.null(gm$avg_mat), !is.null(gm$groups))
  y_transform <- match.arg(y_transform)
  
  genes_all <- rownames(gm$avg_mat)
  agi_in    <- intersect(agi_list, genes_all)
  agi_miss  <- setdiff(agi_list, genes_all)
  
  if (length(agi_miss) > 0)
    warning("Genes not found: ", paste(agi_miss, collapse = ", "))
  if (length(agi_in) == 0)
    stop("No valid genes available for plotting.")
  
  # Reshape to long format
  avg_df <- as.data.frame(gm$avg_mat[agi_in, , drop = FALSE])
  avg_df$AGI <- rownames(avg_df)
  df <- reshape2::melt(avg_df, id.vars = "AGI",
                       variable.name = "Group", value.name = "MeanExpr")
  
  # Fix factor levels
  df$Group <- factor(df$Group, levels = gm$groups)
  df$AGI   <- factor(df$AGI, levels = agi_in)
  
  # Create Base Plot
  p <- ggplot(df, aes(x = Group, y = MeanExpr, fill = AGI)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.minor = element_blank()) +
    labs(title = sprintf("Absolute expression (%s)", gm$grouping_desc),
         x = "Group", y = "Mean expression (log-normalized)")
  
  # Optional: Log Transformation
  if (y_transform == "log10") {
    p <- p +
      aes(y = MeanExpr + pseudo_count) +
      scale_y_continuous(trans = "log10") +
      labs(y = sprintf("Mean expression (log10, +%g)", pseudo_count))
  }
  
  # Save PDF
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  W_PER_GROUP <- 0.45
  MIN_W <- 6.5
  pdf_w <- max(MIN_W, length(gm$groups) * W_PER_GROUP * 1.6)
  pdf_h <- 5.0
  
  out_pdf <- file.path(outdir, paste0(file_stub, ".pdf"))
  ggsave(out_pdf, p, width = pdf_w, height = pdf_h, units = "in", limitsize = FALSE)
  
  message("Saved PDF: ", normalizePath(out_pdf))
  
  invisible(list(pdf = out_pdf, data = df,
                 genes_drawn = agi_in, genes_missing = agi_miss))
}

# ==============================================================================
# Main Execution Example
# ==============================================================================

# 1) Load Seurat Object
x <- readRDS(RDS_PATH)

# 2) Calculate gene metrics
gm <- build_gene_metrics(x, celltype_col = CELLTYPE_COL, stage_col = STAGE_COL)

# 3) Set output directory
outdir <- file.path(dirname(RDS_PATH), "gene_metrics_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 4) Define target genes
# Example: GLK1-PIF4-BIL1 pathway
agi_list <- c("AT2G20570", "AT2G43010", "AT1G75080")

# 5) Generate Absolute Expression Bar Plot (Standard Scale)
res_abs <- plot_bar_abs_expression(agi_list, gm,
                                   outdir = outdir,
                                   file_stub = "AGI_bar_abs",
                                   y_transform = "none")

message("Output directory: ", normalizePath(outdir))