#!/usr/bin/env Rscript
# ============================================================
# Gene Expression Metrics & Dot Plot Visualization
#
# Description:
#   Calculates and visualizes gene expression metrics from a Seurat RDS dataset.
#   - Computes Average Expression (LogNormalized) and Detection Rate (Percent Expressed).
#   - Groups cells by specified metadata (e.g., Cell Type, or Cell Type x Stage).
#   - Handles missing 'counts' slots by using 'data' slots as proxies.
#   - Generates a Dot Plot where:
#       Size  = Detection Rate (0-1)
#       Color = Mean Expression (Plasma scale)
#
# Output:
#   - PDF file containing the dot plot for specified genes (AGI).
# ============================================================

# Load necessary libraries
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
# B) Function: Plot and Save PDF from AGI Vector
#    * Color bar: plasma (Low=Purple -> High=Yellow)
# ==============================================================================
plot_and_save_from_AGI <- function(agi_list, gm, outdir = "gene_metrics_out", file_stub = "AGI_batch") {
  
  stopifnot(is.list(gm), !is.null(gm$avg_mat), !is.null(gm$pct_mat), !is.null(gm$groups))
  
  genes_all <- rownames(gm$avg_mat)
  agi_in <- intersect(agi_list, genes_all)
  agi_miss <- setdiff(agi_list, agi_in)
  
  if (length(agi_miss) > 0) {
    warning("Genes not found: ", paste(agi_miss, collapse = ", "))
  }
  if (length(agi_in) == 0) stop("No valid genes available for plotting.")
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # Reshape to long format
  avg_df <- as.data.frame(gm$avg_mat[agi_in, , drop = FALSE])
  avg_df$AGI <- rownames(avg_df)
  avg_long <- melt(avg_df, id.vars = "AGI", variable.name = "Group", value.name = "MeanExpr")
  
  pct_df <- as.data.frame(gm$pct_mat[agi_in, , drop = FALSE])
  pct_df$AGI <- rownames(pct_df)
  pct_long <- melt(pct_df, id.vars = "AGI", variable.name = "Group", value.name = "Pct")
  
  df <- merge(avg_long, pct_long, by = c("AGI", "Group"))
  
  # Fix factor levels
  df$Group <- factor(df$Group, levels = gm$groups)
  df$AGI   <- factor(df$AGI, levels = rev(agi_in))  # Order as specified (reversed for plot)
  
  # Generate Dot Plot (Color = Mean Expression, Size = Detection Rate)
  p <- ggplot(df, aes(x = Group, y = AGI)) +
    geom_point(aes(size = Pct, color = MeanExpr)) +
    scale_size_continuous(name = "Detection rate", range = c(0, 7), limits = c(0, 1)) +
    scale_color_viridis_c(name = "Mean expr", option = "plasma") +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = sprintf("Gene metrics (%s)", gm$grouping_desc),
      x = "Group", y = "Gene ID",
      caption = "Color: Avg log-normalized expression (plasma) / Size: Detection rate (counts > 0)"
    )
  
  # Calculate PDF dimensions automatically
  W_PER_GROUP <- 0.35
  H_PER_GENE  <- 0.45
  MIN_W <- 6.5
  MIN_H <- 5.0
  pdf_w <- max(MIN_W, length(gm$groups) * W_PER_GROUP * 1.6)
  pdf_h <- max(MIN_H, length(agi_in)    * H_PER_GENE  * 1.6)
  
  dot_pdf <- file.path(outdir, paste0(file_stub, "_dot.pdf"))
  ggsave(dot_pdf, p, width = pdf_w, height = pdf_h, units = "in", limitsize = FALSE)
  
  message("Saved PDF: ", normalizePath(dot_pdf))
  
  invisible(list(
    pdf = dot_pdf,
    data_long = df,
    genes_drawn = agi_in,
    genes_missing = agi_miss
  ))
}

# ==============================================================================
# Main Execution Example
# ==============================================================================

# 1) Load Seurat Object
x <- readRDS(RDS_PATH)

# 2) Calculate metrics (Example: Aggregate by celltype only)
gm <- build_gene_metrics(x, celltype_col = CELLTYPE_COL, stage_col = STAGE_COL)

# 3) Set output directory
outdir <- dirname(RDS_PATH)
outdir <- file.path(outdir, "gene_metrics_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 4) Plot and save using the provided AGI vector
agi_list <- c("AT2G36270", "AT3G56850", "AT3G44460", "AT2G41070", "AT2G35530", "AT1G32150")
res <- plot_and_save_from_AGI(agi_list, gm, outdir = outdir, file_stub = "AGI_batch")

message("Output directory: ", normalizePath(outdir))