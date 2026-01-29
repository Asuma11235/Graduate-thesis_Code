#!/usr/bin/env Rscript

# ============================================================
# Gene Expression Detection Statistics
#
# Description:
#   Calculates expression statistics for a specific target gene across cell types.
#   - Loads a Seurat object.
#   - Aggregates cell sub-types (removes trailing digits, e.g., "type_1" -> "type").
#   - Computes:
#       * Global zero fraction (sparsity).
#       * Per-group statistics: N cells, Non-zero count, Zero fraction, Median, 90th percentile.
#   - Exports a CSV file with these statistics for debugging or analysis.
#
# Input:
#   - Seurat RDS file.
#
# Output:
#   - CSV file containing group-wise expression statistics.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# ==============================================================================
# Configuration
# ==============================================================================
# Input/Output
RDS_PATH      <- "path/to/your/scRNA_data.rds"
OUT_DIR       <- dirname(RDS_PATH)

# Target Gene
GENE          <- "AT5G65300"

# Analysis Settings
CELLTYPE_COL  <- "final_annotation"  # Original metadata column
ASSAY_USE     <- "RNA"
SLOT_USE      <- "data"              # Normalized data slot

# Column name for aggregated cell types
CELLTYPE_SIMPLE_COL <- "celltype_simple"

# ==============================================================================
# 1. Load Data
# ==============================================================================
x <- readRDS(RDS_PATH)
stopifnot(inherits(x, "Seurat"))

message(sprintf("Loaded Seurat: %d cells x %d genes | Assays: %s | Default: %s",
                ncol(x), nrow(x), paste(Assays(x), collapse=", "), DefaultAssay(x)))

# ==============================================================================
# 2. Validation & Normalization
# ==============================================================================
if (!(ASSAY_USE %in% Assays(x))) stop("Assay not found: ", ASSAY_USE)
DefaultAssay(x) <- ASSAY_USE

# Check if data slot exists and is populated
has_data_slot <- tryCatch({
  m <- GetAssayData(x, assay = ASSAY_USE, slot = "data")
  !is.null(m) && nrow(m) > 0 && ncol(m) > 0
}, error = function(e) FALSE)

# Run LogNormalize if necessary
if (!has_data_slot && SLOT_USE == "data") {
  message("No 'data' slot found. Running NormalizeData (LogNormalize)...")
  x <- NormalizeData(x, normalization.method = "LogNormalize", 
                     scale.factor = 10000, verbose = FALSE)
}

# ==============================================================================
# 3. Aggregate Cell Types
# ==============================================================================
# Strategy: Remove trailing digits/underscores (e.g., "cortex_1" -> "cortex")
stopifnot(CELLTYPE_COL %in% colnames(x@meta.data))
raw_ct <- as.character(x@meta.data[[CELLTYPE_COL]])
x[[CELLTYPE_SIMPLE_COL]] <- sub("_[0-9]+$", "", raw_ct)

# Log aggregation results
tab_raw <- sort(table(raw_ct), decreasing = TRUE)
tab_smp <- sort(table(x@meta.data[[CELLTYPE_SIMPLE_COL]]), decreasing = TRUE)

message(sprintf("Raw groups: %d | Collapsed groups: %d", length(tab_raw), length(tab_smp)))
message("Top collapsed groups (n cells):")
print(head(tab_smp, 20))

Idents(x) <- factor(x@meta.data[[CELLTYPE_SIMPLE_COL]])

# ==============================================================================
# 4. Calculate Statistics
# ==============================================================================
mat <- GetAssayData(x, assay = ASSAY_USE, layer = SLOT_USE)

if (!(GENE %in% rownames(mat))) {
  stop(sprintf("Gene '%s' not found in assay '%s' (slot '%s').", GENE, ASSAY_USE, SLOT_USE))
}

# Extract expression vector
v <- as.numeric(mat[GENE, ])
names(v) <- colnames(mat)

# Global Statistics
zero_frac <- mean(v == 0)
qs <- quantile(v, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.99, 1), na.rm = TRUE)

message(sprintf("[%s] Global Zero Fraction = %.3f (Detection Rate = %.3f)", 
                GENE, zero_frac, 1 - zero_frac))
message("Quantiles of expression:")
print(qs)

# Group-wise Statistics
df_debug <- data.frame(
  cell = colnames(x),
  expr = v,
  grp  = as.character(Idents(x)),
  stringsAsFactors = FALSE
)

# Calculate metrics per group
agg_res <- aggregate(expr ~ grp, df_debug, function(z) c(
  n         = length(z),
  nonzero   = sum(z > 0),
  zero_frac = mean(z == 0),
  median    = median(z),
  p90       = as.numeric(quantile(z, 0.9, na.rm = TRUE))
))

# Flatten the aggregate result into a clean data frame
stats_df <- data.frame(
  Group         = agg_res$grp,
  N_Cells       = agg_res$expr[, "n"],
  N_NonZero     = agg_res$expr[, "nonzero"],
  Zero_Fraction = agg_res$expr[, "zero_frac"],
  Median_Expr   = agg_res$expr[, "median"],
  P90_Expr      = agg_res$expr[, "p90"],
  row.names     = NULL
)

# Sort by number of detected cells
stats_df <- stats_df[order(stats_df$N_NonZero, decreasing = TRUE), ]

# ==============================================================================
# 5. Export Results
# ==============================================================================
out_file <- file.path(OUT_DIR, paste0(GENE, "_debug_group_stats.csv"))
write.csv(stats_df, out_file, row.names = FALSE)

message("Statistics saved to: ", normalizePath(out_file))