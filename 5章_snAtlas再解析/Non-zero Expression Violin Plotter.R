#!/usr/bin/env Rscript

# ============================================================
# Non-zero Expression Violin Plotter
#
# Description:
#   Generates a Violin Plot for a specific gene, using ONLY cells where
#   the gene is expressed (expression > 0).
#   - Loads a Seurat object.
#   - Aggregates cell types (removes trailing digits).
#   - Filters out cells with zero expression for the target gene.
#   - Visualizes the distribution of expression levels in the remaining cells.
#
# Input:
#   - Seurat RDS file.
#
# Output:
#   - Violin Plot (PDF & PNG).
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
GENE          <- "AT4G30190"

# Analysis Settings
CELLTYPE_COL  <- "final_annotation" # Original metadata column
ASSAY_USE     <- "RNA"
SLOT_USE      <- "data"             # Normalized data slot

# Column name for aggregated cell types
CELLTYPE_SIMPLE_COL <- "celltype_simple"

# ==============================================================================
# 1. Load Data
# ==============================================================================
x <- readRDS(RDS_PATH)
stopifnot(inherits(x, "Seurat"))

# ==============================================================================
# 2. Validation & Normalization
# ==============================================================================
if (!(ASSAY_USE %in% Assays(x))) stop("Assay not found: ", ASSAY_USE)
DefaultAssay(x) <- ASSAY_USE

# Check if data slot exists and is populated
if (SLOT_USE == "data") {
  has_data_slot <- tryCatch({
    m <- GetAssayData(x, assay = ASSAY_USE, slot = "data")
    !is.null(m) && nrow(m) > 0 && ncol(m) > 0
  }, error = function(e) FALSE)
  
  if (!has_data_slot) {
    message("No 'data' slot found. Running NormalizeData (LogNormalize)...")
    x <- NormalizeData(x, normalization.method = "LogNormalize", 
                       scale.factor = 10000, verbose = FALSE)
  }
}

# ==============================================================================
# 3. Aggregate Cell Types
# ==============================================================================
stopifnot(CELLTYPE_COL %in% colnames(x@meta.data))
raw_ct <- as.character(x@meta.data[[CELLTYPE_COL]])

# Remove trailing digits/underscores (e.g., "type_1" -> "type")
x[[CELLTYPE_SIMPLE_COL]] <- sub("_[0-9]+$", "", raw_ct)

# ==============================================================================
# 4. Filter Non-zero Cells
# ==============================================================================
mat <- GetAssayData(x, assay = ASSAY_USE, slot = SLOT_USE)

if (!(GENE %in% rownames(mat))) stop("Gene not found: ", GENE)

# Extract expression vector
v <- as.numeric(mat[GENE, ])
names(v) <- colnames(mat)

# Identify cells with non-zero expression
keep_cells <- names(v)[v > 0]
n_total <- length(v)
n_keep  <- length(keep_cells)
pct_keep <- 100 * n_keep / n_total

message(sprintf("Non-zero cells for %s: %d / %d (%.2f%%)", GENE, n_keep, n_total, pct_keep))

if (n_keep < 50) {
  warning("Low number of non-zero cells. The violin plot may be unstable: ", n_keep)
}

if (n_keep == 0) {
  stop("No cells express this gene. Cannot plot.")
}

# Subset Seurat object to keep only expressing cells
x_nz <- subset(x, cells = keep_cells)

# Set identity to the aggregated cell type
Idents(x_nz) <- factor(x_nz@meta.data[[CELLTYPE_SIMPLE_COL]])

# ==============================================================================
# 5. Generate Violin Plot
# ==============================================================================
p <- VlnPlot(
  x_nz,
  features = GENE,
  group.by = CELLTYPE_SIMPLE_COL,
  assay = ASSAY_USE,
  layer = SLOT_USE, # 'layer' is preferred over 'slot' in Seurat v5+, but 'slot' works
  pt.size = 0.2
) +
  ggtitle(paste0(GENE, " expression (Non-zero cells only)")) +
  xlab(CELLTYPE_SIMPLE_COL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# ==============================================================================
# 6. Save Outputs
# ==============================================================================
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

out_pdf <- file.path(OUT_DIR, paste0(GENE, "_VlnPlot_nonzero_by_", CELLTYPE_SIMPLE_COL, ".pdf"))
out_png <- file.path(OUT_DIR, paste0(GENE, "_VlnPlot_nonzero_by_", CELLTYPE_SIMPLE_COL, ".png"))

ggsave(out_pdf, plot = p, width = 12, height = 6, useDingbats = FALSE)
ggsave(out_png, plot = p, width = 12, height = 6, dpi = 300)

message("Saved PDF: ", normalizePath(out_pdf))
message("Saved PNG: ", normalizePath(out_png))