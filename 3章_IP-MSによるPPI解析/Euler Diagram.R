#!/usr/bin/env Rscript

# ============================================================
# Euler Diagram Generator (4 Conditions)
#
# Description:
#   Generates an Euler (Venn-like) diagram to visualize overlaps 
#   between four specific conditions: bikinin, DMSO, BR-up, and BR-down.
#   - Input: Excel file with protein accession columns.
#   - Shape: Ellipse approximation.
#   - Color Scheme: Warm colors (bikinin, BR-up) vs Cool colors (DMSO, BR-down).
#
# Output:
#   - PDF and PNG files of the Euler diagram.
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(stringr)
  library(eulerr)
})

# ==============================================================================
# Configuration
# ==============================================================================
# Input/Output paths
IN_XLSX  <- "path/to/your/CoIP_data.xlsx"
OUT_PDF  <- "path/to/output/Venn_BR_conditions.pdf"
OUT_PNG  <- "path/to/output/Venn_BR_conditions.png"

# Data Structure
# Sheet names corresponding to the 4 conditions
SHEETS   <- c("bikinin", "DMSO", "BR-up", "BR-down")
PROT_COL <- "prot_acc"  # Column name for protein IDs

# Color Definitions
# bikinin, BR-up -> Warm colors
# DMSO, BR-down  -> Cool colors
FILL_COLORS <- c(
  "bikinin" = "#E53935", # Red (Warm)
  "BR-up"   = "#FB8C00", # Orange (Warm)
  "DMSO"    = "#3366ff", # Light Blue (Cool)
  "BR-down" = "#00ACC1"  # Teal/Cyan (Cool)
)

# ==============================================================================
# Helper Function: Read Protein Set
# ==============================================================================
read_prot_set <- function(xlsx, sheet, prot_col) {
  # Check if sheet exists
  avail_sheets <- readxl::excel_sheets(xlsx)
  if (!sheet %in% avail_sheets) {
    stop(sprintf("Sheet '%s' not found in file '%s'", sheet, xlsx))
  }
  
  # Read and clean data
  df <- readxl::read_excel(xlsx, sheet = sheet)
  names(df) <- stringr::str_trim(names(df))
  
  if (!prot_col %in% names(df)) {
    stop(sprintf("Column '%s' not found in sheet '%s'", prot_col, sheet))
  }
  
  # Extract unique IDs (non-NA, non-empty)
  v <- as.character(df[[prot_col]])
  v <- stringr::str_trim(v)
  v <- v[!is.na(v) & nzchar(v)]
  unique(v)
}

# ==============================================================================
# Main Process
# ==============================================================================

# 1) Build Sets
# Read protein IDs from each sheet into a named list
sets_4 <- setNames(
  lapply(SHEETS, function(sh) read_prot_set(IN_XLSX, sh, PROT_COL)),
  SHEETS
)

message("=== Set Sizes ===")
print(sapply(sets_4, length))

# 2) Fit Euler Diagram
# Use "ellipse" shape for approximation
fit <- eulerr::euler(sets_4, shape = "ellipse")

# 3) Prepare Colors for Plotting
# Validation: Check if all conditions have a defined color
missing_cols <- setdiff(names(sets_4), names(FILL_COLORS))
if (length(missing_cols) > 0) {
  stop("Missing color definitions for: ", paste(missing_cols, collapse = ", "))
}

# Align colors with the order of labels in the fitted object
set_order <- fit$labels
if (is.null(set_order) || length(set_order) == 0) {
  set_order <- names(sets_4) # Fallback
}

fill_cols_plot <- unname(FILL_COLORS[set_order])

# Validation: Check for NA in colors
if (any(is.na(fill_cols_plot))) {
  stop("NA detected in plot colors. Please check name matching between data and color config.")
}

# 4) Plotting Function
# Defined as a function to reuse for both PDF and PNG devices
plot_euler <- function() {
  plot(
    fit,
    # Shape fill settings (semi-transparent)
    fills = list(fill = fill_cols_plot, alpha = 0.45),
    
    # Border and label settings
    edges  = list(col = "grey30", lwd = 1),
    labels = list(col = "grey10", cex = 1.0),
    quantities = list(col = "grey10", cex = 0.8),
    
    # Legend settings (alpha = 1 ensures the legend keys are opaque)
    legend = list(cex = 1.0, alpha = 1),
    
    main = ""
  )
}

# 5) Save Outputs
# Ensure output directory exists
out_dir <- dirname(OUT_PDF)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Save PDF
pdf(OUT_PDF, width = 9, height = 7)
plot_euler()
invisible(dev.off())

# Save PNG
png(OUT_PNG, width = 1800, height = 1400, res = 200)
plot_euler()
invisible(dev.off())

message("Files saved:")
message(normalizePath(OUT_PDF))
message(normalizePath(OUT_PNG))