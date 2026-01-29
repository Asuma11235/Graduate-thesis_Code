#!/usr/bin/env Rscript

# ============================================================
# Hierarchical Clustering Heatmap (Scaled Abundances)
#
# Description:
#   Generates a heatmap for specific protein abundance data.
#   - Filters columns for specified keywords (e.g., BAK1, BIL7).
#   - Imputes missing values using row-wise medians.
#   - Apply row-wise Z-score scaling.
#   - Performs Hierarchical Clustering (Euclidean distance, Complete linkage).
#
# Input:
#   - Excel file containing "Abundances (Scaled)" columns.
#     Format: "Abundances (Scaled): F?: Sample, Description"
#
# Output:
#   - PNG image of the heatmap with column annotations.
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(pheatmap)
})

# ==============================================================================
# Configuration
# ==============================================================================
# Input/Output paths
IN_XLSX  <- "path/to/your/proteomics_data.xlsx"
OUT_PNG  <- "path/to/output/HC_heatmap_BAK1_BIL7.png"

# Excel Structure
SHEET_IDX           <- 1
ID_COL              <- "Accession"
ABUND_PATTERN       <- "^Abundances \\(Scaled\\):"

# Filtering Logic (Regex for columns to keep)
KEEP_REGEX          <- "(BAK1|BIL7)"

# Clustering & Visualization Parameters
DIST_METHOD         <- "euclidean"
LINKAGE             <- "complete"
SHOW_ROW_NAMES      <- FALSE
SHOW_COL_NAMES      <- TRUE
IMPUTE_ROW_MEDIAN   <- TRUE

# Color Scheme for Groups
GROUP_COLORS <- c(
  BAK1    = "#E53935", # Red
  bikinin = "#5E35B1", # Blue-purple
  DMSO    = "#4FC3F7", # Light blue
  mock    = "#9CCC65"  # Yellow-green
)

# ==============================================================================
# Main Process
# ==============================================================================

# 1) Load Excel Data
df <- read_excel(IN_XLSX, sheet = SHEET_IDX) |>
  as.data.frame(stringsAsFactors = FALSE)

names(df) <- str_trim(names(df))

if (!ID_COL %in% names(df)) stop("ID column not found: ", ID_COL)

# 2) Identify Abundance Columns
abund_cols_all <- grep(ABUND_PATTERN, names(df), value = TRUE)
if (length(abund_cols_all) == 0) stop("No columns matched pattern: ", ABUND_PATTERN)

# 3) Filter for Target Samples (BAK1 / BIL7)
abund_cols <- abund_cols_all[str_detect(abund_cols_all, KEEP_REGEX)]

if (length(abund_cols) == 0) {
  stop("No abundance columns matched filter regex: ", KEEP_REGEX,
       "\nAvailable columns were:\n", paste(abund_cols_all, collapse = "\n"))
}

message("Detected abundance columns (ALL):")
print(abund_cols_all)
message("\nSelected columns (Filtered):")
print(abund_cols)

# Optional warning if column count deviates from expectation (e.g., 4 columns)
if (length(abund_cols) != 4) {
  warning("Number of selected columns is ", length(abund_cols), 
          " (Expected 4). Proceeding with current selection.")
}

# 4) Construct Data Matrix
mat <- df[, abund_cols, drop = FALSE]
mat <- as.matrix(sapply(mat, as.numeric))
rownames(mat) <- as.character(df[[ID_COL]])

# 5) Handle Missing Values (Imputation)
if (IMPUTE_ROW_MEDIAN) {
  for (i in seq_len(nrow(mat))) {
    if (anyNA(mat[i, ])) {
      med <- median(mat[i, ], na.rm = TRUE)
      mat[i, is.na(mat[i, ])] <- med
    }
  }
} else {
  mat[is.na(mat)] <- 0
}

# 6) Row-wise Z-score Scaling
row_z_score <- function(x) {
  s <- sd(x)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x)) / s
}
mat_scaled <- t(apply(mat, 1, row_z_score))

# 7) Extract Clean Column Names
# Removes prefix "Abundances (Scaled): F?: Sample, "
sample_names <- abund_cols |>
  str_remove("^Abundances \\(Scaled\\):\\s*") |>
  str_trim()

colnames(mat_scaled) <- sample_names

# 8) Prepare Annotations
# Define grouping based on sample name keywords
ann_group_vec <- dplyr::case_when(
  str_detect(sample_names, "BAK1")    ~ "BAK1",
  str_detect(sample_names, "Bikinin") ~ "bikinin",
  str_detect(sample_names, "DMSO")    ~ "DMSO",
  str_detect(sample_names, "Mock")    ~ "mock",
  TRUE ~ NA_character_
)

ann_col <- data.frame(
  Group = factor(ann_group_vec, levels = names(GROUP_COLORS)),
  row.names = sample_names
)

ann_colors_list <- list(
  Group = GROUP_COLORS
)

# 9) Generate and Save Heatmap
dir.create(dirname(OUT_PNG), showWarnings = FALSE, recursive = TRUE)

pheatmap(
  mat_scaled,
  color = colorRampPalette(c("green", "black", "red"))(100),
  clustering_distance_rows = DIST_METHOD,
  clustering_distance_cols = DIST_METHOD,
  clustering_method = LINKAGE,
  annotation_col = ann_col,
  annotation_colors = ann_colors_list,
  show_rownames = SHOW_ROW_NAMES,
  show_colnames = SHOW_COL_NAMES,
  fontsize_col = 10,
  border_color = NA,
  filename = OUT_PNG,
  width = 7.5,
  height = 8.5
)

message("Heatmap saved to:\n", normalizePath(OUT_PNG))