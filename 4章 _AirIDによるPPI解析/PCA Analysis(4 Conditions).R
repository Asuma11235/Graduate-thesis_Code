#!/usr/bin/env Rscript

# ============================================================
# PCA Analysis of Protein Membership (4 Conditions)
#
# Description:
#   Performs Principal Component Analysis (PCA) based on protein presence/absence.
#   - Loads protein accession lists from Excel sheets.
#   - Constructs a binary membership matrix (Protein x Condition).
#   - Filters out zero-variance features (proteins present in all or none).
#   - Visualizes samples (mock, DMSO, bikinin, BAK1) in 2D PCA space.
#
# Color Scheme:
#   BAK1    : Red
#   bikinin : Blue-purple
#   DMSO    : Light blue
#   mock    : Yellow-green
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(tidyr)
})

# ==============================================================================
# Configuration
# ==============================================================================
# Input/Output paths
IN_XLSX_AIRID <- "path/to/your/AirID_data.xlsx"
OUT_PCA_PNG   <- "path/to/output/PCA_BR.png"

# Data Structure
SHEETS_ALL    <- c("mock", "DMSO", "bikinin", "BAK1")
PROT_COL      <- "prot_acc"

# Color Definitions
COLS_SYSTEM <- c(
  BAK1    = "#E53935", # Red
  bikinin = "#5E35B1", # Blue-purple
  DMSO    = "#4FC3F7", # Light blue
  mock    = "#9CCC65"  # Yellow-green
)

# ==============================================================================
# Helper Functions
# ==============================================================================
read_sheet_df <- function(xlsx, sheet) {
  df <- read_excel(xlsx, sheet = sheet)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  names(df) <- str_trim(names(df))
  df
}

norm_ids <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x[x == ""] <- NA_character_
  x
}

read_prot_set <- function(xlsx, sheet, prot_col) {
  df <- read_sheet_df(xlsx, sheet)
  if (!prot_col %in% names(df)) {
    stop(sprintf("Column '%s' not found in sheet '%s' (%s)", prot_col, sheet, xlsx))
  }
  v <- norm_ids(df[[prot_col]])
  v <- v[!is.na(v)]
  unique(v)
}

# ==============================================================================
# Main Process
# ==============================================================================

# 1) Load Protein Sets
sets_4 <- setNames(
  lapply(SHEETS_ALL, function(sh) read_prot_set(IN_XLSX_AIRID, sh, PROT_COL)),
  SHEETS_ALL
)

message("=== Set Sizes ===")
print(sapply(sets_4, length))

# 2) Build Membership Matrix (Binary: 0/1)
all_ids <- sort(unique(unlist(sets_4)))
mat <- sapply(sets_4, function(s) as.integer(all_ids %in% s))
rownames(mat) <- all_ids

# 3) PCA Pre-processing
# Transpose so that rows = Samples (Conditions), Cols = Features (Proteins)
X <- t(mat) 

# IMPORTANT: Filter out constant features (variance = 0)
# PCA cannot scale variables with zero variance.
vars <- apply(X, 2, stats::var)
keep <- which(!is.na(vars) & vars > 0)

message("\n=== PCA Feature Filter ===")
message("Total features (prot_acc)      : ", ncol(X))
message("Non-constant features used     : ", length(keep))
message("Dropped constant features      : ", ncol(X) - length(keep))

if (length(keep) < 2) {
  stop("Fewer than 2 non-constant features remaining. PCA cannot proceed.")
}

X_pca <- X[, keep, drop = FALSE]

# 4) Perform PCA
pca <- stats::prcomp(X_pca, center = TRUE, scale. = TRUE)
explained <- (pca$sdev^2) / sum(pca$sdev^2)

# Prepare data for plotting
scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
scores$system <- rownames(scores)

# 5) Validation
# Check if color mapping covers all systems
missing_cols <- setdiff(scores$system, names(COLS_SYSTEM))
if (length(missing_cols) > 0) {
  stop("Missing color definitions for: ", paste(missing_cols, collapse = ", "))
}

# 6) Plot PCA
p_pca <- ggplot(scores, aes(x = PC1, y = PC2, color = system, label = system)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 3.5) +
  scale_color_manual(values = COLS_SYSTEM) +
  labs(
    title = "PCA of 4 conditions (based on prot_acc membership)",
    subtitle = sprintf("PC1: %.1f%% | PC2: %.1f%% variance explained", 
                       explained[1] * 100, explained[2] * 100),
    x = sprintf("PC1 (%.1f%%)", explained[1] * 100),
    y = sprintf("PC2 (%.1f%%)", explained[2] * 100),
    color = "Condition"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

# 7) Save Output
dir.create(dirname(OUT_PCA_PNG), showWarnings = FALSE, recursive = TRUE)
ggsave(OUT_PCA_PNG, p_pca, width = 9.5, height = 6.5, dpi = 200)

message("Saved PCA Plot: ", normalizePath(OUT_PCA_PNG))
message("Done.")