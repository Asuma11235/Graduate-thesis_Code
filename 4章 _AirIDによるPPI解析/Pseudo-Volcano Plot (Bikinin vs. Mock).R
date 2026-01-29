#!/usr/bin/env Rscript

# ============================================================
# Pseudo-Volcano Plot Generator (Bikinin vs. Mock)
#
# Description:
#   Generates a scatter plot comparing Log2 Fold Change vs. Scaled Abundance.
#   - X-axis: Log2 Fold Change (Bikinin vs. Mock).
#   - Y-axis: Scaled Abundance (Conditional selection).
#       * If Log2FC > 0: Uses abundance from Bikinin sample.
#       * If Log2FC < 0: Uses abundance from Mock sample.
#   - Highlights significant hits based on both FC and Abundance thresholds.
#
# Input:
#   - Excel file containing Log2FC and Abundance columns.
#
# Output:
#   - Scatter plot (PNG).
#   - Excel file containing lists of hit proteins.
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(writexl)
})

# ==============================================================================
# Configuration
# ==============================================================================
# Input File Path
IN_XLSX    <- "path/to/your/proteomics_data.xlsx"
SHEET_IDX  <- 1

# Column Names (Must match Excel headers exactly)
COL_LOG2FC <- "Log2FC_bikinin"
COL_POS    <- "Abundances (Scaled): F8: Sample, BIL7-Bikinin" # Sample for positive FC
COL_NEG    <- "Abundances (Scaled): F6: Sample, BIL7-Mock"    # Sample for negative FC

# Thresholds
THR_FC     <- 1    # Absolute Log2FC threshold
THR_AB     <- 100  # Scaled Abundance threshold

# Output Paths (Derived from input filename)
OUT_PNG    <- sub("\\.xlsx$", "_pseudo_volcano_bikinin_mock.png", IN_XLSX)
OUT_HITS   <- sub("\\.xlsx$", "_hits_bikinin_vs_mock.xlsx", IN_XLSX)

# ==============================================================================
# 1. Load Data
# ==============================================================================
if (!file.exists(IN_XLSX)) stop("Input file not found: ", IN_XLSX)

df <- read_excel(IN_XLSX, sheet = SHEET_IDX)

# Validate Column Existence
need_cols <- c(COL_LOG2FC, COL_POS, COL_NEG)
miss_cols <- setdiff(need_cols, colnames(df))
if (length(miss_cols) > 0) {
  stop("Missing columns in Excel:\n  - ", paste(miss_cols, collapse = "\n  - "), call. = FALSE)
}

# ==============================================================================
# 2. Process Data
# ==============================================================================
plot_df <- df %>%
  mutate(
    # Ensure numeric types
    Log2FC = suppressWarnings(as.numeric(.data[[COL_LOG2FC]])),
    Pos_ab = suppressWarnings(as.numeric(.data[[COL_POS]])),
    Neg_ab = suppressWarnings(as.numeric(.data[[COL_NEG]])),
    
    # Determine Series based on direction
    Series = case_when(
      Log2FC > 0 ~ "Bikinin",
      Log2FC < 0 ~ "Mock",
      TRUE       ~ NA_character_
    ),
    
    # Select Abundance based on direction
    Abund_sel = case_when(
      Log2FC > 0 ~ Pos_ab,
      Log2FC < 0 ~ Neg_ab,
      TRUE       ~ NA_real_
    ),
    
    # Define Hits
    Hit = abs(Log2FC) >= THR_FC & Abund_sel >= THR_AB,
    
    # Classify for Plotting
    PlotClass = case_when(
      Hit & Series == "Bikinin" ~ "Hit: Bikinin",
      Hit & Series == "Mock"    ~ "Hit: Mock",
      TRUE                      ~ "Not hit"
    )
  ) %>%
  # Filter invalid rows
  filter(!is.na(Log2FC), !is.na(Abund_sel))

# ==============================================================================
# 3. Generate Plot
# ==============================================================================
p <- ggplot(plot_df, aes(x = Log2FC, y = Abund_sel, color = PlotClass)) +
  geom_point(alpha = 0.7, size = 1.6) +
  
  # Add Threshold Lines
  geom_vline(xintercept = c(-THR_FC, THR_FC), linetype = "dashed") +
  geom_hline(yintercept = THR_AB, linetype = "dashed") +
  
  # Color Settings
  scale_color_manual(values = c(
    "Not hit"      = "grey70",
    "Hit: Bikinin" = "red3",
    "Hit: Mock"    = "blue3"
  )) +
  
  # Labels
  labs(
    title = "Pseudo volcano (Log2FC_bikinin vs scaled abundance)",
    x = "Log2FC_bikinin",
    y = "Scaled Abundance ( + : Bikinin / - : Mock )",
    color = "Category"
  ) +
  theme_bw(base_size = 12)

# Save Plot
dir.create(dirname(OUT_PNG), showWarnings = FALSE, recursive = TRUE)
ggsave(OUT_PNG, p, width = 8.5, height = 6.5, dpi = 200)
message("Saved plot: ", normalizePath(OUT_PNG))

# ==============================================================================
# 4. Export Hit Lists
# ==============================================================================
# Extract hits for Bikinin (Positive FC)
hits_bikinin <- plot_df %>%
  filter(Hit, Series == "Bikinin") %>%
  mutate(Selected_Column = COL_POS)

# Extract hits for Mock (Negative FC)
hits_mock <- plot_df %>%
  filter(Hit, Series == "Mock") %>%
  mutate(Selected_Column = COL_NEG)

# Save to Excel
write_xlsx(
  list(
    Bikinin_hits = hits_bikinin,
    Mock_hits    = hits_mock
  ),
  path = OUT_HITS
)

message("Saved hits: ", normalizePath(OUT_HITS))