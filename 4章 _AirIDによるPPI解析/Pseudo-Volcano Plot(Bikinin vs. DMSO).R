#!/usr/bin/env Rscript

# ============================================================
# Pseudo-Volcano Plot Generator
#
# Description:
#   Generates a scatter plot comparing Log2 Fold Change vs. Scaled Abundance.
#   - X-axis: Log2 Fold Change (Log2FC).
#   - Y-axis: Scaled Abundance (Conditional selection).
#       * If Log2FC > 0: Uses abundance from Sample A (e.g., F8/Bikinin).
#       * If Log2FC < 0: Uses abundance from Sample B (e.g., F7/DMSO).
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
SHEET_IDX  <- 1  # Sheet index or name (e.g., "Sheet1")

# Column Names (Must match Excel headers exactly)
COL_LOG2FC <- "Log2FC"
COL_F8     <- "Abundances (Scaled): F8: Sample, BIL7-Bikinin" # Sample corresponding to positive FC
COL_F7     <- "Abundances (Scaled): F7: Sample, BIL7-DMSO"    # Sample corresponding to negative FC

# Thresholds for Hit Definition
THR_FC     <- 1    # Absolute Log2FC threshold
THR_AB     <- 100  # Scaled Abundance threshold

# Output Paths (Derived from input filename)
OUT_PNG    <- sub("\\.xlsx$", "_pseudo_volcano.png", IN_XLSX)
OUT_HITS   <- sub("\\.xlsx$", "_hits_by_series.xlsx", IN_XLSX)

# ==============================================================================
# 1. Load Data
# ==============================================================================
if (!file.exists(IN_XLSX)) stop("Input file not found: ", IN_XLSX)

df <- read_excel(IN_XLSX, sheet = SHEET_IDX)

# Validate Column Existence
need_cols <- c(COL_LOG2FC, COL_F8, COL_F7)
miss_cols <- setdiff(need_cols, colnames(df))
if (length(miss_cols) > 0) {
  stop("Missing columns in Excel:\n  - ", paste(miss_cols, collapse = "\n  - "), call. = FALSE)
}

# ==============================================================================
# 2. Process Data (Build Pseudo-Volcano Metrics)
# ==============================================================================
plot_df <- df %>%
  mutate(
    # Ensure numeric types
    Log2FC_num = suppressWarnings(as.numeric(.data[[COL_LOG2FC]])),
    F8_num     = suppressWarnings(as.numeric(.data[[COL_F8]])),
    F7_num     = suppressWarnings(as.numeric(.data[[COL_F7]])),
    
    # Determine Series (Direction) based on Log2FC sign
    Series = case_when(
      Log2FC_num > 0 ~ "Bikinin (F8 used)",
      Log2FC_num < 0 ~ "DMSO (F7 used)",
      TRUE           ~ NA_character_
    ),
    
    # Select Abundance based on direction
    Abund_sel  = case_when(
      Log2FC_num > 0 ~ F8_num,
      Log2FC_num < 0 ~ F7_num,
      TRUE           ~ NA_real_
    ),
    
    # Define Hits (Pass both FC and Abundance thresholds)
    Hit = (abs(Log2FC_num) >= THR_FC) & (Abund_sel >= THR_AB),
    
    # Classify for Plotting
    PlotClass = case_when(
      Hit & Series == "Bikinin (F8 used)" ~ "Hit: Bikinin",
      Hit & Series == "DMSO (F7 used)"    ~ "Hit: DMSO",
      TRUE                                ~ "Not hit"
    )
  ) %>%
  # Remove invalid rows
  filter(!is.na(Log2FC_num), !is.na(Abund_sel), !is.na(Series))

# ==============================================================================
# 3. Generate Plot
# ==============================================================================
p <- ggplot(plot_df, aes(x = Log2FC_num, y = Abund_sel, color = PlotClass)) +
  geom_point(alpha = 0.7, size = 1.6) +
  
  # Add Threshold Lines
  geom_vline(xintercept = c(-THR_FC, THR_FC), linetype = "dashed") +
  geom_hline(yintercept = THR_AB, linetype = "dashed") +
  
  # Color Settings
  scale_color_manual(
    values = c(
      "Not hit"       = "grey70",
      "Hit: Bikinin"  = "red3",
      "Hit: DMSO"     = "blue3"
    )
  ) +
  
  # Labels
  labs(
    title = "Pseudo volcano plot",
    subtitle = "x = Log2FC, y = Scaled Abundance (selected by sign)",
    x = "Log2FC (Bikinin / DMSO)",
    y = "Scaled Abundance (Log2FC>0: F8, Log2FC<0: F7)",
    color = "Category"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

# Save Plot
dir.create(dirname(OUT_PNG), showWarnings = FALSE, recursive = TRUE)
ggsave(OUT_PNG, p, width = 8.5, height = 6.5, dpi = 200)
message("Saved plot: ", normalizePath(OUT_PNG))

# ==============================================================================
# 4. Export Hit Lists
# ==============================================================================
# Extract hits for Bikinin (Positive FC)
hits_bikinin <- plot_df %>%
  filter(Hit, Series == "Bikinin (F8 used)") %>%
  mutate(
    Selected_Abundance = Abund_sel,
    Selected_Column = COL_F8
  )

# Extract hits for DMSO (Negative FC)
hits_dmso <- plot_df %>%
  filter(Hit, Series == "DMSO (F7 used)") %>%
  mutate(
    Selected_Abundance = Abund_sel,
    Selected_Column = COL_F7
  )

# Save to Excel
write_xlsx(
  list(
    Bikinin_hits = hits_bikinin,
    DMSO_hits    = hits_dmso
  ),
  path = OUT_HITS
)

message("Saved hits: ", normalizePath(OUT_HITS))