#!/usr/bin/env Rscript

# ============================================================
# Time-Course Gene Expression Plotter
#
# Description:
#   Visualizes gene expression changes across developmental stages (S1-S6).
#   - Loads separate Seurat RDS files for each stage.
#   - Calculates average expression for target genes per stage.
#   - Generates a line plot comparing expression trends.
#   - Exports results to PNG and editable PowerPoint (PPTX) formats.
#
# Dependencies:
#   - Seurat, tidyverse, ggplot2, officer, rvg
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(officer)
  library(rvg)
})

# ==============================================================================
# Configuration
# ==============================================================================
# Input/Output paths
BASE_DIR        <- "path/to/your/scRNA_data"
FILE_PATTERN    <- "S%d_leaf_rename.slim.RDS" # Pattern for files (e.g., %d becomes 1, 2...)
OUTPUT_DIR_NAME <- "TimeCourse_Plot_Output"

# Target Genes and Colors
TARGET_GENES <- c("AT5G65300", "AT1G35210", "AT1G72240", "AT1G22470")

GENE_COLORS <- c(
  "AT5G65300" = "firebrick",
  "AT1G35210" = "orange",
  "AT1G72240" = "forestgreen",
  "AT1G22470" = "dodgerblue"
)

# Stages to analyze
STAGES <- 1:6

# ==============================================================================
# 1. Data Aggregation Loop
# ==============================================================================
results_list <- list()

for (i in STAGES) {
  stage_name <- paste0("S", i)
  rds_path   <- file.path(BASE_DIR, sprintf(FILE_PATTERN, i))
  
  if (!file.exists(rds_path)) {
    warning(sprintf("File not found for stage %s: %s", stage_name, rds_path))
    next
  }
  
  # Load Data
  x <- readRDS(rds_path)
  DefaultAssay(x) <- "RNA"
  
  # Check and Perform Normalization if needed (Compatible with Seurat v4/v5)
  has_data <- tryCatch({
    if (is(x[["RNA"]], "Assay5")) {
      !is.null(Layers(x, search = "data"))
    } else {
      "data" %in% slotNames(x[["RNA"]])
    }
  }, error = function(e) FALSE)
  
  if (!has_data) {
    x <- NormalizeData(x, verbose = FALSE)
  }
  
  # Check gene existence
  valid_genes <- intersect(TARGET_GENES, rownames(x))
  if (length(valid_genes) == 0) {
    rm(x); gc()
    next
  }
  
  # Calculate Average Expression
  # Note: Assumes the entire object represents the stage (taking the first column)
  avg_obj <- tryCatch({
    AverageExpression(x, features = valid_genes, assays = "RNA", slot = "data", verbose = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(avg_obj)) {
    rm(x); gc()
    next
  }
  
  avg_mat <- avg_obj[["RNA"]]
  if (ncol(avg_mat) < 1) {
    rm(x); gc()
    next
  }
  
  # Format Data
  df_tmp <- tibble(
    Gene = rownames(avg_mat),
    Stage = stage_name,
    Expression = as.numeric(avg_mat[, 1]) # Use the first group/column
  ) %>%
    filter(!is.na(Expression))
  
  results_list[[stage_name]] <- df_tmp
  
  # Cleanup memory
  rm(x); gc()
}

# ==============================================================================
# 2. Plotting
# ==============================================================================
if (length(results_list) == 0) stop("No data could be aggregated.")

final_df <- bind_rows(results_list) %>%
  mutate(
    Gene  = trimws(Gene),
    Stage = factor(Stage, levels = paste0("S", STAGES))
  )

# Generate Line Plot
p <- ggplot(final_df, aes(x = Stage, y = Expression, group = Gene, color = Gene)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = GENE_COLORS, drop = FALSE) +
  labs(
    title = "Gene Expression Changes (S1–S6)",
    y = "Average Expression (LogNorm)",
    x = "Developmental Stage"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black")
  )

# ==============================================================================
# 3. Export Results
# ==============================================================================
out_dir <- file.path(BASE_DIR, OUTPUT_DIR_NAME)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

timestamp <- format(Sys.time(), "%H%M%S")

# --- Save as PNG ---
png_path <- file.path(out_dir, paste0("TimeCourse_Expression_", timestamp, ".png"))
ggsave(png_path, p, width = 8, height = 6, dpi = 300)
message("Saved PNG: ", normalizePath(png_path))

# --- Save as PPTX (Editable Vector Graphics) ---
pptx_path <- file.path(out_dir, paste0("TimeCourse_Expression_", timestamp, ".pptx"))

doc <- read_pptx() %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = dml(ggobj = p), location = ph_location(width = 8, height = 6))

print(doc, target = pptx_path)
message("Saved PPTX: ", normalizePath(pptx_path))