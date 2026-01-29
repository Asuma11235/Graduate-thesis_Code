#!/usr/bin/env Rscript

# ============================================================
# Time-Course Gene Expression Plotter (Flower & Silique)
#
# Description:
#   Visualizes gene expression changes across a composite developmental timeline
#   combining Flower stages (Early-Late) and Silique stages (D0-D6).
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
OUTPUT_DIR_NAME <- "TimeCourse_Plot_Output"
FILE_PATTERN    <- "%s_rename.slim.RDS" # e.g., "D0_silique_rename.slim.RDS"

# --- Stage Definitions ---
# 1. Flower Series
STAGES_FLOWER   <- c("Early_flower", "Middle_flower", "Late_flower")

# 2. Silique Series (D0-D6)
STAGES_SILIQUE  <- paste0("D", 0:6, "_silique")

# 3. Combined Ordered List (Plotting Order)
TARGET_STAGES   <- c(STAGES_FLOWER, STAGES_SILIQUE)

# Target Genes and Colors
TARGET_GENES    <- c("AT1G63720", "AT5G52430", "AT4G25620", "AT1G76660")

GENE_COLORS     <- c(
  "AT1G63720" = "firebrick",
  "AT5G52430" = "orange",
  "AT4G25620" = "forestgreen",
  "AT1G76660" = "dodgerblue"
)

# ==============================================================================
# 1. Data Aggregation Loop
# ==============================================================================
results_list <- list()

for (stage_name in TARGET_STAGES) {
  
  # Construct file path
  file_name <- sprintf(FILE_PATTERN, stage_name)
  rds_path  <- file.path(BASE_DIR, file_name)
  
  # Skip if file does not exist
  if (!file.exists(rds_path)) {
    warning("File not found: ", rds_path)
    next
  }
  
  # Load Data
  x <- readRDS(rds_path)
  DefaultAssay(x) <- "RNA"
  
  # Check and Perform Normalization (Compatible with Seurat v4/v5)
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
    # Ensure stages are ordered as defined in TARGET_STAGES
    Stage = factor(Stage, levels = TARGET_STAGES)
  )

# Generate Line Plot
p <- ggplot(final_df, aes(x = Stage, y = Expression, group = Gene, color = Gene)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = GENE_COLORS, drop = FALSE) +
  labs(
    title = "Gene Expression (Flower & Silique)",
    y = "Average Expression (LogNorm)",
    x = "Developmental Stage"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate labels for readability
  )

# ==============================================================================
# 3. Export Results
# ==============================================================================
out_dir <- file.path(BASE_DIR, OUTPUT_DIR_NAME)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

timestamp <- format(Sys.time(), "%H%M%S")

# --- Save as PNG ---
png_path <- file.path(out_dir, paste0("TimeCourse_FlowerSilique_", timestamp, ".png"))
ggsave(png_path, p, width = 10, height = 6, dpi = 300)
message("Saved PNG: ", normalizePath(png_path))

# --- Save as PPTX (Editable Vector Graphics) ---
pptx_path <- file.path(out_dir, paste0("TimeCourse_FlowerSilique_", timestamp, ".pptx"))

doc <- read_pptx() %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = dml(ggobj = p), location = ph_location(width = 10, height = 6))

print(doc, target = pptx_path)
message("Saved PPTX: ", normalizePath(pptx_path))