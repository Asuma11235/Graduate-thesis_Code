#!/usr/bin/env Rscript

# ============================================================
# PCA + Hierarchical Clustering for 8 Systems
# (AirID 4 + CoIP 4)
#
# Workflow:
#  1. Load protein accession IDs (prot_acc) from Excel sheets.
#  2. Build a binary membership matrix (prot_acc x 8 systems).
#  3. Perform PCA (prcomp) and plot PC1/PC2.
#  4. Perform Hierarchical Clustering (1 - Jaccard index) and plot Dendrogram/Heatmap.
#
# Outputs:
#  - PCA Scatter Plot (PNG)
#  - Dendrogram (PNG)
#  - Jaccard Similarity Heatmap with clustering order (PNG)
# ============================================================

# ----------------------------
# 1) Configuration
# ----------------------------

# Input Files
IN_XLSX_AIRID <- "./input_data/AirID_Dataset.xlsx"
SHEETS_AIRID  <- c("DMSO", "bikinin", "mock", "BAK1")

IN_XLSX_COIP  <- "./input_data/CoIP_Dataset.xlsx"
SHEETS_COIP   <- c("DMSO", "bikinin", "BR-up", "BR-down")

PROT_COL      <- "prot_acc"

# Output Files
OUT_DIR         <- "./results_pca_hc"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_PCA_PNG     <- file.path(OUT_DIR, "PCA_AirID4_CoIP4.png")
OUT_DENDRO_PNG  <- file.path(OUT_DIR, "HC_dendrogram_AirID4_CoIP4.png")
OUT_HEATMAP_PNG <- file.path(OUT_DIR, "HC_heatmap_Jaccard_AirID4_CoIP4.png")

# ----------------------------
# 2) Libraries
# ----------------------------
pkgs <- c("readxl", "stringr", "dplyr", "ggplot2", "tibble", "tidyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# ----------------------------
# 3) Helper Functions
# ----------------------------
read_sheet_df <- function(xlsx, sheet) {
  df <- readxl::read_excel(xlsx, sheet = sheet)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  names(df) <- stringr::str_trim(names(df))
  df
}

norm_ids <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x[x == ""] <- NA_character_
  x
}

read_prot_set <- function(xlsx, sheet, prot_col = PROT_COL) {
  df <- read_sheet_df(xlsx, sheet)
  if (!prot_col %in% names(df)) {
    stop(sprintf("Column '%s' not found in sheet '%s' (%s)", prot_col, sheet, xlsx))
  }
  v <- norm_ids(df[[prot_col]])
  v <- v[!is.na(v)]
  unique(v)
}

# ----------------------------
# 4) Data Loading & Matrix Construction
# ----------------------------
message("Loading datasets...")

airid_sets <- setNames(
  lapply(SHEETS_AIRID, function(sh) read_prot_set(IN_XLSX_AIRID, sh)),
  paste0("AirID_", SHEETS_AIRID)
)

coip_sets <- setNames(
  lapply(SHEETS_COIP, function(sh) read_prot_set(IN_XLSX_COIP, sh)),
  paste0("CoIP_", SHEETS_COIP)
)

sets_8 <- c(airid_sets, coip_sets)
set_names <- names(sets_8)

message("=== Set Sizes ===")
print(sapply(sets_8, length))

# Build membership matrix (Rows: prot_acc, Cols: 8 systems)
all_ids <- sort(unique(unlist(sets_8)))
mat <- sapply(sets_8, function(s) as.integer(all_ids %in% s))
rownames(mat) <- all_ids

# ----------------------------
# 5) PCA Analysis
# ----------------------------
# Transpose matrix so that systems are rows (observations)
pca <- stats::prcomp(t(mat), center = TRUE, scale. = TRUE)
explained <- (pca$sdev^2) / sum(pca$sdev^2)

scores <- as.data.frame(pca$x[, 1:2])
scores$system <- rownames(scores)
scores$type <- ifelse(grepl("^AirID_", scores$system), "AirID", "CoIP")

# Color scheme: Warm (AirID) vs Cool (CoIP)
cols_type <- c(AirID = "#E76F51", CoIP = "#1F77B4")

p_pca <- ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2, color = type, label = system)) +
  ggplot2::geom_point(size = 4) +
  ggplot2::geom_text(vjust = -0.8, size = 3.5, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = cols_type) +
  ggplot2::labs(
    title = "PCA of 8 Systems (Protein Membership)",
    subtitle = sprintf("PC1: %.1f%% | PC2: %.1f%% variance explained", explained[1]*100, explained[2]*100),
    color = "Method"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid = ggplot2::element_blank())

ggplot2::ggsave(OUT_PCA_PNG, p_pca, width = 8, height = 6, dpi = 300)
message("Saved PCA Plot: ", OUT_PCA_PNG)

# ----------------------------
# 6) Hierarchical Clustering (Jaccard)
# ----------------------------
jaccard_index <- function(a, b) {
  inter <- length(intersect(a, b))
  uni   <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

# Calculate pairwise Jaccard similarity matrix
jmat <- matrix(0, nrow = length(set_names), ncol = length(set_names),
               dimnames = list(set_names, set_names))

for (i in seq_along(set_names)) {
  for (j in seq_along(set_names)) {
    a <- set_names[i]; b <- set_names[j]
    jmat[a, b] <- jaccard_index(sets_8[[a]], sets_8[[b]])
  }
}

# Distance = 1 - Jaccard Similarity
dmat <- 1 - jmat
diag(dmat) <- 0

hc <- stats::hclust(stats::as.dist(dmat), method = "average")
ord <- hc$labels[hc$order]

# ----------------------------
# 7) Plot Dendrogram
# ----------------------------
png(OUT_DENDRO_PNG, width = 1200, height = 700, res = 160)
par(mar = c(8, 4, 4, 2) + 0.1)
plot(hc, main = "Hierarchical Clustering (Distance = 1 - Jaccard)",
     xlab = "", sub = "", cex = 0.9)
dev.off()
message("Saved Dendrogram: ", OUT_DENDRO_PNG)

# ----------------------------
# 8) Plot Heatmap
# ----------------------------
# Convert matrix to long format for ggplot
j_df <- as.data.frame(jmat) %>%
  tibble::rownames_to_column("A") %>%
  tidyr::pivot_longer(cols = -A, names_to = "B", values_to = "jacc")

# Order factors based on clustering result
j_df$A <- factor(j_df$A, levels = ord)
j_df$B <- factor(j_df$B, levels = ord)

p_heat <- ggplot2::ggplot(j_df, ggplot2::aes(x = B, y = A, fill = jacc)) +
  ggplot2::geom_tile(color = "grey75", linewidth = 0.35) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", jacc)), size = 3.2) +
  ggplot2::scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, 1), name = "Jaccard") +
  ggplot2::labs(
    title = "Jaccard Similarity Heatmap",
    subtitle = "Ordered by Hierarchical Clustering",
    x = NULL, y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

ggplot2::ggsave(OUT_HEATMAP_PNG, p_heat, width = 9, height = 7, dpi = 300)
message("Saved Heatmap: ", OUT_HEATMAP_PNG)

message("Analysis Complete.")