#!/usr/bin/env Rscript

# ============================================================
# PCA & Hierarchical Clustering Analysis (Sample-wise)
#
# Description:
#   Performs sample-level quality control and visualization using protein abundance data.
#   - Parses sample metadata (System: N/C, Treatment: DMSO/PG) from column names.
#   - Conducts PCA (Principal Component Analysis).
#   - Conducts Hierarchical Clustering based on Pearson correlation distance.
#
# Input:
#   - Excel file containing "Abundances (Scaled)" columns.
#
# Output:
#   - PCA Plot (PNG)
#   - Dendrogram (PNG)
#   - Correlation Heatmap (PNG)
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(stringr)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})

# ==============================================================================
# Configuration
# ==============================================================================
# Input/Output paths
IN_XLSX         <- "path/to/your/AmPR2_data.xlsx"
OUT_PCA_PNG     <- "path/to/output/PCA_AirID_scaledAbundance.png"
OUT_DENDRO_PNG  <- "path/to/output/HC_samples_dendrogram.png"
OUT_HEATMAP_PNG <- "path/to/output/HC_samples_heatmap.png"

# Excel Structure
SHEET_IDX       <- 1
PROT_COL        <- "Accession"

# Sample Annotation & Color Definitions
# System: N or C
# Treatment: DMSO or PG
COLS_CONDITION <- c(
  "N-DMSO" = "#90CAF9", # Light Blue
  "N-PG"   = "#1565C0", # Dark Blue
  "C-DMSO" = "#EF9A9A", # Light Red
  "C-PG"   = "#B71C1C"  # Dark Red
)

stopifnot(file.exists(IN_XLSX))

# ==============================================================================
# 1. Load Data
# ==============================================================================
df <- read_excel(IN_XLSX, sheet = SHEET_IDX) %>%
  as.data.frame(stringsAsFactors = FALSE)
names(df) <- str_trim(names(df))

if (!PROT_COL %in% names(df)) stop("Column not found: ", PROT_COL)

# Identify Abundance Columns
abund_cols <- grep("^Abundances \\(Scaled\\):", names(df), value = TRUE)

message("=== Column Check ===")
message("Protein ID Column : ", PROT_COL)
message("Abundance Columns : ", length(abund_cols))
# print(abund_cols) # Uncomment to debug

if (length(abund_cols) < 4) stop("Too few 'Abundances (Scaled)' columns found.")

# ==============================================================================
# 2. Parse Sample Metadata
# ==============================================================================
sample_info <- tibble(sample_col = abund_cols) %>%
  mutate(
    system = case_when(
      grepl("Sample,\\s*N-", sample_col) ~ "N",
      grepl("Sample,\\s*C-", sample_col) ~ "C",
      TRUE ~ NA_character_
    ),
    treatment = case_when(
      grepl("PG", sample_col)   ~ "PG",
      grepl("DMSO", sample_col) ~ "DMSO",
      TRUE ~ NA_character_
    ),
    condition = paste(system, treatment, sep = "-")
  )

# Validation
if (any(is.na(sample_info$system)) || any(is.na(sample_info$treatment))) {
  stop("Failed to parse system/treatment from column names:\n",
       paste(abund_cols, collapse = "\n"))
}

# Color Validation
missing_def <- setdiff(unique(sample_info$condition), names(COLS_CONDITION))
if (length(missing_def) > 0) {
  stop("Color not defined for conditions: ", paste(missing_def, collapse = ", "))
}

# ==============================================================================
# 3. Data Cleaning & Matrix Construction
# ==============================================================================
# Safely convert abundance columns to numeric (handling strings/commas)
df2 <- df %>%
  select(all_of(PROT_COL), all_of(abund_cols)) %>%
  mutate(
    prot_acc = as.character(.data[[PROT_COL]]),
    prot_acc = str_trim(prot_acc),
    prot_acc = na_if(prot_acc, "")
  ) %>%
  filter(!is.na(prot_acc)) %>%
  distinct(prot_acc, .keep_all = TRUE) %>%
  mutate(
    across(
      all_of(abund_cols),
      ~ {
        x <- as.character(.x)
        x <- str_trim(x)
        x[x == ""] <- NA_character_
        x <- gsub(",", "", x, fixed = TRUE) # Handle '1,000' format
        suppressWarnings(as.numeric(x))
      }
    )
  )

# Matrix: Rows = Proteins, Cols = Samples
mat_ps <- as.matrix(df2[, abund_cols, drop = FALSE])
rownames(mat_ps) <- df2$prot_acc
colnames(mat_ps) <- abund_cols # Explicitly set

# Transpose for PCA (Rows = Samples, Cols = Proteins)
X <- t(mat_ps)

message("=== Matrix Dimensions ===")
message("Samples x Proteins: ", nrow(X), " x ", ncol(X))

# ==============================================================================
# 4. Filter NA and Constant Features
# ==============================================================================
# Remove proteins with ANY missing values (to ensure robust PCA/Correlation)
na_rate_prot <- colMeans(is.na(X))
keep_prot <- which(na_rate_prot == 0)

message("\n=== Filtering ===")
message("Total proteins         : ", ncol(X))
message("Proteins without NA    : ", length(keep_prot))
message("Dropped (contain NA)   : ", ncol(X) - length(keep_prot))

if (length(keep_prot) < 2) stop("Too few proteins remain after removing NA-containing columns.")

X <- X[, keep_prot, drop = FALSE]

# Remove constant proteins (Variance = 0)
vars <- apply(X, 2, var)
keep_var <- which(!is.na(vars) & vars > 0)

message("Non-constant proteins  : ", length(keep_var))
message("Dropped constant proteins: ", ncol(X) - length(keep_var))

if (length(keep_var) < 2) stop("Too few variable proteins remain for PCA.")

X_pca <- X[, keep_var, drop = FALSE]

# ==============================================================================
# 5. PCA Analysis
# ==============================================================================
pca <- prcomp(X_pca, center = TRUE, scale. = TRUE)
explained <- (pca$sdev^2) / sum(pca$sdev^2)

scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
scores$sample_col <- rownames(scores)
scores <- left_join(scores, sample_info, by = "sample_col")

# Shorten labels for plotting (e.g., extract "F123")
scores$label <- str_extract(scores$sample_col, "F\\d+")

# Plot PCA
p_pca <- ggplot(scores, aes(x = PC1, y = PC2, color = condition, shape = system)) +
  geom_point(size = 4) +
  geom_text(aes(label = label), vjust = -0.8, size = 3) +
  scale_color_manual(values = COLS_CONDITION) +
  labs(
    title = "PCA of Samples (Scaled Abundances)",
    subtitle = sprintf("PC1: %.1f%% | PC2: %.1f%%", explained[1] * 100, explained[2] * 100),
    color = "Condition",
    shape = "System"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

# Save PCA
dir.create(dirname(OUT_PCA_PNG), showWarnings = FALSE, recursive = TRUE)
ggsave(OUT_PCA_PNG, p_pca, width = 10, height = 7, dpi = 200)
message("Saved PCA Plot: ", normalizePath(OUT_PCA_PNG))

# ==============================================================================
# 6. Hierarchical Clustering (Sample-wise)
# ==============================================================================
# Distance metric: 1 - Pearson Correlation
cor_mat <- cor(t(X_pca), method = "pearson")
dmat <- 1 - cor_mat
diag(dmat) <- 0 # Ensure self-distance is 0

hc <- hclust(as.dist(dmat), method = "average")

# Save Dendrogram
png(OUT_DENDRO_PNG, width = 1200, height = 700, res = 160)
par(mar = c(10, 4, 4, 2) + 0.1)
plot(hc, main = "HC of Samples (Distance = 1 - Pearson r)", xlab = "", sub = "", cex = 0.8)
dev.off()
message("Saved Dendrogram: ", normalizePath(OUT_DENDRO_PNG))

# ==============================================================================
# 7. Correlation Heatmap
# ==============================================================================
# Reorder correlation matrix based on clustering result
ord <- hc$labels[hc$order]
cor_df <- as.data.frame(cor_mat) %>%
  rownames_to_column("A") %>%
  pivot_longer(cols = -A, names_to = "B", values_to = "cor")

cor_df$A <- factor(cor_df$A, levels = ord)
cor_df$B <- factor(cor_df$B, levels = ord)

p_heat <- ggplot(cor_df, aes(x = B, y = A, fill = cor)) +
  geom_tile(color = "grey75", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", cor)), size = 2.6) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limit = c(-1, 1)) +
  labs(
    title = "Sample Correlation Heatmap",
    subtitle = "Ordered by hierarchical clustering",
    x = NULL, y = NULL, fill = "Pearson r"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  )

# Save Heatmap
ggsave(OUT_HEATMAP_PNG, p_heat, width = 11, height = 9, dpi = 200)
message("Saved Heatmap: ", normalizePath(OUT_HEATMAP_PNG))

message("Done.")