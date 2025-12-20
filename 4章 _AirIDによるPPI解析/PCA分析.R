# ============================================================
# PCA + Hierarchical Clustering for 8 systems (AirID 4 + CoIP 4)
# - Reload prot_acc from Excel sheets
# - Build 0/1 membership matrix (prot_acc × 8 systems)
# - PCA (prcomp) and plot PC1/PC2
# - Hierarchical clustering (1 - Jaccard) and dendrogram / heatmap
#
# Outputs:
#  - PCA scatter (PNG)
#  - Dendrogram (PNG)
#  - Jaccard heatmap with clustering order (PNG)
#
# Notes:
#  - Designed to be robust in Bioconductor-heavy environments:
#    use explicit namespaces (dplyr::, etc.) to avoid select() conflicts.
# ============================================================

# ----------------------------
# INPUT
# ----------------------------
IN_XLSX_AIRID <- "/Users/ryo/Downloads/251220_AirID.xlsx"
SHEETS_AIRID  <- c("DMSO","bikinin","mock","BAK1")

IN_XLSX_COIP  <- "/Users/ryo/Downloads/251108_CoIP.xlsx"
SHEETS_COIP   <- c("DMSO","bikinin","BR-up","BR-down")

PROT_COL <- "prot_acc"

# ----------------------------
# OUTPUT
# ----------------------------
OUT_PCA_PNG      <- "/Users/ryo/Downloads/PCA_AirID4_CoIP4.png"
OUT_DENDRO_PNG   <- "/Users/ryo/Downloads/HC_dendrogram_AirID4_CoIP4.png"
OUT_HEATMAP_PNG  <- "/Users/ryo/Downloads/HC_heatmap_Jaccard_AirID4_CoIP4.png"

# ----------------------------
# Packages
# ----------------------------
pkgs <- c("readxl","stringr","dplyr","ggplot2")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(readxl)
library(stringr)
library(dplyr)
library(ggplot2)

# ----------------------------
# Helpers: read prot_acc sets
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
# Build 8 sets
# ----------------------------
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

cat("=== set sizes ===\n")
print(sapply(sets_8, length))

# ----------------------------
# Membership matrix (prot_acc × 8)
# ----------------------------
all_ids <- sort(unique(unlist(sets_8)))
mat <- sapply(sets_8, function(s) as.integer(all_ids %in% s))
rownames(mat) <- all_ids

# ----------------------------
# PCA (systems as rows => transpose)
# ----------------------------
pca <- stats::prcomp(t(mat), center = TRUE, scale. = TRUE)
explained <- (pca$sdev^2) / sum(pca$sdev^2)

scores <- as.data.frame(pca$x[, 1:2])
scores$system <- rownames(scores)
scores$type <- ifelse(grepl("^AirID_", scores$system), "AirID", "CoIP")

# Color scheme: warm for AirID, cool for CoIP
cols_type <- c(AirID = "#E76F51", CoIP = "#1F77B4")

p_pca <- ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2, color = type, label = system)) +
  ggplot2::geom_point(size = 4) +
  ggplot2::geom_text(vjust = -0.8, size = 3.5) +
  ggplot2::scale_color_manual(values = cols_type) +
  ggplot2::labs(
    title = "PCA of 8 systems (based on prot_acc membership)",
    subtitle = sprintf("PC1: %.1f%% | PC2: %.1f%% variance explained", explained[1]*100, explained[2]*100),
    color = "Method"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid = ggplot2::element_blank())

ggplot2::ggsave(OUT_PCA_PNG, p_pca, width = 9.5, height = 6.5, dpi = 200)
cat("Saved PCA:", OUT_PCA_PNG, "\n")

# ----------------------------
# Hierarchical clustering using 1 - Jaccard similarity
# ----------------------------
jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  uni   <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

# build Jaccard matrix among 8 sets
jmat <- matrix(0, nrow = length(set_names), ncol = length(set_names),
               dimnames = list(set_names, set_names))

for (i in seq_along(set_names)) {
  for (j in seq_along(set_names)) {
    a <- set_names[i]; b <- set_names[j]
    jmat[a, b] <- jaccard(sets_8[[a]], sets_8[[b]])
  }
}

dmat <- 1 - jmat
diag(dmat) <- 0

hc <- stats::hclust(stats::as.dist(dmat), method = "average")
ord <- hc$labels[hc$order]

# ----------------------------
# Dendrogram plot (base)
# ----------------------------
png(OUT_DENDRO_PNG, width = 1200, height = 700, res = 160)
par(mar = c(8, 4, 4, 2) + 0.1)
plot(hc, main = "Hierarchical clustering (distance = 1 - Jaccard)", xlab = "", sub = "", cex = 0.9)
dev.off()
cat("Saved dendrogram:", OUT_DENDRO_PNG, "\n")

# ----------------------------
# Heatmap of Jaccard (ordered by clustering)
# ----------------------------
j_df <- as.data.frame(jmat) %>%
  tibble::rownames_to_column("A") %>%
  tidyr::pivot_longer(cols = -A, names_to = "B", values_to = "jacc")

j_df$A <- factor(j_df$A, levels = ord)
j_df$B <- factor(j_df$B, levels = ord)

p_heat <- ggplot2::ggplot(j_df, ggplot2::aes(x = B, y = A, fill = jacc)) +
  ggplot2::geom_tile(color = "grey75", linewidth = 0.35) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", jacc)), size = 3.2) +
  ggplot2::labs(
    title = "Jaccard similarity heatmap (ordered by hierarchical clustering)",
    x = NULL, y = NULL, fill = "Jaccard"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
  )

ggplot2::ggsave(OUT_HEATMAP_PNG, p_heat, width = 9.5, height = 7.5, dpi = 200)
cat("Saved heatmap:", OUT_HEATMAP_PNG, "\n")

cat("Done.\n")

