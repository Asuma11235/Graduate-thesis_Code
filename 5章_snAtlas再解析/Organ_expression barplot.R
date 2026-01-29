#!/usr/bin/env Rscript

# ============================================================
# Pseudo-bulk CPM Calculator & Visualizer
#
# Description:
#   Calculates "Pseudo-bulk" CPM (Counts Per Million) for a specific target gene
#   across multiple Seurat RDS datasets.
#   - Aggregates counts across all cells in each dataset.
#   - Normalizes by total library size.
#   - handles Seurat v4/v5 assay structure differences.
#   - Provides a fallback approximation if raw counts are missing.
#
# Output:
#   - Bar plot (PDF) comparing CPM across datasets.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

# =========================
# 1) Configuration
# =========================

# Input datasets (Label = File Path)
# NOTE: Update paths to match your directory structure.
INPUTS <- list(
  leaf   = "./input_data/Leaf_scRNA.RDS",
  stem   = "./input_data/Stem_scRNA.RDS",
  root   = "./input_data/Root_scRNA.RDS",
  flower = "./input_data/Flower_scRNA.RDS"
)

# Target Gene ID (AGI code or Symbol)
TARGET_GENE <- "AT1G63720"

# Output Settings
OUTDIR <- "./results_pseudobulk"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
OUTPDF <- file.path(OUTDIR, paste0(TARGET_GENE, "_pseudoBulk_CPM_bar.pdf"))

# =========================
# 2) Helper Functions
# =========================

# Null coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Safe nnzero function
nnz <- function(m) tryCatch(Matrix::nnzero(m), error = function(e) NA_integer_)

# Wrapper for Seurat v4/v5 compatibility (GetAssayData)
get_assay_data_safe <- function(x, assay = NULL, layer = c("counts", "data", "scale.data")) {
  layer <- match.arg(layer)
  
  # Try v5 syntax (layer) first
  ok <- try({
    return(Seurat::GetAssayData(x, assay = assay %||% Seurat::DefaultAssay(x), layer = layer))
  }, silent = TRUE)
  
  if (!inherits(ok, "try-error")) return(ok)
  
  # Fallback to v4 syntax (slot)
  slot <- switch(layer, counts = "counts", data = "data", "scale.data" = "scale.data")
  Seurat::GetAssayData(x, assay = assay %||% Seurat::DefaultAssay(x), slot = slot)
}

# Automatically select the best assay containing data
choose_best_assay <- function(x, prefer = c("RNA", "SCT", "integrated", "Spatial")) {
  asy <- Seurat::Assays(x)
  # Prioritize preferred assays, then check others
  ord <- unique(c(intersect(prefer, asy), setdiff(asy, prefer)))
  
  for (a in ord) {
    cnt <- tryCatch(get_assay_data_safe(x, assay = a, layer = "counts"), error = function(e) NULL)
    dat <- tryCatch(get_assay_data_safe(x, assay = a, layer = "data"),   error = function(e) NULL)
    
    if (!is.null(cnt) && !is.na(nnz(cnt)) && nnz(cnt) > 0) return(a)
    if (!is.null(dat) && !is.na(nnz(dat)) && nnz(dat) > 0) return(a)
  }
  ord[1]
}

# =========================
# 3) Pseudo-bulk Calculation
# =========================

calc_pseudobulk_cpm <- function(obj, gene_id, prefer_assay = c("RNA", "SCT", "integrated", "Spatial")) {
  stopifnot(inherits(obj, "Seurat"))
  Seurat::DefaultAssay(obj) <- choose_best_assay(obj, prefer = prefer_assay)
  
  # Identify gene row (case-insensitive matching)
  rn <- rownames(obj)
  gid <- toupper(gene_id)
  if (!(gid %in% toupper(rn))) {
    stop(sprintf("Target gene not found: %s", gene_id))
  }
  
  # Retrieve actual row name
  gene_rowname <- rn[match(gid, toupper(rn))]
  
  # Strategy 1: Use raw counts (True pseudo-bulk)
  cnt <- get_assay_data_safe(obj, layer = "counts")
  if (!is.null(cnt) && !is.na(nnz(cnt)) && nnz(cnt) > 0) {
    lib_size <- sum(cnt)
    if (lib_size <= 0) stop("Counts retrieved but total sum is 0.")
    
    gene_sum <- sum(cnt[gene_rowname, , drop = TRUE])
    cpm <- 1e6 * as.numeric(gene_sum) / as.numeric(lib_size)
    return(list(cpm = cpm, source = "counts"))
  }
  
  # Strategy 2: Fallback to 'data' slot (Approximation)
  warning("[fallback] Counts slot empty. approximating CPM from 'data' slot (using expm1).")
  dat <- get_assay_data_safe(obj, layer = "data")
  if (is.null(dat) || is.na(nnz(dat)) || nnz(dat) == 0) {
    stop("Both counts and data slots are empty.")
  }
  
  # Approximation: Sum of expm1(data)
  gene_sum <- sum(expm1(dat[gene_rowname, , drop = TRUE]))
  lib_size <- sum(expm1(dat))
  
  if (lib_size <= 0) stop("Sum of expm1(data) is 0.")
  
  cpm <- 1e6 * as.numeric(gene_sum) / as.numeric(lib_size)
  return(list(cpm = cpm, source = "data(approx)"))
}

# =========================
# 4) Main Execution
# =========================

message("Starting Pseudo-bulk CPM calculation...")

vals <- lapply(names(INPUTS), function(series) {
  p <- INPUTS[[series]]
  message("Processing: ", series, " (", basename(p), ")")
  
  obj <- readRDS(p)
  res <- calc_pseudobulk_cpm(obj, TARGET_GENE)
  
  data.frame(
    Series = series,
    Gene = TARGET_GENE,
    CPM = res$cpm,
    Source = res$source,
    stringsAsFactors = FALSE
  )
})

df <- do.call(rbind, vals)

# Optional: To use Log scale, uncomment below
# df$log1pCPM <- log1p(df$CPM)

# =========================
# 5) Plotting & Saving
# =========================

p <- ggplot(df, aes(x = Series, y = CPM, fill = Series)) +
  geom_col(width = 0.7) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 1)
  ) +
  labs(
    title = sprintf("Pseudo-bulk CPM of %s", TARGET_GENE),
    x = "Dataset",
    y = "CPM (counts per million)",
    caption = "Counts summed across all cells per dataset, then library-size normalized."
  )

ggsave(OUTPDF, p, width = 6.5, height = 4.5, units = "in", limitsize = FALSE)

message("\n=== Pseudo-bulk CPM Results ===")
print(df)
message("\nSaved PDF: ", normalizePath(OUTPDF))