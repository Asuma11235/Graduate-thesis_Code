#!/usr/bin/env Rscript

# ============================================================
# Co-expression Analysis (Binomial & Fisher's Exact Test)
#
# Description:
#   Analyzes gene co-expression patterns in scRNA-seq data (Seurat Object).
#   1. Subset Analysis (Size 2-4):
#      - Uses Binomial test to evaluate enrichment/depletion of simultaneous expression
#      - Compares observed co-expression counts vs. expected counts under independence.
#   2. Pairwise Analysis (Size 2):
#      - Uses Fisher's Exact Test to evaluate association (Odds Ratio).
#      - Determines Co-expression (Odds > 1) or Mutual Exclusion (Odds < 1).
#
#   * Analyses are performed both Globally (all cells) and Per Cell Type.
#
# Output:
#   - Excel file containing statistical results for all subsets and pairs.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(openxlsx)
})

# ==============================================================================
# Configuration
# ==============================================================================
# Input Data
RDS_PATH        <- "path/to/your/scRNA_data.rds"
OUTDIR          <- dirname(RDS_PATH)

# Target Genes
GENES           <- c("AT1G63720", "AT5G52430", "AT4G25620", 
                     "AT1G76660", "AT1G75080", "AT4G18710")

# Analysis Settings
CELLTYPE_COL    <- "final_annotation"
ASSAY_USE       <- "RNA"
SLOT_USE        <- "data"  # Assumes log-normalized data
EXPR_THRESHOLD  <- 0       # Threshold to define "Expressed" (e.g., >0)

# Statistical Settings
P_ADJUST_METHOD <- "BH"    # "BH", "bonferroni", "holm", etc.
SUBSET_SIZES    <- 2:4     # Size of gene subsets to analyze (e.g., pairs, triplets)
MIN_CELLS_GROUP <- 50      # Minimum cells required to analyze a specific group

# ==============================================================================
# Helper Functions
# ==============================================================================

# Clean cell type names (remove trailing digits/underscores)
collapse_celltype <- function(x) sub("_[0-9]+$", "", as.character(x))

# Generate all combinations of size k
all_subsets <- function(vec, k) {
  combn(vec, k, simplify = FALSE)
}

# Binomial Test for Co-expression of a Subset S
subset_binom_test <- function(expr_bool_mat, S) {
  # expr_bool_mat: cells x genes (logical matrix)
  N <- nrow(expr_bool_mat)
  
  # Probability of expression for each gene
  ps <- colMeans(expr_bool_mat[, S, drop = FALSE])
  
  # Expected probability of co-expression under independence assumption
  p_indep <- prod(ps)
  
  # Observed count of cells expressing ALL genes in S
  k_obs <- sum(rowSums(expr_bool_mat[, S, drop = FALSE]) == length(S))
  exp_k <- N * p_indep
  
  # Enrichment: P(X >= k_obs)
  p_enrich <- pbinom(k_obs - 1, size = N, prob = p_indep, lower.tail = FALSE)
  
  # Depletion: P(X <= k_obs)
  p_deplete <- pbinom(k_obs, size = N, prob = p_indep, lower.tail = TRUE)
  
  data.frame(
    subset = paste(S, collapse = "|"),
    k = length(S),
    N = N,
    observed_all = k_obs,
    expected_all = exp_k,
    fold = ifelse(exp_k > 0, k_obs / exp_k, NA_real_),
    p_indep = p_indep,
    p_enrich = p_enrich,
    p_deplete = p_deplete,
    stringsAsFactors = FALSE
  )
}

# Fisher's Exact Test for a Pair of Genes
pair_fisher <- function(expr_bool_mat, g1, g2) {
  a <- expr_bool_mat[, g1]
  b <- expr_bool_mat[, g2]
  tab <- table(a, b) # rows: a F/T, cols: b F/T
  
  # Fisher test requires a 2x2 table
  if (!all(dim(tab) == c(2, 2))) {
    # Handle cases where a gene is all-TRUE or all-FALSE
    return(data.frame(
      g1 = g1, g2 = g2, N = nrow(expr_bool_mat),
      n11 = sum(a & b), n10 = sum(a & !b), n01 = sum(!a & b), n00 = sum(!a & !b),
      odds_ratio = NA_real_, p_greater = NA_real_, p_less = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  # Test for positive association (Co-expression)
  ft_greater <- fisher.test(tab, alternative = "greater")
  # Test for negative association (Mutual Exclusion)
  ft_less    <- fisher.test(tab, alternative = "less")
  
  data.frame(
    g1 = g1, g2 = g2, N = nrow(expr_bool_mat),
    n11 = sum(a & b), n10 = sum(a & !b), n01 = sum(!a & b), n00 = sum(!a & !b),
    odds_ratio = unname(ft_greater$estimate),
    p_greater = ft_greater$p.value,
    p_less = ft_less$p.value,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# 1. Load Data & Preprocessing
# ==============================================================================
x <- readRDS(RDS_PATH)
stopifnot(inherits(x, "Seurat"))

if (!(ASSAY_USE %in% Assays(x))) stop("Assay not found: ", ASSAY_USE)
DefaultAssay(x) <- ASSAY_USE

# Normalize if 'data' slot is missing
if (SLOT_USE == "data") {
  has_data_slot <- tryCatch({
    m <- GetAssayData(x, assay = ASSAY_USE, slot = "data")
    !is.null(m) && nrow(m) > 0 && ncol(m) > 0
  }, error = function(e) FALSE)
  
  if (!has_data_slot) {
    message("No 'data' slot found. Running NormalizeData (LogNormalize)...")
    x <- NormalizeData(x, normalization.method = "LogNormalize", 
                       scale.factor = 10000, verbose = FALSE)
  }
}

# Metadata handling
stopifnot(CELLTYPE_COL %in% colnames(x@meta.data))
x$celltype_simple <- collapse_celltype(x@meta.data[[CELLTYPE_COL]])

# Extract Expression Matrix
mat <- GetAssayData(x, assay = ASSAY_USE, slot = SLOT_USE)

missing_genes <- setdiff(GENES, rownames(mat))
if (length(missing_genes) > 0) {
  stop("Genes not found in assay: ", paste(missing_genes, collapse = ", "))
}

cells <- colnames(mat)

# Create Logical Matrix (True if expression > threshold)
expr_bool <- sapply(GENES, function(g) as.numeric(mat[g, cells]) > EXPR_THRESHOLD)
expr_bool <- as.matrix(expr_bool)
rownames(expr_bool) <- cells
colnames(expr_bool) <- GENES

# ==============================================================================
# 2. Evaluate All Subsets (Global)
# ==============================================================================
sub_list <- unlist(lapply(SUBSET_SIZES, function(k) all_subsets(GENES, k)), recursive = FALSE)

overall_sub <- do.call(rbind, lapply(sub_list, function(S) subset_binom_test(expr_bool, S)))
overall_sub$p_enrich_adj  <- p.adjust(overall_sub$p_enrich,  method = P_ADJUST_METHOD)
overall_sub$p_deplete_adj <- p.adjust(overall_sub$p_deplete, method = P_ADJUST_METHOD)

# Significance flags
overall_sub$call_enriched <- overall_sub$p_enrich_adj < 0.05
overall_sub$call_depleted <- overall_sub$p_deplete_adj < 0.05

# ==============================================================================
# 3. Evaluate All Subsets (Per Cell Type)
# ==============================================================================
grp <- x@meta.data[cells, "celltype_simple", drop = TRUE]
grp_levels <- sort(unique(grp))

bygrp_rows <- list()
for (g in grp_levels) {
  idx <- which(grp == g)
  if (length(idx) < MIN_CELLS_GROUP) next 
  
  m <- expr_bool[idx, , drop = FALSE]
  tmp <- do.call(rbind, lapply(sub_list, function(S) subset_binom_test(m, S)))
  tmp$group <- g
  bygrp_rows[[g]] <- tmp
}

bygrp_sub <- if (length(bygrp_rows) > 0) do.call(rbind, bygrp_rows) else data.frame()

if (nrow(bygrp_sub) > 0) {
  # Adjust p-values within each group
  bygrp_sub$p_enrich_adj <- ave(bygrp_sub$p_enrich, bygrp_sub$group,
                                FUN = function(p) p.adjust(p, method = P_ADJUST_METHOD))
  bygrp_sub$p_deplete_adj <- ave(bygrp_sub$p_deplete, bygrp_sub$group,
                                 FUN = function(p) p.adjust(p, method = P_ADJUST_METHOD))
  
  bygrp_sub$call_enriched <- bygrp_sub$p_enrich_adj < 0.05
  bygrp_sub$call_depleted <- bygrp_sub$p_deplete_adj < 0.05
}

# ==============================================================================
# 4. Pairwise Fisher Tests (Global + Per Cell Type)
# ==============================================================================
pairs <- combn(GENES, 2, simplify = FALSE)

# Global Pairs
overall_pair <- do.call(rbind, lapply(pairs, function(p) pair_fisher(expr_bool, p[1], p[2])))
overall_pair$p_greater_adj  <- p.adjust(overall_pair$p_greater, method = P_ADJUST_METHOD)
overall_pair$p_less_adj     <- p.adjust(overall_pair$p_less,    method = P_ADJUST_METHOD)
overall_pair$call_coexpress <- overall_pair$p_greater_adj < 0.05
overall_pair$call_exclusive <- overall_pair$p_less_adj    < 0.05

# Per Cell Type Pairs
bygrp_pair_rows <- list()
for (g in grp_levels) {
  idx <- which(grp == g)
  if (length(idx) < MIN_CELLS_GROUP) next
  
  m <- expr_bool[idx, , drop = FALSE]
  tmp <- do.call(rbind, lapply(pairs, function(p) pair_fisher(m, p[1], p[2])))
  tmp$group <- g
  
  # Adjust p-values within each group
  tmp$p_greater_adj <- p.adjust(tmp$p_greater, method = P_ADJUST_METHOD)
  tmp$p_less_adj    <- p.adjust(tmp$p_less,    method = P_ADJUST_METHOD)
  tmp$call_coexpress <- tmp$p_greater_adj < 0.05
  tmp$call_exclusive <- tmp$p_less_adj    < 0.05
  
  bygrp_pair_rows[[g]] <- tmp
}

bygrp_pair <- if (length(bygrp_pair_rows) > 0) do.call(rbind, bygrp_pair_rows) else data.frame()

# ==============================================================================
# 5. Export Results to Excel
# ==============================================================================
# Note: Filename may be long depending on the number of genes
out_xlsx <- file.path(OUTDIR, paste0("coexpression_significance_", paste(GENES, collapse = "_"), ".xlsx"))

wb <- createWorkbook()

# Sheet: README
addWorksheet(wb, "README")
writeData(wb, "README", data.frame(
  item = c("RDS_PATH", "assay_use", "slot_use", "expr_threshold", 
           "p_adjust_method", "subset_sizes", "grouping"),
  value = c(RDS_PATH, ASSAY_USE, SLOT_USE, EXPR_THRESHOLD, 
            P_ADJUST_METHOD, paste(SUBSET_SIZES, collapse = ","), 
            "celltype_simple = final_annotation minus _digits")
))

# Sheet: Global Subsets
addWorksheet(wb, "overall_subsets_binom")
writeData(wb, "overall_subsets_binom",
          overall_sub[order(overall_sub$p_enrich_adj, overall_sub$p_deplete_adj), ])

# Sheet: Group Subsets
addWorksheet(wb, "bygroup_subsets_binom")
writeData(wb, "bygroup_subsets_binom",
          if (nrow(bygrp_sub) > 0) bygrp_sub[order(bygrp_sub$group, bygrp_sub$p_enrich_adj), ] else bygrp_sub)

# Sheet: Global Pairs
addWorksheet(wb, "overall_pairs_fisher")
writeData(wb, "overall_pairs_fisher",
          overall_pair[order(pmin(overall_pair$p_greater_adj, overall_pair$p_less_adj)), ])

# Sheet: Group Pairs
addWorksheet(wb, "bygroup_pairs_fisher")
writeData(wb, "bygroup_pairs_fisher",
          if (nrow(bygrp_pair) > 0) bygrp_pair[order(bygrp_pair$group, pmin(bygrp_pair$p_greater_adj, bygrp_pair$p_less_adj)), ] else bygrp_pair)

# Sheet: Significant Summary
addWorksheet(wb, "SIGNIFICANT_SUMMARY")

sig_sub_overall  <- subset(overall_sub, call_enriched | call_depleted)
sig_pair_overall <- subset(overall_pair, call_coexpress | call_exclusive)

writeData(wb, "SIGNIFICANT_SUMMARY", "Overall subsets (binomial):", startRow = 1, startCol = 1)
writeData(wb, "SIGNIFICANT_SUMMARY", sig_sub_overall, startRow = 2, startCol = 1)

r0 <- 4 + nrow(sig_sub_overall)
writeData(wb, "SIGNIFICANT_SUMMARY", "Overall pairs (Fisher):", startRow = r0, startCol = 1)
writeData(wb, "SIGNIFICANT_SUMMARY", sig_pair_overall, startRow = r0 + 1, startCol = 1)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("Results saved to: ", normalizePath(out_xlsx))