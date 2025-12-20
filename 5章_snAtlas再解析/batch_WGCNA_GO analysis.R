#!/usr/bin/env Rscript

# ============================================================
# WGCNA Module GO Enrichment Analysis
#
# Description:
#   Identifies WGCNA modules highly correlated with a specific trait (e.g., BIL7)
#   and performs Gene Ontology (GO) enrichment analysis using clusterProfiler.
#
# Workflow:
#   1. Load WGCNA module-trait correlations and gene-module assignments.
#   2. Filter modules based on correlation threshold.
#   3. Perform GO enrichment (enrichGO) for target modules against the background universe.
#   4. Output results as CSV tables and dotplots.
# ============================================================

# ----------------------------
# Prerequisites (Uncomment to install)
# ----------------------------
# install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.At.tair.db", "dplyr", "ggplot2"))

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.At.tair.db)
  library(dplyr)
  library(ggplot2)
})

# =========================
# 0) Configuration
# =========================

# Directory containing WGCNA output subdirectories (e.g., "./WGCNA_Results")
BASE_OUTDIR <- "./WGCNA_Results" 

# Annotation Database for Arabidopsis thaliana
ORGDB       <- org.At.tair.db
KEYTYPE     <- "TAIR"           # Assumes AGI codes
ONTOLOGY    <- "BP"             # Choose from "BP", "MF", "CC"

# Correlation threshold for selecting modules of interest
# (Target modules with |cor_trait| >= threshold)
COR_THRESHOLD <- 0.6

# q-value cutoff for GO enrichment (after multiple testing correction)
QVALUE_CUTOFF <- 0.05

# =========================
# 1) Data Loading Functions
# =========================

read_wgcna_results <- function(organ) {
  organ_dir <- file.path(BASE_OUTDIR, organ)
  
  # Define file paths (Adjust naming convention if necessary)
  module_trait_file <- file.path(organ_dir, paste0("Module_trait_correlation_BIL7_", organ, ".csv"))
  all_genes_file    <- file.path(organ_dir, paste0("All_genes_module_and_kME_", organ, ".csv"))
  
  if (!file.exists(module_trait_file)) stop("File not found: ", module_trait_file)
  if (!file.exists(all_genes_file))    stop("File not found: ", all_genes_file)
  
  mt_df <- read.csv(module_trait_file, stringsAsFactors = FALSE)
  ag_df <- read.csv(all_genes_file,    stringsAsFactors = FALSE)
  
  list(
    module_trait = mt_df,
    all_genes    = ag_df,
    organ_dir    = organ_dir
  )
}

# =========================
# 2) GO Analysis Function (Single Module)
# =========================

run_go_for_module <- function(gene_vec,
                              universe_vec,
                              module_label,
                              organ,
                              outdir,
                              ont = ONTOLOGY,
                              q_cutoff = QVALUE_CUTOFF) {
  
  # Clean inputs: remove NAs and duplicates
  gene_vec     <- unique(na.omit(gene_vec))
  universe_vec <- unique(na.omit(universe_vec))
  
  # Skip if too few genes
  if (length(gene_vec) < 5) {
    warning("[", organ, "] Module ", module_label, " has fewer than 5 genes. Skipping GO.")
    return(NULL)
  }
  
  # Perform GO Enrichment
  ego <- tryCatch({
    enrichGO(
      gene          = gene_vec,
      universe      = universe_vec,
      OrgDb         = ORGDB,
      keyType       = KEYTYPE,
      ont           = ont,
      pAdjustMethod = "BH",
      qvalueCutoff  = q_cutoff,
      readable      = TRUE  # Maps gene IDs to symbols if available
    )
  }, error = function(e) {
    warning("enrichGO failed for ", module_label, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    warning("[", organ, "] Module ", module_label, " has no significant GO terms.")
    return(NULL)
  }
  
  # Save Results Table
  ego_df  <- as.data.frame(ego)
  out_csv <- file.path(outdir, paste0("GO_", ont, "_", organ, "_", module_label, ".csv"))
  write.csv(ego_df, out_csv, row.names = FALSE)
  message("[", organ, "] Saved GO table: ", out_csv)
  
  # Save Dotplot (Top 20 categories)
  dp <- dotplot(ego, showCategory = 20) +
    ggtitle(paste0("GO ", ont, " - ", organ, " - ", module_label)) +
    theme_bw() +
    theme(
      plot.title   = element_text(size = 10),
      axis.text.x  = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 10),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9),
      legend.title = element_text(size = 8),
      legend.text  = element_text(size = 7)
    )
  
  out_pdf <- file.path(outdir, paste0("GO_", ont, "_", organ, "_", module_label, "_dotplot.pdf"))
  ggsave(out_pdf, dp, width = 7, height = 6)
  message("[", organ, "] Saved GO dotplot: ", out_pdf)
  
  invisible(ego)
}

# =========================
# 3) Main Analysis Wrapper (Per Organ)
# =========================

run_go_for_organ <- function(organ) {
  message("========== [", organ, "] Starting GO Analysis ==========")
  
  # Load Data
  dat       <- read_wgcna_results(organ)
  mt_df     <- dat$module_trait
  ag_df     <- dat$all_genes
  organ_dir <- dat$organ_dir
  
  # Select Target Modules based on correlation threshold
  target_modules <- mt_df %>%
    dplyr::filter(abs(cor_BIL7) >= COR_THRESHOLD) %>%
    dplyr::arrange(desc(abs(cor_BIL7)))
  
  # Fallback: If no module passes the threshold, select the top correlated module
  if (nrow(target_modules) == 0) {
    warning("[", organ, "] No module passes |cor_BIL7| >= ", COR_THRESHOLD,
            ". Using only the top correlated module.")
    target_modules <- mt_df[order(-abs(mt_df$cor_BIL7)), , drop = FALSE][1, , drop = FALSE]
  }
  
  message("[", organ, "] Target modules identified:")
  print(target_modules[, c("module", "module_color", "cor_BIL7", "p_BIL7")])
  
  # Define Universe (Background gene set)
  universe_genes <- ag_df$gene
  
  # Iterate through target modules
  for (i in seq_len(nrow(target_modules))) {
    mod_label <- target_modules$module[i]       # e.g., "MEturquoise"
    mod_color <- target_modules$module_color[i] # e.g., "turquoise"
    
    message("[", organ, "] Processing module: ", mod_label, " (color=", mod_color, ")")
    
    # Extract genes belonging to the module
    genes_in_mod <- ag_df %>%
      dplyr::filter(module == mod_color & !is.na(gene)) %>%
      dplyr::pull(gene)
    
    if (length(genes_in_mod) == 0) {
      warning("[", organ, "] No genes found for module color ", mod_color)
      next
    }
    
    # Run GO Enrichment
    run_go_for_module(
      gene_vec     = genes_in_mod,
      universe_vec = universe_genes,
      module_label = mod_label,
      organ        = organ,
      outdir       = organ_dir,
      ont          = ONTOLOGY,
      q_cutoff     = QVALUE_CUTOFF
    )
  }
  
  message("========== [", organ, "] GO Analysis Completed ==========")
}

# =========================
# 4) Execution Example
# =========================

# Example usage:
# run_go_for_organ("S3_leaf_GLK1")

message("=== Script loaded successfully. Ready to run analysis. ===")