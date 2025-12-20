#!/usr/bin/env Rscript

# ============================================================
# scRNA-seq Weighted Gene Co-expression Network Analysis (WGCNA)
#
# Description:
#   1. Loads a Seurat RDS file.
#   2. Aggregates expression data into pseudo-bulk (by celltype/stage).
#   3. Performs WGCNA to identify gene modules.
#   4. Correlates modules with a specific trait gene (e.g., BIL7).
#   5. Identifies hub genes within the most correlated module.
#
# Input:  Seurat RDS file
# Output: WGCNA results (Module-trait correlations, Hub genes, kME tables)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(reshape2)
  library(WGCNA)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# =========================
# 0) Configuration
# =========================

# Input File Path
RDS_PATH   <- "./input_data/scRNA_dataset.RDS" 

# Output Identification Name
ORGAN_NAME <- "S3_leaf_GLK1"  # e.g., "whole_plant", "leaf", "merged"

# Metadata Columns
CELLTYPE_COL <- "integrated_annotation"  # Column for cell type annotation
STAGE_COL    <- NULL                     # Column for stage (set NULL if not used)

# Target Trait Gene (AGI Code)
TRAIT_GENE <- "AT2G20570"  # Target gene to correlate modules with

# Output Directory
BASE_OUTDIR <- "./WGCNA_Results"
dir.create(BASE_OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Filtering Threshold
MIN_MEAN_EXPR <- 0.05   # Minimum mean expression threshold for filtering genes

# =========================
# A) Gene Metrics Function
# =========================
# Aggregates single-cell data into pseudo-bulk average expression matrices

build_gene_metrics <- function(x, celltype_col, stage_col = NULL) {
  stopifnot(inherits(x, "Seurat"))
  DefaultAssay(x) <- "RNA"
  
  # Ensure data slot is populated
  if (!"data" %in% slotNames(x[["RNA"]])) {
    message("No 'data' slot found. Running NormalizeData (LogNormalize, scale.factor=1e4).")
    x <- NormalizeData(x, normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = FALSE)
  }
  
  # Set Identities and grouping column
  stopifnot(celltype_col %in% colnames(x@meta.data))
  Idents(x) <- factor(x@meta.data[[celltype_col]])
  
  if (!is.null(stage_col)) {
    stopifnot(stage_col %in% colnames(x@meta.data))
    x$.__group__ <- interaction(x@meta.data[[stage_col]],
                                as.character(Idents(x)),
                                drop = TRUE, sep = "|")
    grouping_desc <- paste0(stage_col, " x ", celltype_col)
  } else {
    x$.__group__ <- factor(as.character(Idents(x)))
    grouping_desc <- celltype_col
  }
  groups <- levels(x$.__group__)
  
  # Calculate Average Expression (slot='data')
  avg_obj <- AverageExpression(
    x, group.by = ".__group__", assays = DefaultAssay(x),
    slot = "data", verbose = FALSE
  )
  
  # Handle output format (Seurat v3 returns list, v5 might return matrix directly)
  if (is.list(avg_obj)) {
    if (DefaultAssay(x) %in% names(avg_obj)) {
      avg_mat <- avg_obj[[DefaultAssay(x)]]
    } else {
      avg_mat <- avg_obj[[1]]
    }
  } else {
    avg_mat <- avg_obj
  }
  
  # Ensure matrix format
  if (inherits(avg_mat, "data.frame") || inherits(avg_mat, "tbl_df")) avg_mat <- as.matrix(avg_mat)
  if (inherits(avg_mat, "Matrix") || inherits(avg_mat, "dgCMatrix")) avg_mat <- as.matrix(avg_mat)
  if (!is.matrix(avg_mat)) stop("AverageExpression could not be converted to matrix: ", paste(class(avg_mat), collapse=", "))
  
  if (is.null(rownames(avg_mat))) rownames(avg_mat) <- rownames(x)
  if (is.null(colnames(avg_mat))) colnames(avg_mat) <- groups
  
  # Calculate Detection Rate (Percent Expressed)
  # Uses 'counts > 0' if available, otherwise 'data > 0'
  cnt <- GetAssayData(x, assay = DefaultAssay(x), slot = "counts")
  
  if (Matrix::nnzero(cnt) == 0) {
    message("Counts slot is empty. Using data > 0 as detection proxy.")
    dat <- GetAssayData(x, assay = DefaultAssay(x), slot = "data")
    pct_mat <- sapply(groups, function(g) {
      cells_g <- colnames(x)[x$.__group__ == g]
      Matrix::rowMeans(dat[, cells_g, drop = FALSE] > 0)
    })
  } else {
    pct_mat <- sapply(groups, function(g) {
      cells_g <- colnames(x)[x$.__group__ == g]
      Matrix::rowMeans(cnt[, cells_g, drop = FALSE] > 0)
    })
  }
  
  if (is.null(dim(pct_mat))) {
    pct_mat <- matrix(pct_mat, ncol = 1, dimnames = list(rownames(x), groups))
  } else {
    rownames(pct_mat) <- rownames(x)
    colnames(pct_mat) <- groups
  }
  
  list(
    x = x,
    grouping_desc = grouping_desc,
    groups = groups,
    avg_mat = avg_mat,   # Rows=Genes, Cols=Groups
    pct_mat = pct_mat    # Rows=Genes, Cols=Groups (0-1)
  )
}

# =========================
# B) WGCNA Analysis Function
# =========================

run_wgcna_for_organ <- function(
    gm,
    organ_name,
    trait_gene,
    base_outdir,
    min_mean_expr = 0.05
) {
  message("========== Processing: ", organ_name, " ==========")
  outdir <- file.path(base_outdir, organ_name)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  expr_mat <- gm$avg_mat  # genes x groups
  message("[", organ_name, "] Expression Matrix Dimensions: ", paste(dim(expr_mat), collapse = " x "))
  
  # ---- Prepare Trait Data (Expression of Target Gene) ----
  if (!trait_gene %in% rownames(expr_mat)) {
    stop("TRAIT_GENE (", trait_gene, ") not found in expression matrix for ", organ_name)
  }
  trait_vec <- as.numeric(expr_mat[trait_gene, ])
  names(trait_vec) <- colnames(expr_mat)
  
  # ---- Prepare datExpr (Samples x Genes) ----
  datExpr0 <- t(as.matrix(expr_mat))
  datTraits <- data.frame(BIL7 = trait_vec)
  rownames(datTraits) <- rownames(datExpr0)
  
  # ---- Gene Filtering ----
  message("[", organ_name, "] Filtering genes...")
  gene_means <- apply(datExpr0, 2, mean)
  keep_genes <- gene_means > min_mean_expr
  datExpr <- datExpr0[, keep_genes]
  message("[", organ_name, "] Genes retained for WGCNA: ", ncol(datExpr))
  
  # ---- Sample/Gene Quality Check ----
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0)
      message("[", organ_name, "] Removing genes: ", paste(colnames(datExpr)[!gsg$goodGenes], collapse = ", "))
    if (sum(!gsg$goodSamples) > 0)
      message("[", organ_name, "] Removing samples: ", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", "))
    
    datExpr   <- datExpr[gsg$goodSamples, gsg$goodGenes]
    datTraits <- datTraits[gsg$goodSamples, , drop = FALSE]
  }
  
  # ---- Soft Power Selection ----
  message("[", organ_name, "] Picking soft threshold...")
  powers <- c(1:20)
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  
  R2 <- sft$fitIndices[, "SFT.R.sq"]
  softPower <- NA
  for (p in powers) {
    idx <- which(sft$fitIndices[, "Power"] == p)
    if (length(idx) == 1 && R2[idx] > 0.85) {
      softPower <- p
      break
    }
  }
  if (is.na(softPower)) {
    warning("[", organ_name, "] No power found with R^2 > 0.85. Defaulting to power = 6.")
    softPower <- 6
  }
  message("[", organ_name, "] Selected softPower: ", softPower)
  
  # ---- Module Detection (blockwiseModules) ----
  message("[", organ_name, "] Detecting modules...")
  net <- blockwiseModules(
    datExpr,
    power = softPower,
    TOMType = "signed",
    minModuleSize = 50,
    reassignThreshold = 0,
    mergeCutHeight = 0.30,
    numericLabels = FALSE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    verbose = 3
  )
  
  moduleColors <- labels2colors(net$colors)
  MEs0 <- net$MEs
  MEs  <- orderMEs(MEs0)
  
  # ---- Module-Trait Correlation ----
  message("[", organ_name, "] Calculating Module-Trait correlations...")
  MEtraitCor    <- cor(MEs, datTraits$BIL7, use = "p")
  MEtraitPvalue <- corPvalueStudent(MEtraitCor, nrow(datExpr))
  
  mt_df <- data.frame(
    module     = rownames(MEtraitCor),
    cor_BIL7   = as.numeric(MEtraitCor[, 1]),
    p_BIL7     = as.numeric(MEtraitPvalue[, 1])
  )
  mt_df$module_color <- gsub("^ME", "", mt_df$module)
  mt_df <- mt_df[order(-abs(mt_df$cor_BIL7)), ]
  
  write.csv(mt_df,
            file = file.path(outdir, paste0("Module_trait_correlation_BIL7_", organ_name, ".csv")),
            row.names = FALSE)
  message("[", organ_name, "] Saved: Module trait correlation table")
  
  # ---- Analyze the Top Correlated Module ----
  best_module <- mt_df$module[1]       # e.g., "MEturquoise"
  best_color  <- mt_df$module_color[1] # e.g., "turquoise"
  message("[", organ_name, "] Most correlated module: ", best_module, " (color=", best_color, ")")
  
  modGenes <- moduleColors == best_color
  genes_in_mod <- colnames(datExpr)[modGenes]
  
  # Calculate Eigengene-based Connectivity (kME)
  kME_all <- as.data.frame(signedKME(datExpr, MEs, outputColumnName = "kME"))
  kME_colname <- paste0("kME", gsub("^ME", "", best_module))
  
  if (!kME_colname %in% colnames(kME_all)) {
    stop("[", organ_name, "] kME column not found: ", kME_colname)
  }
  
  kME_mod <- kME_all[modGenes, kME_colname, drop = FALSE]
  
  hub_df <- data.frame(
    gene = rownames(kME_mod),
    kME  = kME_mod[, 1]
  )
  hub_df <- hub_df[order(-hub_df$kME), ]
  
  # Save Top Hub Genes
  n_top <- min(50, nrow(hub_df))
  hub_top <- hub_df[seq_len(n_top), ]
  
  write.csv(
    hub_top,
    file = file.path(outdir, paste0("Hub_genes_BIL7_", organ_name, "_", best_module, "_top", n_top, ".csv")),
    row.names = FALSE
  )
  message("[", organ_name, "] Saved: Top ", n_top, " hub genes")
  
  # Save All Genes with Module Assignment and kME
  all_genes_df <- data.frame(
    gene   = colnames(datExpr),
    module = moduleColors
  )
  all_genes_df <- merge(all_genes_df, hub_df, by = "gene", all.x = TRUE)
  write.csv(all_genes_df,
            file = file.path(outdir, paste0("All_genes_module_and_kME_", organ_name, ".csv")),
            row.names = FALSE)
  
  message("[", organ_name, "] Saved: All genes module assignment table")
  message("========== [", organ_name, "] Analysis Complete ==========")
  
  invisible(list(
    datExpr      = datExpr,
    datTraits    = datTraits,
    net          = net,
    moduleColors = moduleColors,
    MEs          = MEs,
    moduleTrait  = mt_df,
    hub_top      = hub_top
  ))
}

# =========================
# C) Execution
# =========================

message("[1] Loading RDS file: ", RDS_PATH)
obj <- readRDS(RDS_PATH)

message("[2] Building gene metrics for: ", ORGAN_NAME)
gm <- build_gene_metrics(obj, celltype_col = CELLTYPE_COL, stage_col = STAGE_COL)

message("[3] Running WGCNA...")
res <- run_wgcna_for_organ(
  gm            = gm,
  organ_name    = ORGAN_NAME,
  trait_gene    = TRAIT_GENE,
  base_outdir   = BASE_OUTDIR,
  min_mean_expr = MIN_MEAN_EXPR
)

message("=== All processes finished successfully. ===")