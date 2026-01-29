#!/usr/bin/env Rscript

# ============================================================
# WGCNA on Pseudobulk Bins (Merged Datasets)
#
# Description:
#   Performs WGCNA on a merged dataset composed of multiple scRNA-seq samples (e.g., S1, S2, S3).
#   1. Loads and merges multiple Seurat objects.
#   2. Filters cells expressing a specific Index Gene.
#   3. Sorts and bins cells based on Index Gene expression (Quantile binning).
#   4. Constructs pseudobulk expression profiles for each bin.
#   5. Performs WGCNA to identify co-expression modules correlated with the Index Gene.
#
# Input:
#   - Multiple Seurat RDS files.
#
# Output:
#   - Module assignments, Eigengenes, Correlational stats, and Plots.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(WGCNA)
  library(Matrix)
  library(matrixStats)
})

options(stringsAsFactors = FALSE)
# Note: Threading support depends on OS/R config
enableWGCNAThreads(nThreads = 10)

# ==============================================================================
# Configuration
# ==============================================================================
# Input Files
SEURAT_RDS_1 <- "path/to/S1_leaf_data.rds"
SEURAT_RDS_2 <- "path/to/S2_leaf_data.rds"
SEURAT_RDS_3 <- "path/to/S3_leaf_data.rds"

# Output Directory
OUT_DIR      <- "path/to/output/WGCNA_youngLeaf_binned"

# Analysis Target
INDEX_GENE   <- "AT1G63720" # Gene used for sorting/binning
N_BINS       <- 50          # Number of pseudobulk bins

# Data Extraction
ASSAY_USE    <- NULL        # NULL = DefaultAssay
SLOT_USE     <- "data"      # Normalized data

# WGCNA Parameters
NET_TYPE         <- "signed"
COR_TYPE         <- "pearson"
MAX_P_OUTLIERS   <- 1
MAD_Q            <- 0.25
MIN_MAD          <- 0.01
MIN_MODULE_SIZE  <- 50
MERGE_CUT_HEIGHT <- 0.3

# ==============================================================================
# 1. Setup and Load Data
# ==============================================================================
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUT_DIR)

# Verify inputs
input_files <- c(SEURAT_RDS_1, SEURAT_RDS_2, SEURAT_RDS_3)
if (!all(file.exists(input_files))) {
  stop("One or more input RDS files not found.")
}

# Load Seurat Objects
obj1 <- readRDS(SEURAT_RDS_1)
obj2 <- readRDS(SEURAT_RDS_2)
obj3 <- readRDS(SEURAT_RDS_3)

# Rename cells to ensure uniqueness during merge
obj1 <- RenameCells(obj1, add.cell.id = "S1_leaf")
obj2 <- RenameCells(obj2, add.cell.id = "S2_leaf")
obj3 <- RenameCells(obj3, add.cell.id = "S3_leaf")

# Add sample metadata
obj1$sample_id <- "S1_leaf"
obj2$sample_id <- "S2_leaf"
obj3$sample_id <- "S3_leaf"

# Merge Objects
message("Merging datasets...")
obj <- merge(obj1, y = list(obj2, obj3))

# Cleanup individual objects
rm(obj1, obj2, obj3); gc()

# Prepare Expression Matrix
if (is.null(ASSAY_USE)) ASSAY_USE <- DefaultAssay(obj)

expr <- GetAssayData(obj, assay = ASSAY_USE, slot = SLOT_USE)
# Remove genes with 0 total counts
expr <- expr[Matrix::rowSums(expr) > 0, , drop = FALSE] * 1.0

if (!(INDEX_GENE %in% rownames(expr))) {
  stop("Index gene not found in merged dataset: ", INDEX_GENE)
}

# ==============================================================================
# 2. Filter Cells by Index Gene Expression
# ==============================================================================
idx_expr <- expr[INDEX_GENE, ]

# Keep cells with non-zero expression of index gene
cells_use <- which(idx_expr > 0)
message("Cells expressing ", INDEX_GENE, ": ", length(cells_use))

if (length(cells_use) < N_BINS) {
  stop("Not enough cells expressing index gene to form ", N_BINS, " bins.")
}

expr <- expr[, cells_use, drop = FALSE]
idx_expr <- idx_expr[cells_use]

# ==============================================================================
# 3. Binning Based on Expression Quantiles
# ==============================================================================
bin_id <- cut(
  idx_expr,
  breaks = quantile(idx_expr, probs = seq(0, 1, length.out = N_BINS + 1)),
  include.lowest = TRUE,
  labels = paste0("bin", seq_len(N_BINS))
)

message("Bin distribution:")
print(table(bin_id))

# ==============================================================================
# 4. Construct Pseudobulk Matrix (Genes x Bins)
# ==============================================================================
cell_idx_list <- split(seq_len(ncol(expr)), bin_id)

# Calculate mean expression per bin
pb <- sapply(names(cell_idx_list), function(b) {
  Matrix::rowMeans(expr[, cell_idx_list[[b]], drop = FALSE])
})
pb <- as.matrix(pb) # Genes x Bins
colnames(pb) <- names(cell_idx_list)

write.csv(pb, "pseudobulk_by_bins.csv")

# ==============================================================================
# 5. Gene Filtering by MAD
# ==============================================================================
data0 <- t(pb) # Bins x Genes for WGCNA

m.mad <- matrixStats::colMads(data0, constant = 1)
thr <- max(as.numeric(quantile(m.mad, probs = MAD_Q)), MIN_MAD)
keep <- which(m.mad > thr)

# Ensure Index Gene is retained
keep <- union(keep, which(colnames(data0) == INDEX_GENE))

dataExpr <- as.data.frame(data0[, keep, drop = FALSE])

# Check for outliers
gsg <- goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK) {
  dataExpr <- dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nSamples <- nrow(dataExpr)
nGenes   <- ncol(dataExpr)
message("WGCNA Input: ", nSamples, " bins, ", nGenes, " genes")

# Define Trait (Index Gene Expression Profile)
trait <- data.frame(
  IndexGeneExpr = dataExpr[, INDEX_GENE],
  row.names = rownames(dataExpr)
)

# ==============================================================================
# 6. WGCNA Execution
# ==============================================================================
# Sample Clustering
sampleTree <- hclust(dist(dataExpr), method = "average")
pdf("sample_clustering.pdf")
plot(sampleTree, main = "Sample Clustering", sub = "", xlab = "")
dev.off()

# Pick Soft Threshold
disableWGCNAThreads()
powers <- c(1:10, seq(12, 30, 2))
sft <- pickSoftThreshold(dataExpr, powerVector = powers, networkType = NET_TYPE, verbose = 5)
enableWGCNAThreads(nThreads = 10)

power <- sft$powerEstimate
if (is.na(power)) {
  power <- 6
  warning("Soft threshold power estimation failed. Using default: ", power)
}
message("Chosen power: ", power)

# Network Construction & Module Detection
net <- blockwiseModules(
  dataExpr,
  power = power,
  maxBlockSize = nGenes,
  TOMType = NET_TYPE,
  minModuleSize = MIN_MODULE_SIZE,
  reassignThreshold = 0,
  mergeCutHeight = MERGE_CUT_HEIGHT,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  corType = COR_TYPE,
  maxPOutliers = MAX_P_OUTLIERS,
  verbose = 3
)

if (length(unique(net$colors)) == 1 && unique(net$colors) == 0) {
  writeLines("No modules detected.", "NO_MODULES.txt")
  message("No modules detected. Exiting.")
  quit(save = "no")
}

moduleColors <- labels2colors(net$colors)

# Dendrogram Plot
pdf("module_colors.pdf")
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  addGuide = TRUE,
  hang = 0.03
)
dev.off()

# ==============================================================================
# 7. Export Results
# ==============================================================================
# Module Eigengenes
MEs <- orderMEs(net$MEs)
write.csv(MEs, "MEs.csv")

# Correlation with Trait
modTraitCor <- cor(MEs, trait$IndexGeneExpr, use = "p")
modTraitP   <- corPvalueStudent(as.matrix(modTraitCor), nSamples)

write.csv(modTraitCor, "ME_trait_cor.csv")
write.csv(modTraitP,   "ME_trait_pvalue.csv")

# Gene-Module Map
gene_module <- data.frame(
  gene = colnames(dataExpr),
  module = labels2colors(net$colors),
  module_label = net$colors,
  stringsAsFactors = FALSE
)

write.csv(gene_module, "gene_module_map.csv", row.names = FALSE)

message("Analysis Complete. Results saved to: ", normalizePath(OUT_DIR))