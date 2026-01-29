#!/usr/bin/env Rscript

# ============================================================
# WGCNA on Pseudobulk Bins (Single-Gene Gradient)
#
# Description:
#   Performs Weighted Gene Co-expression Network Analysis (WGCNA) on single-cell data
#   sorted along the expression gradient of a specific "Index Gene".
#   1. Filters cells expressing the Index Gene.
#   2. Sorts and bins cells based on Index Gene expression (Quantile binning).
#   3. Constructs pseudobulk expression profiles for each bin.
#   4. Filters genes based on Median Absolute Deviation (MAD).
#   5. Constructs WGCNA network and identifies modules.
#   6. Correlates modules with the Index Gene expression pattern.
#
# Input:
#   - Seurat RDS file.
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
# Note: Threading support depends on the OS and R installation
enableWGCNAThreads(nThreads = 10)

# ==============================================================================
# Configuration
# ==============================================================================
# Input/Output
SEURAT_RDS      <- "path/to/your/scRNA_data.rds"
OUT_DIR         <- "path/to/output/WGCNA_Analysis"
INDEX_GENE      <- "AT5G52430"  # Gene used for sorting/binning cells

# Data Extraction Settings
ASSAY_USE       <- NULL         # NULL = Use DefaultAssay
SLOT_USE        <- "data"       # Normalized data (e.g., LogNormalize)
N_BINS          <- 52           # Number of pseudobulk bins

# WGCNA Parameters
NET_TYPE        <- "signed"
COR_TYPE        <- "pearson"
MAX_P_OUTLIERS  <- 1
MAD_Q           <- 0.25         # Quantile for MAD filtering
MIN_MAD         <- 0.01         # Minimum MAD threshold
MIN_MODULE_SIZE <- 50
MERGE_CUT_HEIGHT<- 0.3

# ==============================================================================
# 1. Setup and Load Data
# ==============================================================================
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUT_DIR)

if (!file.exists(SEURAT_RDS)) stop("Seurat file not found: ", SEURAT_RDS)

obj <- readRDS(SEURAT_RDS)

if (is.null(ASSAY_USE)) ASSAY_USE <- DefaultAssay(obj)

# Extract expression matrix
expr <- GetAssayData(obj, assay = ASSAY_USE, slot = SLOT_USE)

# Remove genes with 0 counts across all cells (Global filter)
expr <- expr[Matrix::rowSums(expr) > 0, , drop = FALSE] * 1.0

if (!(INDEX_GENE %in% rownames(expr))) {
  stop("Index gene not found in the dataset: ", INDEX_GENE)
}

# ==============================================================================
# 2. Filter Cells by Index Gene Expression
# ==============================================================================
# Extract expression of the index gene
idx_expr <- expr[INDEX_GENE, ]

# Keep only cells where the index gene is expressed (> 0)
cells_use <- which(idx_expr > 0)
message("Cells expressing ", INDEX_GENE, ": ", length(cells_use))

if (length(cells_use) < N_BINS) {
  stop("Not enough cells expressing the index gene to form ", N_BINS, " bins.")
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
# 5. Gene Filtering by MAD (Median Absolute Deviation)
# ==============================================================================
data0 <- t(pb) # Transpose to Bins x Genes for WGCNA

m.mad <- matrixStats::colMads(data0, constant = 1)
thr <- max(as.numeric(quantile(m.mad, probs = MAD_Q)), MIN_MAD)
keep <- which(m.mad > thr)

# Ensure Index Gene is kept for correlation analysis
keep <- union(keep, which(colnames(data0) == INDEX_GENE))

dataExpr <- as.data.frame(data0[, keep, drop = FALSE])

# Check for samples/genes with too many missing values
gsg <- goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK) {
  dataExpr <- dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nSamples <- nrow(dataExpr)
nGenes   <- ncol(dataExpr)
message("WGCNA Input: ", nSamples, " bins, ", nGenes, " genes")

# Define Trait: Mean expression of the index gene in each bin
trait <- data.frame(
  IndexGeneExpr = dataExpr[, INDEX_GENE],
  row.names = rownames(dataExpr)
)

# ==============================================================================
# 6. WGCNA Execution
# ==============================================================================
# Sample Clustering Check
sampleTree <- hclust(dist(dataExpr), method = "average")
pdf("sample_clustering.pdf")
plot(sampleTree, main = "Sample Clustering", sub = "", xlab = "")
dev.off()

# Pick Soft Threshold
disableWGCNAThreads() # Often required for pickSoftThreshold stability
powers <- c(1:10, seq(12, 30, 2))
sft <- pickSoftThreshold(dataExpr, powerVector = powers, networkType = NET_TYPE, verbose = 5)
enableWGCNAThreads(nThreads = 10)

power <- sft$powerEstimate
if (is.na(power)) {
  power <- 6
  warning("Soft threshold power could not be estimated. Using default: ", power)
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

# Plot Dendrogram
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
# Module Eigengenes (MEs)
MEs <- orderMEs(net$MEs)
write.csv(MEs, "MEs.csv")

# Module-Trait Correlation
modTraitCor <- cor(MEs, trait$IndexGeneExpr, use = "p")
modTraitP   <- corPvalueStudent(as.matrix(modTraitCor), nSamples)

write.csv(modTraitCor, "ME_trait_cor.csv")
write.csv(modTraitP,   "ME_trait_pvalue.csv")

# Gene-Module Mapping
gene_module <- data.frame(
  gene = colnames(dataExpr),
  module = labels2colors(net$colors),
  module_label = net$colors,
  stringsAsFactors = FALSE
)

write.csv(gene_module, "gene_module_map.csv", row.names = FALSE)

message("Analysis Complete. Results saved to: ", normalizePath(OUT_DIR))