

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
# 0) ★設定★
# =========================

# ★★★ 変更箇所: 1つのRDSファイルを指定してください ★★★
RDS_PATH   <- "/Users/ryo/Desktop/scAtlas/S3_leaf_rename.slim.RDS" 

# ★★★ 変更箇所: 出力ファイル名やフォルダ名に使われる名前を指定してください ★★★
ORGAN_NAME <- "S3_leaf_GLK1"  # 例: "whole_plant", "leaf", "merged" など

celltype_col <- "integrated_annotation"  # メタデータの細胞型列
stage_col    <- NULL                     # 必要なければ NULL

# BIL7 の AGI（ここを本物に変えてください）
TRAIT_GENE <- "AT2G20570"  # 例: "AT5Gxxxxx" など

BASE_OUTDIR <- "/Users/ryo/Desktop/scAtlas"
dir.create(BASE_OUTDIR, showWarnings = FALSE, recursive = TRUE)

MIN_MEAN_EXPR <- 0.05   # 低発現遺伝子フィルタのしきい値

# =========================
# A) Gene metrics 関数（変更なし）
# =========================

build_gene_metrics <- function(x, celltype_col, stage_col = NULL) {
  stopifnot(inherits(x, "Seurat"))
  DefaultAssay(x) <- "RNA"
  if (!"data" %in% slotNames(x[["RNA"]])) {
    message("No 'data' slot found. NormalizeData(LogNormalize, scale.factor=1e4).")
    x <- NormalizeData(x, normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = FALSE)
  }
  
  # Idents とグループ列
  stopifnot(celltype_col %in% colnames(x@meta.data))
  Idents(x) <- factor(x@meta.data[[celltype_col]])
  if (!is.null(stage_col)) {
    stopifnot(stage_col %in% colnames(x@meta.data))
    x$.__group__ <- interaction(x@meta.data[[stage_col]],
                                as.character(Idents(x)),
                                drop = TRUE, sep = "|")
    grouping_desc <- paste0(stage_col, " × ", celltype_col)
  } else {
    x$.__group__ <- factor(as.character(Idents(x)))
    grouping_desc <- celltype_col
  }
  groups <- levels(x$.__group__)
  
  # 平均発現（slot='data'）
  avg_obj <- AverageExpression(
    x, group.by = ".__group__", assays = DefaultAssay(x),
    slot = "data", verbose = FALSE
  )
  if (is.list(avg_obj)) {
    if (DefaultAssay(x) %in% names(avg_obj)) {
      avg_mat <- avg_obj[[DefaultAssay(x)]]
    } else {
      avg_mat <- avg_obj[[1]]
    }
  } else {
    avg_mat <- avg_obj
  }
  if (inherits(avg_mat, "data.frame") || inherits(avg_mat, "tbl_df")) avg_mat <- as.matrix(avg_mat)
  if (inherits(avg_mat, "Matrix") || inherits(avg_mat, "dgCMatrix")) avg_mat <- as.matrix(avg_mat)
  if (!is.matrix(avg_mat)) stop("AverageExpression を matrix にできません: ", paste(class(avg_mat), collapse=", "))
  
  if (is.null(rownames(avg_mat))) rownames(avg_mat) <- rownames(x)
  if (is.null(colnames(avg_mat))) colnames(avg_mat) <- groups
  
  # 検出率（counts>0）。counts が空なら data>0 を代用
  cnt <- GetAssayData(x, assay = DefaultAssay(x), slot = "counts")
  if (Matrix::nnzero(cnt) == 0) {
    message("counts slot is empty. Using data>0 as detection proxy.")
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
    avg_mat = avg_mat,   # 行=遺伝子, 列=group
    pct_mat = pct_mat    # 行=遺伝子, 列=group（0~1）
  )
}

# =========================
# B) WGCNA を実行する関数（変更なし）
# =========================

run_wgcna_for_organ <- function(
    gm,
    organ_name,
    trait_gene,
    base_outdir,
    min_mean_expr = 0.05
) {
  message("========== [", organ_name, "] ==========")
  outdir <- file.path(base_outdir, organ_name)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  expr_mat <- gm$avg_mat  # genes x groups（ここでは group=celltype）
  message("[", organ_name, "] expr_mat: genes x groups = ", 
          paste(dim(expr_mat), collapse = " x "))
  
  # ---- trait (BIL7 平均発現) ----
  if (!trait_gene %in% rownames(expr_mat)) {
    stop("TRAIT_GENE (", trait_gene, ") が ", organ_name,
         " の expr_mat に見つかりません。AGI 名を確認してください。")
  }
  trait_vec <- as.numeric(expr_mat[trait_gene, ])
  names(trait_vec) <- colnames(expr_mat)  # group 名（celltype）
  
  # ---- datExpr: samples x genes ----
  datExpr0 <- t(as.matrix(expr_mat))  # groups x genes
  datTraits <- data.frame(BIL7 = trait_vec)
  rownames(datTraits) <- rownames(datExpr0)
  
  # ---- フィルタリング ----
  message("[", organ_name, "] Gene filtering...")
  gene_means <- apply(datExpr0, 2, mean)
  keep_genes <- gene_means > min_mean_expr
  datExpr <- datExpr0[, keep_genes]
  message("[", organ_name, "] Genes kept for WGCNA: ", ncol(datExpr))
  
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0)
      message("[", organ_name, "] Removing genes: ",
              paste(colnames(datExpr)[!gsg$goodGenes], collapse = ", "))
    if (sum(!gsg$goodSamples) > 0)
      message("[", organ_name, "] Removing samples: ",
              paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", "))
    datExpr  <- datExpr[gsg$goodSamples, gsg$goodGenes]
    datTraits <- datTraits[gsg$goodSamples, , drop = FALSE]
  }
  
  # ---- power 選択 ----
  message("[", organ_name, "] pickSoftThreshold...")
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
    warning("[", organ_name, "] R^2 > 0.85 を満たす power が見つからなかったため、power=6 を仮に使用します。")
    softPower <- 6
  }
  message("[", organ_name, "] Selected softPower: ", softPower)
  
  # ---- モジュール検出 ----
  message("[", organ_name, "] blockwiseModules...")
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
  
  # ---- モジュールと BIL7 trait の相関 ----
  message("[", organ_name, "] ME-trait correlation...")
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
  message("[", organ_name, "] Saved: Module_trait_correlation_BIL7_", organ_name, ".csv")
  
  # ---- BIL7 と最も相関が強いモジュールの hub 抽出 ----
  best_module <- mt_df$module[1]        # 例: "MEturquoise"
  best_color  <- mt_df$module_color[1]  # 例: "turquoise"
  message("[", organ_name, "] Most BIL7-correlated module: ",
          best_module, " (color=", best_color, ")")
  
  modGenes <- moduleColors == best_color
  genes_in_mod <- colnames(datExpr)[modGenes]
  
  # eigengene membership (kME)
  kME_all <- as.data.frame(signedKME(datExpr, MEs, outputColumnName = "kME"))
  kME_colname <- paste0("kME", gsub("^ME", "", best_module))  # "MEturquoise" -> "kMEturquoise"
  
  if (!kME_colname %in% colnames(kME_all)) {
    stop("[", organ_name, "] kME column not found: ", kME_colname,
         " ; available: ", paste(colnames(kME_all), collapse = ", "))
  }
  
  kME_mod <- kME_all[modGenes, kME_colname, drop = FALSE]
  
  hub_df <- data.frame(
    gene = rownames(kME_mod),
    kME  = kME_mod[, 1]
  )
  hub_df <- hub_df[order(-hub_df$kME), ]
  
  n_top <- min(50, nrow(hub_df))
  hub_top <- hub_df[seq_len(n_top), ]
  
  write.csv(
    hub_top,
    file = file.path(outdir,
                     paste0("Hub_genes_BIL7_", organ_name, "_", best_module, "_top", n_top, ".csv")),
    row.names = FALSE
  )
  message("[", organ_name, "] Saved: Hub_genes_BIL7_", organ_name, "_", best_module,
          "_top", n_top, ".csv")
  
  # 全遺伝子＋モジュール色＋kME も出しておく
  all_genes_df <- data.frame(
    gene   = colnames(datExpr),
    module = moduleColors
  )
  all_genes_df <- merge(all_genes_df, hub_df, by = "gene", all.x = TRUE)
  write.csv(all_genes_df,
            file = file.path(outdir, paste0("All_genes_module_and_kME_", organ_name, ".csv")),
            row.names = FALSE)
  
  message("[", organ_name, "] Saved: All_genes_module_and_kME_", organ_name, ".csv")
  message("========== [", organ_name, "] DONE ==========")
  
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
# C) 実行処理（1つのRDSのみ処理）
# =========================

message("[1] Read RDS: ", RDS_PATH)
obj <- readRDS(RDS_PATH)

message("[2] Build gene metrics for: ", ORGAN_NAME)
gm <- build_gene_metrics(obj, celltype_col = celltype_col, stage_col = stage_col)

# WGCNA 実行
res <- run_wgcna_for_organ(
  gm          = gm,
  organ_name  = ORGAN_NAME,
  trait_gene  = TRAIT_GENE,
  base_outdir = BASE_OUTDIR,
  min_mean_expr = MIN_MEAN_EXPR
)

message("=== All done. ===")

