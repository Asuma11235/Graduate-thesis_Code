
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(reshape2)  # melt 用
})

# =========================
# 0) ★設定（必要なら編集）★
# =========================
RDS_PATH     <- "/Users/ryo/Desktop/scAtlas/stem_rename.slim.RDS"
celltype_col <- "integrated_annotation"  # 例: "final_annotation" / "integrated_annotation"
stage_col    <- NULL                      # 例: "stage"（無ければ NULL）
# =========================


# -------------------------
# A) Gene metrics 構築関数
# -------------------------
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

# -----------------------------------------
# 絶対値（非相対化）の棒グラフを作成・保存する関数（新規）
# ・横軸=Group、縦軸=平均発現（slot='data'のAverageExpression）
# ・色=AGI（遺伝子）、並列バー
# ・必要に応じてログ軸も指定可能（y_transform="log10", pseudo_countで0回避）
# -----------------------------------------
plot_bar_abs_expression <- function(agi_list,
                                    gm,
                                    outdir   = "gene_metrics_out",
                                    file_stub = "AGI_bar_abs",
                                    y_transform = c("none", "log10"),
                                    pseudo_count = 0) {
  stopifnot(is.list(gm), !is.null(gm$avg_mat), !is.null(gm$groups))
  y_transform <- match.arg(y_transform)
  
  genes_all <- rownames(gm$avg_mat)
  agi_in    <- intersect(agi_list, genes_all)
  agi_miss  <- setdiff(agi_list, genes_all)
  
  if (length(agi_miss) > 0)
    warning("見つからない遺伝子: ", paste(agi_miss, collapse = ", "))
  if (length(agi_in) == 0)
    stop("描画可能な遺伝子がありません。")
  
  # long に整形
  avg_df <- as.data.frame(gm$avg_mat[agi_in, , drop = FALSE])
  avg_df$AGI <- rownames(avg_df)
  df <- reshape2::melt(avg_df, id.vars = "AGI",
                       variable.name = "Group", value.name = "MeanExpr")
  # 並びを固定
  df$Group <- factor(df$Group, levels = gm$groups)
  df$AGI   <- factor(df$AGI, levels = agi_in)
  
  # プロット本体
  p <- ggplot(df, aes(x = Group, y = MeanExpr, fill = AGI)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.minor = element_blank()) +
    labs(title = sprintf("Absolute expression (%s)", gm$grouping_desc),
         x = "Group", y = "Mean expression (log-normalized)")
  
  # 任意：ログ軸
  if (y_transform == "log10") {
    p <- p +
      aes(y = MeanExpr + pseudo_count) +
      scale_y_continuous(trans = "log10") +
      labs(y = sprintf("Mean expression (log10, +%g)", pseudo_count))
  }
  
  # 保存
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  W_PER_GROUP <- 0.45
  MIN_W <- 6.5
  pdf_w <- max(MIN_W, length(gm$groups) * W_PER_GROUP * 1.6)
  pdf_h <- 5.0
  out_pdf <- file.path(outdir, paste0(file_stub, ".pdf"))
  ggsave(out_pdf, p, width = pdf_w, height = pdf_h, units = "in", limitsize = FALSE)
  message("Saved PDF: ", normalizePath(out_pdf))
  
  invisible(list(pdf = out_pdf, data = df,
                 genes_drawn = agi_in, genes_missing = agi_miss))
}
# 1) 読み込み
x <- readRDS(RDS_PATH)

# 2) Gene metrics
gm <- build_gene_metrics(x, celltype_col = celltype_col, stage_col = stage_col)

# 3) 出力フォルダ
outdir <- file.path(dirname(RDS_PATH), "gene_metrics_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 4) 対象遺伝子（例）
agi_list <- c("AT2G20570","AT2G43010","AT1G75080")

# 5) 絶対値の棒グラフ（通常軸）
res_abs <- plot_bar_abs_expression(agi_list, gm,
                                   outdir = outdir,
                                   file_stub = "AGI_bar_abs",
                                   y_transform = "none")



# GLK1-PIF4-BIL1("AT2G20570","AT2G43010","AT1G75080")
# BIL7 familiy("AT1G63720","AT5G52430","AT4G25620","AT1G76660")

