


# BIL7 familiy("AT1G63720","AT5G52430","AT4G25620","AT1G76660")
# AmPR familiy("AT4G30850","AT2G24150","AT5G20270","AT4G37680","AT4G38320")
# ナツキさん関連("AT2G36270","AT3G56850","AT3G44460","AT2G41070","AT2G35530","AT1G32150")


ABI5: AT2G36270
AREB3: AT3G56850
DPBF2: AT3G44460
DPBF4: AT2G41070
bZIP16: AT2G35530
bZIP68: AT1G32150





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
RDS_PATH     <- "/Users/ryo/Downloads/D5_silique_rename.slim.RDS"
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
# B) AGI ベクタから描画 & PDF 保存する関数
#    ※ カラーバーは plasma（低=紫→高=黄）
# -----------------------------------------
plot_and_save_from_AGI <- function(agi_list, gm, outdir = "gene_metrics_out", file_stub = "AGI_batch") {
  
  stopifnot(is.list(gm), !is.null(gm$avg_mat), !is.null(gm$pct_mat), !is.null(gm$groups))
  genes_all <- rownames(gm$avg_mat)
  agi_in <- intersect(agi_list, genes_all)
  agi_miss <- setdiff(agi_list, agi_in)
  if (length(agi_miss) > 0) {
    warning("見つからない遺伝子: ", paste(agi_miss, collapse = ", "))
  }
  if (length(agi_in) == 0) stop("描画可能な遺伝子がありません。")
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # 長い形に整形
  avg_df <- as.data.frame(gm$avg_mat[agi_in, , drop = FALSE])
  avg_df$AGI <- rownames(avg_df)
  avg_long <- melt(avg_df, id.vars = "AGI", variable.name = "Group", value.name = "MeanExpr")
  
  pct_df <- as.data.frame(gm$pct_mat[agi_in, , drop = FALSE])
  pct_df$AGI <- rownames(pct_df)
  pct_long <- melt(pct_df, id.vars = "AGI", variable.name = "Group", value.name = "Pct")
  
  df <- merge(avg_long, pct_long, by = c("AGI", "Group"))
  # 並びを固定
  df$Group <- factor(df$Group, levels = gm$groups)
  df$AGI   <- factor(df$AGI, levels = rev(agi_in))  # 上から指定順に
  
  # ドット図（色＝平均発現: plasma、サイズ＝検出率）
  p <- ggplot(df, aes(x = Group, y = AGI)) +
    geom_point(aes(size = Pct, color = MeanExpr)) +
    scale_size_continuous(name = "Detection rate", range = c(0, 7), limits = c(0, 1)) +
    scale_color_viridis_c(name = "Mean expr", option = "plasma") +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = sprintf("Gene metrics (%s)", gm$grouping_desc),
      x = "Group", y = "AGI",
      caption = "Color: avg log-normalized expression (plasma; low=purple → high=yellow) / Size: detection rate (counts>0)"
    )
  
  # PDF サイズ自動
  W_PER_GROUP <- 0.35
  H_PER_GENE  <- 0.45
  MIN_W <- 6.5
  MIN_H <- 5.0
  pdf_w <- max(MIN_W, length(gm$groups) * W_PER_GROUP * 1.6)
  pdf_h <- max(MIN_H, length(agi_in)    * H_PER_GENE  * 1.6)
  
  dot_pdf <- file.path(outdir, paste0(file_stub, "_dot.pdf"))
  ggsave(dot_pdf, p, width = pdf_w, height = pdf_h, units = "in", limitsize = FALSE)
  
  message("Saved PDF: ", normalizePath(dot_pdf))
  
  invisible(list(
    pdf = dot_pdf,
    data_long = df,
    genes_drawn = agi_in,
    genes_missing = agi_miss
  ))
}


# =========================
# 実行例（あなたのフォーマット通り）
# =========================

# 1) .RDS 読み込み
x <- readRDS(RDS_PATH)

# 2) 例A: celltype だけで集計
gm <- build_gene_metrics(x, celltype_col = celltype_col, stage_col = stage_col)

# 3) .RDSと同じ場所に保存用フォルダを作る
outdir <- dirname(RDS_PATH)
outdir <- file.path(outdir, "gene_metrics_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 4) 用意した AGI ベクタでまとめて描画・保存
agi_list <- c("AT2G36270","AT3G56850","AT3G44460","AT2G41070","AT2G35530","AT1G32150")
res <- plot_and_save_from_AGI(agi_list, gm, outdir = outdir, file_stub = "AGI_batch")

message("出力先: ", normalizePath(outdir))




# BIL7 familiy("AT1G63720","AT5G52430","AT4G25620","AT1G76660")
# AmPR familiy("AT4G30850","AT2G24150","AT5G20270","AT4G37680","AT4G38320")
# ナツキさん関連("AT2G36270","AT3G56850","AT3G44460","AT2G41070","AT2G35530","AT1G32150")

























# 全遺伝子 x 全細胞 でスケーリング→処理落ち
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
RDS_PATH     <- "/Users/ryo/Downloads/D0_silique_rename.slim.RDS"
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
  
  # --- 平均発現（slot='scale.data'）版。必要なら scale.data を自動生成 ---
  # ※ scale.data は正負あり（平均0に中心化済み）。解釈は z スコアの平均です。
  
  # 1) scale.data が無い／空／必要遺伝子が足りないときは埋める
  sd_mat <- tryCatch(GetAssayData(x, assay = DefaultAssay(x), slot = "scale.data"),
                     error = function(e) NULL)
  need_scale <- is.null(sd_mat) || nrow(sd_mat) == 0 ||
    !all(rownames(x) %in% rownames(sd_mat))
  
  if (need_scale) {
    message("[avg-from-scale] scale.data が不足 → ScaleData() を実行して埋めます（全遺伝子）")
    vf_backup <- VariableFeatures(x)
    # 全遺伝子を対象にスケーリング（重い場合は、必要遺伝子だけに変えてOK）
    x <- ScaleData(x, features = rownames(x), verbose = FALSE)
    # 可変遺伝子の設定は元に戻す（下流の解析に影響させないため）
    VariableFeatures(x) <- vf_backup
  }
  
  # 2) scale.data からグループ平均を計算
  avg_obj <- AverageExpression(
    x,
    group.by = ".__group__",
    assays   = DefaultAssay(x),
    slot     = "scale.data",
    verbose  = FALSE
  )
  
  # 3) 返り値の型を吸収して matrix に整形
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
  if (inherits(avg_mat, "Matrix") || inherits(avg_mat, "dgCMatrix"))  avg_mat <- as.matrix(avg_mat)
  if (!is.matrix(avg_mat)) stop("AverageExpression（scale.data）を matrix にできません: ",
                                paste(class(avg_mat), collapse=", "))
  
  # 4) 行名・列名のフォールバック
  if (is.null(rownames(avg_mat))) rownames(avg_mat) <- rownames(x)
  if (is.null(colnames(avg_mat))) colnames(avg_mat) <- groups
  
  
  # 検出率（counts>0）。counts が空なら data>0 を代用
  cnt <- GetAssayData(x, assay = DefaultAssay(x), slot = "counts")
  if (Matrix::nnzero(cnt) == 0) {
    message("counts slot is empty. Using data>0 as detection proxy.")
    dat <- GetAssayData(x, assay = DefaultAssay(x), slot = "scale.data")
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
# B) AGI ベクタから描画 & PDF 保存する関数
#    ※ カラーバーは plasma（低=紫→高=黄）
# -----------------------------------------
plot_and_save_from_AGI <- function(agi_list, gm, outdir = "gene_metrics_out", file_stub = "AGI_batch") {
  
  stopifnot(is.list(gm), !is.null(gm$avg_mat), !is.null(gm$pct_mat), !is.null(gm$groups))
  genes_all <- rownames(gm$avg_mat)
  agi_in <- intersect(agi_list, genes_all)
  agi_miss <- setdiff(agi_list, agi_in)
  if (length(agi_miss) > 0) {
    warning("見つからない遺伝子: ", paste(agi_miss, collapse = ", "))
  }
  if (length(agi_in) == 0) stop("描画可能な遺伝子がありません。")
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # 長い形に整形
  avg_df <- as.data.frame(gm$avg_mat[agi_in, , drop = FALSE])
  avg_df$AGI <- rownames(avg_df)
  avg_long <- melt(avg_df, id.vars = "AGI", variable.name = "Group", value.name = "MeanExpr")
  
  pct_df <- as.data.frame(gm$pct_mat[agi_in, , drop = FALSE])
  pct_df$AGI <- rownames(pct_df)
  pct_long <- melt(pct_df, id.vars = "AGI", variable.name = "Group", value.name = "Pct")
  
  df <- merge(avg_long, pct_long, by = c("AGI", "Group"))
  # 並びを固定
  df$Group <- factor(df$Group, levels = gm$groups)
  df$AGI   <- factor(df$AGI, levels = rev(agi_in))  # 上から指定順に
  
  # ドット図（色＝平均発現: plasma、サイズ＝検出率）
  p <- ggplot(df, aes(x = Group, y = AGI)) +
    geom_point(aes(size = Pct, color = MeanExpr)) +
    scale_size_continuous(name = "Detection rate", range = c(0, 7), limits = c(0, 1)) +
    scale_color_viridis_c(name = "Mean expr", option = "plasma") +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = sprintf("Gene metrics (%s)", gm$grouping_desc),
      x = "Group", y = "AGI",
      caption = "Color: avg log-normalized expression (plasma; low=purple → high=yellow) / Size: detection rate (counts>0)"
    )
  
  # PDF サイズ自動
  W_PER_GROUP <- 0.35
  H_PER_GENE  <- 0.45
  MIN_W <- 6.5
  MIN_H <- 5.0
  pdf_w <- max(MIN_W, length(gm$groups) * W_PER_GROUP * 1.6)
  pdf_h <- max(MIN_H, length(agi_in)    * H_PER_GENE  * 1.6)
  
  dot_pdf <- file.path(outdir, paste0(file_stub, "_dot.pdf"))
  ggsave(dot_pdf, p, width = pdf_w, height = pdf_h, units = "in", limitsize = FALSE)
  
  message("Saved PDF: ", normalizePath(dot_pdf))
  
  invisible(list(
    pdf = dot_pdf,
    data_long = df,
    genes_drawn = agi_in,
    genes_missing = agi_miss
  ))
}


# =========================
# 実行例（あなたのフォーマット通り）
# =========================

# 1) .RDS 読み込み
x <- readRDS(RDS_PATH)

# 2) 例A: celltype だけで集計
gm <- build_gene_metrics(x, celltype_col = celltype_col, stage_col = stage_col)

# 3) .RDSと同じ場所に保存用フォルダを作る
outdir <- dirname(RDS_PATH)
outdir <- file.path(outdir, "gene_metrics_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 4) 用意した AGI ベクタでまとめて描画・保存
agi_list <- c("AT2G36270","AT3G56850","AT3G44460","AT2G41070","AT2G35530","AT1G32150")
res <- plot_and_save_from_AGI(agi_list, gm, outdir = outdir, file_stub = "AGI_batch")

message("出力先: ", normalizePath(outdir))

