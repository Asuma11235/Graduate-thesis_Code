#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

# =========================
# 入力設定
# =========================
inputs <- list(
  leaf   = "/Users/ryo/Desktop/scAtlas/S3_leaf_rename.slim.RDS",
  stem   = "/Users/ryo/Desktop/scAtlas/stem_rename.slim.RDS",
  root   = "/Users/ryo/Desktop/scAtlas/D11_root_rename.slim.RDS",
  flower = "/Users/ryo/Desktop/scAtlas/Middle_flower_rename.slim.RDS"
)
TARGET_GENE <- "AT1G63720"

OUTDIR <- file.path("/Users/ryo/Desktop/scAtlas", "pseudobulk_out")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
OUTPDF <- file.path(OUTDIR, paste0(TARGET_GENE, "_pseudoBulk_CPM_bar.pdf"))

# =========================
# 共通ユーティリティ
# =========================
`%||%` <- function(a, b) if (!is.null(a)) a else b
nnz <- function(m) tryCatch(Matrix::nnzero(m), error = function(e) NA_integer_)

# v5/v4 互換ラッパー
GA <- function(x, assay = NULL, layer = c("counts","data","scale.data")) {
  layer <- match.arg(layer)
  ok <- try({
    return(GetAssayData(x, assay = assay %||% DefaultAssay(x), layer = layer))
  }, silent = TRUE)
  if (!inherits(ok, "try-error")) return(ok)
  slot <- switch(layer, counts="counts", data="data", "scale.data"="scale.data")
  GetAssayData(x, assay = assay %||% DefaultAssay(x), slot = slot)
}

choose_best_assay <- function(x, prefer = c("RNA","SCT","integrated","Spatial")) {
  asy <- Assays(x)
  ord <- unique(c(intersect(prefer, asy), setdiff(asy, prefer)))
  for (a in ord) {
    cnt <- tryCatch(GA(x, assay = a, layer = "counts"), error = function(e) NULL)
    dat <- tryCatch(GA(x, assay = a, layer = "data"),   error = function(e) NULL)
    if (!is.null(cnt) && !is.na(nnz(cnt)) && nnz(cnt) > 0) return(a)
    if (!is.null(dat) && !is.na(nnz(dat)) && nnz(dat) > 0) return(a)
  }
  ord[1]
}

# =========================
# 純バルク風CPM計算（1オブジェクト）
# =========================
pseudobulk_cpm_for_gene <- function(obj, gene_id, prefer_assay = c("RNA","SCT","integrated","Spatial")) {
  stopifnot(inherits(obj, "Seurat"))
  DefaultAssay(obj) <- choose_best_assay(obj, prefer = prefer_assay)
  
  # gene行の特定（大文字合わせ）
  rn <- rownames(obj)
  gid <- toupper(gene_id)
  if (!(gid %in% toupper(rn))) {
    stop(sprintf("対象遺伝子が見つかりません: %s", gene_id))
  }
  # 実際の行名（大文字照合→元の表記へ）
  gene_rowname <- rn[match(gid, toupper(rn))]
  
  # まず counts を使う（真のバルク風）
  cnt <- GA(obj, layer = "counts")
  if (!is.null(cnt) && !is.na(nnz(cnt)) && nnz(cnt) > 0) {
    # 合算ライブラリサイズ
    lib_size <- sum(cnt)
    if (lib_size <= 0) stop("counts は取得できたが合計が0です。")
    gene_sum <- sum(cnt[gene_rowname, , drop = TRUE])
    cpm <- 1e6 * as.numeric(gene_sum) / as.numeric(lib_size)
    return(list(cpm = cpm, source = "counts"))
  }
  
  # フォールバック：data（LogNormalize想定）を近似的にcountへ戻す
  warning("[fallback] counts が空のため data から近似CPMを計算します（expm1でカウント近似）")
  dat <- GA(obj, layer = "data")
  if (is.null(dat) || is.na(nnz(dat)) || nnz(dat) == 0) {
    stop("counts も data も空です。")
  }
  # data は log1p(1e4 * count / lib_size_cell) などのため、厳密な逆変換はできないが、
  # 可視化目的の近似として exmp1 を用いる
  # 注意: ここでの lib_size_approx は全要素のexpm1合計
  gene_sum <- sum(expm1(dat[gene_rowname, , drop = TRUE]))
  lib_size <- sum(expm1(dat))
  if (lib_size <= 0) stop("data近似の合計が0です。")
  cpm <- 1e6 * as.numeric(gene_sum) / as.numeric(lib_size)
  return(list(cpm = cpm, source = "data(approx)"))
}

# =========================
# 全データを回してCPM表を作成
# =========================
vals <- lapply(names(inputs), function(series) {
  p <- inputs[[series]]
  obj <- readRDS(p)
  res <- pseudobulk_cpm_for_gene(obj, TARGET_GENE)
  data.frame(Series = series, Gene = TARGET_GENE,
             CPM = res$cpm, Source = res$source, stringsAsFactors = FALSE)
})
df <- do.call(rbind, vals)

# 数値をログ化した見やすい軸にしたい場合は以下の列を使う（今回は生CPMを描画）
# df$log1pCPM <- log1p(df$CPM)

# =========================
# プロット & 保存
# =========================
p <- ggplot(df, aes(x = Series, y = CPM, fill = Series)) +
  geom_col(width = 0.7) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 15, hjust = 1)) +
  labs(
    title = sprintf("Pseudo-bulk CPM of %s", TARGET_GENE),
    x = "Series (dataset)",
    y = "CPM (counts per million)",
    caption = "Counts summed across all cells per dataset, then library-size normalized (CPM)."
  )

ggsave(OUTPDF, p, width = 6.5, height = 4.5, units = "in", limitsize = FALSE)

message("=== Pseudo-bulk CPM ===")
print(df)
message("Saved PDF: ", normalizePath(OUTPDF))

