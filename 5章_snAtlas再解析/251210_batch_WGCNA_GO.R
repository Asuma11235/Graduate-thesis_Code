install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.At.tair.db"))







suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.At.tair.db)
  library(dplyr)
  library(ggplot2)
})

# =========================
# 0) ★設定★
# =========================

BASE_OUTDIR <- "/Users/ryo/Desktop/scAtlas"  # WGCNA 出力先
ORGDB       <- org.At.tair.db       # Arabidopsis 用 OrgDb
KEYTYPE     <- "TAIR"               # AGI を想定
ONTOLOGY    <- "BP"                 # "BP", "MF", "CC" から選択

# BIL7 とどの程度相関したモジュールを GO 対象にするか
COR_THRESHOLD <- 0.6   # |cor_BIL7| >= 0.6 のモジュールを対象にする（必要なら変更）

# GO の q-value 閾値（多重検定補正後）
QVALUE_CUTOFF <- 0.05

# =========================
# 1) モジュール情報 & 全遺伝子情報の読込関数
# =========================

read_wgcna_results <- function(organ) {
  organ_dir <- file.path(BASE_OUTDIR, organ)
  
  module_trait_file <- file.path(organ_dir,
                                 paste0("Module_trait_correlation_BIL7_",
                                        organ, ".csv"))
  all_genes_file    <- file.path(organ_dir,
                                 paste0("All_genes_module_and_kME_",
                                        organ, ".csv"))
  
  if (!file.exists(module_trait_file)) {
    stop("File not found: ", module_trait_file)
  }
  if (!file.exists(all_genes_file)) {
    stop("File not found: ", all_genes_file)
  }
  
  mt_df  <- read.csv(module_trait_file, stringsAsFactors = FALSE)
  ag_df  <- read.csv(all_genes_file, stringsAsFactors = FALSE)
  
  list(
    module_trait = mt_df,
    all_genes    = ag_df,
    organ_dir    = organ_dir
  )
}

# =========================
# 2) GO 解析（1モジュール分）関数
# =========================

run_go_for_module <- function(gene_vec,
                              universe_vec,
                              module_label,
                              organ,
                              outdir,
                              ont = ONTOLOGY,
                              q_cutoff = QVALUE_CUTOFF) {
  # NA や重複を除く
  gene_vec     <- unique(na.omit(gene_vec))
  universe_vec <- unique(na.omit(universe_vec))
  
  if (length(gene_vec) < 5) {
    warning("[", organ, "] Module ", module_label,
            " has fewer than 5 genes. Skipping GO.")
    return(NULL)
  }
  
  ego <- enrichGO(
    gene          = gene_vec,
    universe      = universe_vec,
    OrgDb         = ORGDB,
    keyType       = KEYTYPE,
    ont           = ont,
    pAdjustMethod = "BH",
    qvalueCutoff  = q_cutoff,
    readable      = TRUE   # 遺伝子名も付与
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    warning("[", organ, "] Module ", module_label,
            " has no significant GO terms.")
    return(NULL)
  }
  
  # 結果をデータフレームにして保存
  ego_df <- as.data.frame(ego)
  
  out_csv <- file.path(outdir,
                       paste0("GO_", ont, "_", organ, "_", module_label, ".csv"))
  write.csv(ego_df, out_csv, row.names = FALSE)
  message("[", organ, "] Saved GO table: ", out_csv)
  
  # 簡易 dotplot も保存（上位20件くらい）
  # dotplot（読みやすいフォントサイズに調整）
  dp <- dotplot(ego, showCategory = 20) +
    ggtitle(paste0("GO ", ont, " - ", organ, " - ", module_label)) +
    theme_bw() +
    theme(
      plot.title      = element_text(size = 10),
      axis.text.x     = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y     = element_text(size = 10),
      axis.title.x    = element_text(size = 9),
      axis.title.y    = element_text(size = 9),
      legend.title    = element_text(size = 8),
      legend.text     = element_text(size = 7)
    )
  
  out_pdf <- file.path(outdir,
                       paste0("GO_", ont, "_", organ, "_", module_label, "_dotplot.pdf"))
  ggsave(out_pdf, dp, width = 7, height = 6)   # height を少し増やすと重なりを防げる
  
  message("[", organ, "] Saved GO dotplot: ", out_pdf)
  
  invisible(ego)
}

# =========================
# 3) organ ごとに GO を実行
# =========================

run_go_for_organ <- function(organ) {
  message("========== [", organ, "] GO analysis ==========")
  dat <- read_wgcna_results(organ)
  mt_df <- dat$module_trait
  ag_df <- dat$all_genes
  organ_dir <- dat$organ_dir
  
  # 対象とするモジュール（|cor_BIL7| >= COR_THRESHOLD）
  target_modules <- mt_df %>%
    dplyr::filter(abs(cor_BIL7) >= COR_THRESHOLD) %>%
    dplyr::arrange(desc(abs(cor_BIL7)))
  
  # もし何も閾値を満たさなければ「最も相関が大きいモジュール」だけを解析
  if (nrow(target_modules) == 0) {
    warning("[", organ, "] No module passes |cor_BIL7| >= ", COR_THRESHOLD,
            ". Using only the top module.")
    target_modules <- mt_df[order(-abs(mt_df$cor_BIL7)), , drop = FALSE][1, , drop = FALSE]
  }
  
  message("[", organ, "] Target modules:")
  print(target_modules[, c("module", "module_color", "cor_BIL7", "p_BIL7")])
  
  # 全遺伝子（universe）
  universe_genes <- ag_df$gene
  
  # モジュールごとに GO
  for (i in seq_len(nrow(target_modules))) {
    mod_label <- target_modules$module[i]        # 例: "MEturquoise"
    mod_color <- target_modules$module_color[i]  # 例: "turquoise"
    
    message("[", organ, "] Running GO for module: ", mod_label,
            " (color=", mod_color, ")")
    
    genes_in_mod <- ag_df %>%
      dplyr::filter(module == mod_color & !is.na(gene)) %>%
      dplyr::pull(gene)
    
    if (length(genes_in_mod) == 0) {
      warning("[", organ, "] No genes found for module color ", mod_color)
      next
    }
    
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
  
  message("========== [", organ, "] GO analysis DONE ==========")
}

# 実行
run_go_for_organ("S3_leaf_GLK1")


message("=== All GO analyses finished. ===")

