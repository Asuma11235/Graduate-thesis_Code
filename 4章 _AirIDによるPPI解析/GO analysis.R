

#!/usr/bin/env Rscript

# ============================================================
# GO enrichment 可視化 拡張版（白背景・堅牢化・デバッグ付き）
#  - enrichGO + bar
#  - (A) GOplot: GOBubble（不能時は dotplot 代替）/ GOChord（堅牢化）
#  - (B) GO階層ツリー（sig内の親子） + （任意）enrichplot::goplot DAG
# ============================================================

# ---------------------------
# 0) Install & Load
# ---------------------------
install_if_missing <- function(pkgs, bioc=FALSE){
  for(p in pkgs){
    if (!requireNamespace(p, quietly=TRUE)){
      if (bioc){
        if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
        BiocManager::install(p, ask=FALSE, update=FALSE)
      } else install.packages(p)
    }
    suppressPackageStartupMessages(library(p, character.only=TRUE))
  }
}

install_if_missing(c(
  "readxl","dplyr","stringr","ggplot2","tidyr","purrr","tibble",
  "igraph","ggraph","ggforce","GOplot","circlize"
))
install_if_missing(c(
  "clusterProfiler","DOSE","enrichplot","org.At.tair.db","AnnotationDbi","GO.db"
), bioc=TRUE)

# 白背景を全体のデフォルトに
ggplot2::theme_set(ggplot2::theme_bw(base_size = 11))

# --- 保存を堅牢化するユーティリティ（ggsave フォールバック＋ログ） ---
safe_save_png <- function(path, plot, width=7, height=6, dpi=300, bg="white"){
  dir.create(dirname(path), showWarnings=FALSE, recursive=TRUE)
  ok <- TRUE
  # まず ggsave を試す（環境によりこれで十分）
  tryCatch({
    ggplot2::ggsave(filename = path, plot = plot, width = width, height = height, dpi = dpi,
                    device = function(...) grDevices::png(..., bg = bg))
  }, error=function(e){ ok <<- FALSE; message("[save] ggsave failed: ", e$message) })
  if (ok && file.exists(path)) {
    sz <- try(file.info(path)$size, silent=TRUE)
    message("[save] ggsave OK: ", path, " (", ifelse(inherits(sz,"try-error")||is.na(sz),"?", sz), " bytes)")
    return(invisible(TRUE))
  }
  # フォールバック：明示的にデバイスを開く
  message("[save] fallback device path: ", path)
  devok <- FALSE
  if (requireNamespace("Cairo", quietly=TRUE)) {
    try({
      Cairo::CairoPNG(filename = path, width = width*dpi, height = height*dpi, dpi = dpi, bg = bg)
      devok <- TRUE; message("[save] Cairo::CairoPNG open")
    }, silent=FALSE)
  }
  if (!devok) {
    try({
      grDevices::png(path, width = width*dpi, height = height*dpi, res = dpi, type="cairo-png", bg = bg)
      devok <- TRUE; message("[save] png(type='cairo-png') open")
    }, silent=FALSE)
  }
  if (!devok) {
    grDevices::png(path, width = width*dpi, height = height*dpi, res = dpi, bg = bg)
    devok <- TRUE; message("[save] png(default) open")
  }
  print(plot)
  grDevices::dev.off()
  if (file.exists(path)) {
    sz <- try(file.info(path)$size, silent=TRUE)
    message("[save] fallback OK: ", path, " (", ifelse(inherits(sz,"try-error")||is.na(sz),"?", sz), " bytes)")
  } else {
    message("[save] fallback FAILED: ", path)
  }
  invisible(TRUE)
}

ggsave_white <- function(filename, plot, width=7, height=6, dpi=300){
  # （既存関数は互換維持。内部は safe_save_png に委譲）
  safe_save_png(filename, plot, width=width, height=height, dpi=dpi, bg="white")
}

# ---------------------------
# 1) Settings
# ---------------------------
IN_XLSX   <- "/Users/ryo/Downloads/251108_CoIP.xlsx"
SHEETS    <- c("DMSO","bikinin","BR-up","BR-down", "Col-up", "Col-down")
OUT_DIR   <- "251108_xlsxglobal_final_!!"; dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

GO_ONTOLOGY  <- "CC"      # "BP" / "MF" / "CC"
PV_CUTOFF    <- 0.05
QV_CUTOFF    <- 0.05

USE_QUANTILE <- TRUE       # TRUE: 分位, FALSE: 固定閾値
TOP_PCT      <- 100
ABS_CUTOFF   <- 0.1

COL_SYM  <- "SYMBOL"
COL_LFC  <- "Log2FC"

DRAW_DAG <- TRUE   # enrichplot::goplot を出すか

# ★★ 新規: 背景選択モード ★★
# "sheet"          : 既存の挙動（当該シートでTAIRにマップできた全遺伝子）
# "global"         : org.At.tair.db に存在する全 TAIR を背景
# "excel_universe" : IN_XLSX 内 "universe" シートの SYMBOL 列から作る背景

BG_MODE <- "excel_universe"  # ← "sheet" / "global" / "excel_universe"

UNIVERSE_SHEET_NAME <- "universe"  # excel_universe モードで使うシート名
UNIVERSE_SYMBOL_COL <- "SYMBOL"    # universe シート内の列名

# ---------------------------
# 2) Helpers（※既存の関数は流用）
# ---------------------------
# ---- 可視化の自動スケール用（線形マップ＋クリップ） ----
.map_linear <- function(x, from_min, from_max, to_min, to_max) {
  if (is.na(x) || !is.finite(x)) return((to_min + to_max)/2)
  x <- max(from_min, min(from_max, x))
  r <- (x - from_min) / (from_max - from_min)
  to_min + r * (to_max - to_min)
}


clean_chr <- function(x){
  x <- as.character(x); x <- trimws(x); x[x==""] <- NA; x
}
extract_agi_anycase <- function(x){
  stringr::str_extract(toupper(as.character(x)), "AT[1-5MC]G\\d{5}")
}
map_symbol_to_tair_any <- function(sym_or_alias_vec){
  k <- unique(stats::na.omit(toupper(trimws(sym_or_alias_vec))))
  if (length(k) == 0) return(setNames(character(0), character(0)))
  kt <- AnnotationDbi::keytypes(org.At.tair.db)
  m <- AnnotationDbi::mapIds(org.At.tair.db, keys=k, keytype="SYMBOL",
                             column="TAIR", multiVals="first")
  need <- names(m)[is.na(m)]
  if (length(need) && "SYNONYM" %in% kt){
    m2 <- AnnotationDbi::mapIds(org.At.tair.db, keys=need, keytype="SYNONYM",
                                column="TAIR", multiVals="first")
    m[need] <- ifelse(is.na(m[need]), m2[need], m[need])
  }
  if (length(need) && "ALIAS" %in% kt){
    m3 <- AnnotationDbi::mapIds(org.At.tair.db, keys=need, keytype="ALIAS",
                                column="TAIR", multiVals="first")
    m[need] <- ifelse(is.na(m[need]), m3[need], m[need])
  }
  m
}

# enrichGO（軽いデバッグ出力のみ追加）
run_enrichGO_TAIR <- function(genes_tair, universe_tair, label){
  genes_tair <- unique(genes_tair[!is.na(genes_tair) & genes_tair!=""])
  message(sprintf("[%s] input=%d, universe=%d", label, length(genes_tair), length(universe_tair)))
  if (length(genes_tair) < 5L){
    message(sprintf("[%s] 入力遺伝子が少ないためスキップ (%d genes).", label, length(genes_tair)))
    return(NULL)
  }
  ego <- tryCatch({
    clusterProfiler::enrichGO(
      gene          = genes_tair,
      universe      = universe_tair,
      OrgDb         = org.At.tair.db,
      keyType       = "TAIR",
      ont           = GO_ONTOLOGY,
      pAdjustMethod = "BH",
      pvalueCutoff  = PV_CUTOFF,
      qvalueCutoff  = QV_CUTOFF,
      readable      = FALSE
    )
  }, error=function(e){
    message(sprintf("[%s] enrichGO error: %s", label, e$message)); NULL
  })
  if (is.null(ego) || nrow(as.data.frame(ego))==0){
    message(sprintf("[%s] 有意ヒットなし。", label)); return(NULL)
  }
  message(sprintf("[%s] enrichGO: %d rows", label, nrow(as.data.frame(ego))))
  ego
}

plot_bar_save <- function(ego, title, out_png){
  dfp <- as.data.frame(ego); if (nrow(dfp) == 0) return(invisible(NULL))
  dfp$GeneRatioNum <- suppressWarnings(sapply(dfp$GeneRatio, function(x) if (nzchar(x)) eval(parse(text=x)) else NA_real_))
  dfp <- dfp[!is.na(dfp$GeneRatioNum), , drop=FALSE]
  if (!nrow(dfp)) return(invisible(NULL))
  dfp <- dfp[order(-dfp$GeneRatioNum, dfp$p.adjust), , drop=FALSE]
  topn <- min(15, nrow(dfp)); if (topn < 1) return(invisible(NULL))
  df_top <- dfp[seq_len(topn), , drop=FALSE]
  df_top <- df_top[order(df_top$GeneRatioNum), , drop=FALSE]
  df_top$Description <- factor(df_top$Description, levels=df_top$Description)
  
  p <- ggplot2::ggplot(df_top,
                       ggplot2::aes(x=Description, y=GeneRatioNum, fill=p.adjust)) +
    ggplot2::geom_col(color="black") +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_gradient(low="#d95159", high="#6baee0", name="adj. P") +
    ggplot2::labs(y="GeneRatio", x=NULL, title=title) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(color="black", fill=NA)
    )
  safe_save_png(out_png, p, width=6, height=4, dpi=300)  # ← ggsave を堅牢保存に置換
}

# ---- GOplot 前処理（circle_dat 入力）
diagnose_circ <- function(terms, genes){
  g_in_terms <- unique(unlist(strsplit(terms$genes, ",\\s*")))
  hit <- intersect(g_in_terms, genes$ID)
  message(sprintf("[GOplot] terms側の遺伝子数: %d / genesテーブルID数: %d / 一致数: %d",
                  length(g_in_terms), nrow(genes), length(hit)))
}

prep_GOplot_inputs <- function(ego, de_df, ontology_tag){
  if (is.null(ego)) return(NULL)
  egodf <- as.data.frame(ego); if (!nrow(egodf)) return(NULL)
  
  terms <- egodf %>%
    dplyr::transmute(
      category = ontology_tag,
      ID       = .data$ID,
      term     = .data$Description,
      genes    = gsub("/", ", ", .data$geneID),   # TAIR をカンマ列挙
      adj_pval = .data$p.adjust
    )
  
  genes <- de_df %>%
    dplyr::mutate(ID = toupper(ID),
                  logFC = suppressWarnings(as.numeric(logFC))) %>%
    dplyr::filter(!is.na(ID), is.finite(logFC)) %>%
    dplyr::distinct(ID, .keep_all = TRUE)
  
  terms$genes <- vapply(strsplit(terms$genes, ",\\s*"), function(v){
    paste(toupper(v), collapse = ", ")
  }, FUN.VALUE = character(1))
  
  diagnose_circ(terms, genes)
  
  circ <- tryCatch(GOplot::circle_dat(terms, genes), error=function(e){ message("[circle_dat] ", e$message); NULL })
  if (is.null(circ) || !NROW(circ)) return(NULL)
  
  circ <- circ %>% dplyr::filter(!is.na(ID), !is.na(term))
  if (!NROW(circ)) return(NULL)
  
  list(terms=terms, genes=genes, circ=circ)
}

# ---- GOplot 出力（Bubbleは代替dotplot、Chordは circlize で安定出力＋ラベル付き）
save_GOplot_figs <- function(circ, out_stub, ego=NULL, de_df=NULL){
  # 代替Bubble
  save_alt_dotplot <- function(ego_obj, stub){
    if (is.null(ego_obj)) return(invisible(NULL))
    egodf <- as.data.frame(ego_obj); if (!NROW(egodf)) return(invisible(NULL))
    p_alt <- enrichplot::dotplot(ego_obj, showCategory = min(15, nrow(egodf))) +
      ggplot2::theme_bw(base_size=12)
    safe_save_png(paste0(stub, "_GOplot_Bubble_alt_dotplot.png"), p_alt, width=7, height=6, dpi=300)
    message("[GOplot] Bubble 代替: dotplot を出力 → ", paste0(stub, "_GOplot_Bubble_alt_dotplot.png"))
  }
  
  # ===== Bubble =====
  if (is.null(circ) || !NROW(circ)) {
    message("[GOplot] circ が空 → Bubble代替")
    save_alt_dotplot(ego, out_stub)
  } else {
    n_pos <- sum(suppressWarnings(as.numeric(circ$logFC)) > 0, na.rm=TRUE)
    n_neg <- sum(suppressWarnings(as.numeric(circ$logFC)) < 0, na.rm=TRUE)
    if (n_pos == 0 || n_neg == 0) {
      message(sprintf("[GOplot] 片側のみ (pos=%d, neg=%d) → Bubble代替", n_pos, n_neg))
      save_alt_dotplot(ego, out_stub)
    } else {
      lab_n <- max(0, min(5, nrow(circ)))
      p1 <- tryCatch(GOplot::GOBubble(circ, labels = lab_n) + ggplot2::theme_bw(12),
                     error=function(e){ message("[GOBubble] ", e$message); NULL })
      if (!is.null(p1)) {
        safe_save_png(paste0(out_stub, "_GOplot_Bubble.png"), p1, width=7, height=6, dpi=300)
      } else {
        save_alt_dotplot(ego, out_stub)
      }
    }
  }
  
  # ===== Chord（ID整合+頻度優先で sel を作る 安定版）=====
  if (is.null(circ) || !NROW(circ)) { message("[GOplot] circ 空 → Chord スキップ"); return(invisible(NULL)) }
  
  gene_col <- if ("genes" %in% names(circ)) "genes" else if ("gene" %in% names(circ)) "gene" else NA_character_
  if (is.na(gene_col)) { message("[GOplot] circ に gene 列が無い → Chord スキップ"); return(invisible(NULL)) }
  
  # 1) circ に現れる gene リスト（これを絶対基準にする）
  genes_in_circ <- unique(toupper(as.character(circ[[gene_col]])))
  genes_in_circ <- genes_in_circ[!is.na(genes_in_circ) & nzchar(genes_in_circ)]
  
  # 2) logFC テーブル（de_df があれば補完）。キーは "GENE" = 大文字揃え
  de_tbl <- NULL
  if (!is.null(de_df)) {
    de_tbl <- de_df %>%
      dplyr::transmute(GENE = toupper(as.character(ID)),
                       logFC = suppressWarnings(as.numeric(logFC))) %>%
      dplyr::filter(!is.na(GENE), is.finite(logFC)) %>%
      dplyr::distinct(GENE, .keep_all = TRUE)
  }
  
  # 3) gene の term 出現頻度
  freq_tbl <- circ %>%
    dplyr::mutate(GENE = toupper(as.character(.data[[gene_col]]))) %>%
    dplyr::filter(!is.na(GENE) & nzchar(GENE)) %>%
    dplyr::count(GENE, name = "freq")
  
  # 4) sel 候補
  sel_pool <- freq_tbl$GENE
  
  # 5) sel のスコア = (頻度, |logFC|)
  if (!is.null(de_tbl)) {
    sel_df <- freq_tbl %>%
      dplyr::left_join(de_tbl, by = c("GENE")) %>%
      dplyr::mutate(absLFC = abs(logFC))
  } else {
    sel_df <- freq_tbl %>% dplyr::mutate(absLFC = NA_real_)
  }
  sel_df <- sel_df %>% dplyr::arrange(dplyr::desc(freq), dplyr::desc(absLFC), GENE)
  sel <- utils::head(sel_df$GENE, 20L)
  if (length(sel) < 2L) { message("[GOplot] 候補遺伝子 < 2 → Chord スキップ"); return(invisible(NULL)) }
  
  # 6) chdat
  chdat <- tryCatch(GOplot::chord_dat(circ, genes = sel),
                    error=function(e){ message("[chord_dat] ", e$message); NULL })
  if (is.null(chdat) || !NROW(chdat)) { message("[GOplot] chord_dat が空 → Chord スキップ"); return(invisible(NULL)) }
  
  # 7) 行名・列名の安全取得
  rn <- rownames(chdat)
  if (!length(rn) || anyNA(rn) || any(rn=="")) {
    first_col <- if (ncol(chdat) >= 1) chdat[[1]] else NULL
    if (!is.null(first_col)) rn <- as.character(first_col)
  }
  if (!length(rn) || length(rn) != nrow(chdat)) rn <- paste0("gene_", seq_len(nrow(chdat)))
  cn <- colnames(chdat)
  term_cols <- if (ncol(chdat) >= 2) 2:ncol(chdat) else integer(0)
  if (!length(term_cols)) { message("[circlize] term 列が見つからない → スキップ"); return(invisible(NULL)) }
  tn <- cn[term_cols]; if (!length(tn)) tn <- paste0("term_", seq_along(term_cols))
  
  mat0 <- suppressWarnings(as.matrix(chdat[, term_cols, drop=FALSE]))
  storage.mode(mat0) <- "numeric"; mat0[is.na(mat0)] <- 0
  
  idx <- which(mat0 > 0, arr.ind = TRUE)
  if (!nrow(idx)) { message("[circlize] リンクなし（sel と term の交差=0） → sel の決め方を再検討"); return(invisible(NULL)) }
  
  edges_df <- data.frame(
    gene   = rn[idx[,"row"]],
    term   = tn[idx[,"col"]],
    weight = 1,
    stringsAsFactors = FALSE
  )
  message("[Chord-debug] genes_in_circ: ", length(genes_in_circ))
  message("[Chord-debug] sel (head): ", paste(head(sel), collapse=", "))
  message("[Chord-debug] chdat dim: ", paste(dim(chdat), collapse=" x "))
  message("[Chord-debug] terms: ", length(tn))
  message("[Chord-debug] edges: ", nrow(idx))
  
  # ラベルの用意
  genes_sectors <- unique(edges_df$gene)
  terms_sectors <- unique(edges_df$term)
  
  tair2sym <- tryCatch({
    AnnotationDbi::mapIds(org.At.tair.db,
                          keys = genes_sectors,
                          keytype = "TAIR",
                          column = "SYMBOL",
                          multiVals = "first")
  }, error = function(e) setNames(rep(NA_character_, length(genes_sectors)), genes_sectors))
  gene_labels <- ifelse(is.na(unname(tair2sym[genes_sectors])) | unname(tair2sym[genes_sectors])=="",
                        genes_sectors, unname(tair2sym[genes_sectors]))
  
  ego_df <- tryCatch(as.data.frame(ego), error=function(e) NULL)
  go_desc_map <- if (!is.null(ego_df) && "ID" %in% names(ego_df) && "Description" %in% names(ego_df)) {
    setNames(as.character(ego_df$Description), ego_df$ID)
  } else setNames(rep(NA_character_, length(terms_sectors)), terms_sectors)
  term_labels0 <- unname(go_desc_map[terms_sectors])
  term_labels  <- ifelse(is.na(term_labels0) | term_labels0=="", terms_sectors, term_labels0)
  
  .wrap <- function(x, width=18){
    vapply(x, function(s){
      s <- as.character(s); if (!nzchar(s)) return(s)
      out <- character()
      while (nchar(s) > width) {
        cutpos <- max(gregexpr("[ \\-_/]", substr(s, 1, width))[[1]])
        if (is.finite(cutpos) && cutpos > 0) {
          out <- c(out, substr(s, 1, cutpos))
          s <- trimws(substr(s, cutpos + 1, nchar(s)))
        } else {
          out <- c(out, substr(s, 1, width))
          s <- substr(s, width + 1, nchar(s))
        }
      }
      paste(c(out, s), collapse = "\n")
    }, character(1))
  }
  gene_labels_wrapped <- .wrap(gene_labels, width=16)
  term_labels_wrapped <- .wrap(term_labels, width=22)
  
  sectors <- c(genes_sectors, terms_sectors)
  label_map <- c(setNames(gene_labels_wrapped, genes_sectors),
                 setNames(term_labels_wrapped, terms_sectors))
  
  # === circlize描画：デバイス堅牢化＋詳細ログ（既存のまま） ===
  f_png2 <- paste0(out_stub, "_GOplot_Chord_alt_circlize.png")
  message(sprintf("[circlize] sectors: genes=%d / terms=%d; edges=%d",
                  length(unique(edges_df$gene)), length(unique(edges_df$term)), nrow(edges_df)))
  try(circlize::circos.clear(), silent = TRUE)
  
  open_device <- function(path) {
    ok <- FALSE
    if (requireNamespace("Cairo", quietly = TRUE)) {
      try({
        Cairo::CairoPNG(filename = path, width = 9*300, height = 6*300, dpi = 300)
        ok <<- TRUE
        message("[device] Cairo::CairoPNG でオープン")
      }, silent = FALSE)
      if (ok) return(invisible(TRUE))
    }
    try({
      grDevices::png(path, width = 9*300, height = 6*300, res = 300, type = "cairo-png")
      ok <<- TRUE
      message("[device] png(type='cairo-png') でオープン")
    }, silent = FALSE)
    if (ok) return(invisible(TRUE))
    grDevices::png(path, width = 9*300, height = 6*300, res = 300)
    message("[device] png(通常) でオープン")
    invisible(TRUE)
  }
  
  open_device(f_png2)
  ok_draw <- TRUE
  try({
    circlize::circos.clear()
    # ★ キャンバス余白と余裕を広げる
    circlize::circos.par(
      gap.after = c(rep(2, length(unique(edges_df$gene)) - 1), 10,
                    rep(2, length(unique(edges_df$term)) - 1), 10),
      track.margin = c(0.01, 0.04),                 # ★ 外側に少し広め（上側の余白）
      cell.padding = c(0, 0, 0, 0),                 # ★ 余白の“謎の食い込み”防止
      canvas.xlim = c(-1.06, 1.06),                 # ★ 左右に少しだけ拡張
      canvas.ylim = c(-1.08, 1.12),                 # ★ 上側を広めに拡張（見切れ対策の本丸）
      points.overflow.warning = FALSE               # ★ クリップ警告を抑制（挙動安定）
    )
    
    circlize::chordDiagram(
      x = edges_df[, c("gene","term","weight")],
      order = c(unique(edges_df$gene), unique(edges_df$term)),
      directional = 0,
      annotationTrack = "grid",
      transparency = 0.2,
      link.sort = TRUE,
      link.largest.ontop = TRUE,
      preAllocateTracks = list(track.height = 0.18) # ★ ラベル用トラックを少し厚めに
    )
    
    # ★ テキストを“外縁から2mm内側”に押し込んで描く（= 物理長さ基準）
    circlize::circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      sector_name <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      ylim <- circlize::get.cell.meta.data("ylim")
      lab <- label_map[[sector_name]]; if (is.null(lab)) lab <- sector_name
      
      # 2mm をデータ座標へ変換（track.index=1 / 現在セルに対して）
      inset <- circlize::convert_y(9.5, "mm", sector.index = sector_name, track.index = 1)
      
      circlize::circos.text(
        x = mean(xlim),
        y = ylim[2] - inset,                         # ★ 外縁ぴったり→内側へ 2mm
        labels = lab,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),                             # ★ ベースライン中央寄せ
        cex = 0.6
      )
    })
    circlize::circos.clear()
  }, silent = FALSE) -> res_draw
  grDevices::dev.off()
  
  
  if (!file.exists(f_png2)) {
    message("[circlize] 出力PNGが存在しません: ", f_png2)
  } else {
    sz <- try(file.info(f_png2)$size, silent = TRUE)
    message("[circlize] 出力ファイル: ", f_png2, " / size=", if (inherits(sz, "try-error")) "NA" else as.character(sz), " bytes")
    if (!inherits(res_draw, "try-error") && !isTRUE(inherits(sz, "try-error")) && !is.na(sz) && sz > 2048) {
      message("[circlize] OK: PNGサイズ十分（>2KB）")
    } else {
      message("[circlize] 小さすぎ/描画エラーの可能性 → 上記ログのエラー/警告を確認してください")
    }
  }
  invisible(NULL)
}  # save_GOplot_figs()

# ---- (B) GO階層ツリー（sig内の親子のみでエッジ）
get_GO_parents <- function(go_ids, ontology=c("BP","MF","CC")){
  ontology <- match.arg(ontology)
  parents_map <- switch(ontology,
                        BP = GO.db::GOBPPARENTS,
                        MF = GO.db::GOMFPARENTS,
                        CC = GO.db::GOCCPARENTS)
  out <- setNames(vector("list", length(go_ids)), go_ids)
  for (i in seq_along(go_ids)){
    ch <- go_ids[i]
    p <- suppressWarnings(AnnotationDbi::mget(ch, parents_map, ifnotfound=NA)[[1]])
    if (all(is.na(p))) out[[i]] <- character(0)
    else out[[i]] <- unique(unname(as.character(p)))
  }
  out
}

plot_go_sig_tree <- function(ego, title, out_png, ontology_tag){
  egodf <- as.data.frame(ego); if (!nrow(egodf)) return(invisible(NULL))
  sig_ids <- unique(egodf$ID); if (length(sig_ids) < 2L){ message("[tree] sig<2"); return(invisible(NULL)) }
  
  # 親子（sig内のみ）
  parents_list <- get_GO_parents(sig_ids, ontology=ontology_tag)
  edges <- purrr::map2(sig_ids, parents_list, function(child, parents){
    if (!length(parents)) return(NULL)
    parents_in_sig <- intersect(parents, sig_ids)
    if (!length(parents_in_sig)) return(NULL)
    tibble::tibble(from=parents_in_sig, to=child)
  }) %>% dplyr::bind_rows()
  if (is.null(edges) || !nrow(edges)){ message("[tree] 親子なし"); return(invisible(NULL)) }
  
  nodes <- egodf %>%
    dplyr::select(ID, Description, p.adjust, GeneRatio) %>%
    dplyr::distinct() %>%
    dplyr::mutate(GeneRatioNum = sapply(GeneRatio, function(x) if (nzchar(x)) eval(parse(text=x)) else NA_real_)) %>%
    dplyr::rename(name=ID)
  
  # ---- 量に応じた自動スケール ----
  n_nodes <- nrow(nodes)
  n_edges <- nrow(edges)
  
  lab_size   <- .map_linear(n_nodes,  10, 200, 12.0, 5.5)   # ノード数が多いほど文字小さく
  node_size0 <- .map_linear(n_nodes,  1, 1000, 0.4, 0.4)   # 点自体の最小サイズ
  edge_width <- .map_linear(n_edges,  20,2000, 2.4,0.5)  # エッジ本数が多いほど細く
  edge_alpha <- .map_linear(n_edges,  20,2000, 0.80,0.25)  # 同上：透明度上げる
  box_size   <- .map_linear(n_nodes,  10, 200, 0.30,0.15)  # ラベルボックス枠線
  repel_force<- .map_linear(n_nodes,  10, 200,  1.0, 2.5)  # 混雑時に少し強めに弾く
  
  g <- igraph::graph_from_data_frame(edges, directed=TRUE, vertices=nodes)
  
  p <- ggraph::ggraph(g, layout="tree") +
    # エッジ（量に応じて幅/透明度調整）
    ggraph::geom_edge_link(
      linewidth = edge_width, alpha = edge_alpha,
      arrow = grid::arrow(length=unit(3,"mm"), type="closed"),
      end_cap = ggraph::circle(2.5,'mm'), colour = "grey30"
    ) +
    # ノード（点の下限サイズを与えて見失わないように）
    ggraph::geom_node_point(
      ggplot2::aes(size = pmax(GeneRatioNum, node_size0), alpha = -log10(p.adjust)),
      show.legend = FALSE
    ) +
    # ラベル（フォントを量で可変。重なりは抑制しつつ極力残す）
    ggraph::geom_node_label(
      ggplot2::aes(label = Description),
      size = lab_size, label.size = box_size, label.padding = unit(0.10, "lines"),
      fill = "white", colour = "black", segment.colour = "grey40",
      # ggrepel由来の引数：混雑に応じて調整
      force = repel_force, max.overlaps = Inf
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(
      title = title,
      subtitle = paste0("GO ", ontology_tag, " (sig terms only, parent-child)")
    )
  
  safe_save_png(out_png, p, width=6, height=4, dpi=300)
}


plot_goplot_DAG <- function(ego, title, out_png){
  if (is.null(ego)) return(invisible(NULL))
  egodf <- as.data.frame(ego); if (!nrow(egodf)) return(invisible(NULL))
  
  # ノード数＝有意語数（goplotはここが主な混雑要因）
  n_nodes <- nrow(egodf)
  base_text <- .map_linear(n_nodes, 10, 200, 5.5, 5.5)     # 図全体の文字
  lab_size  <- .map_linear(n_nodes, 10, 200, 1.8, 2.0)  # ノードラベル
  edge_w    <- .map_linear(n_nodes, 10, 200, 1.0, 0.3)  # エッジ線幅（近似）
  edge_a    <- .map_linear(n_nodes, 10, 200, 0.8, 0.3)  # エッジ透明度
  node_size0 <- .map_linear(n_nodes,  1, 1000, 0.4, 0.4)   # 点自体の最小サイズ
  
  p <- tryCatch({
    enrichplot::goplot(ego, showCategory = min(30, n_nodes))
  }, error=function(e){ message("[goplot] ", e$message); NULL })
  if (is.null(p)) return(invisible(NULL))
  
  # 可能なら edge 幅/透明度のスケールを付与（無視される場合もあり）
  suppressWarnings({
    p <- p +
      ggplot2::theme_bw(base_size = base_text) +
      ggplot2::theme(
        plot.title  = ggplot2::element_text(size = base_text + 1),
        legend.text = ggplot2::element_text(size = base_text - 1),
        legend.title= ggplot2::element_text(size = base_text)
      )
    
    # 一部の enrichplot バージョンでは edge_* aesthetics が有効
    if (requireNamespace("ggraph", quietly = TRUE)) {
      p <- p +
        ggraph::scale_edge_width(range = c(edge_w, edge_w), guide = "none") +
        ggraph::scale_edge_alpha(range = c(edge_a, edge_a), guide = "none")
    }
    
    # --- ノード点サイズの自動調整（n_nodes に応じて） ---
    node_min <- .map_linear(n_nodes, 10, 200, 3.5, 1.2)
    node_max <- .map_linear(n_nodes, 10, 200, 7.0, 2.0)
    
    # size が「美術マッピング」されているレイヤがあるか判定（aes_evalは使わない）
    has_size_mapping <- any(vapply(p$layers, function(lyr) {
      nm <- tryCatch(names(lyr$mapping), error = function(e) character(0))
      "size" %in% nm
    }, logical(1)))
    
    if (has_size_mapping) {
      # マッピングされているなら scale_size でレンジを制御
      p <- p + ggplot2::scale_size(range = c(node_min, node_max), guide = "none")
    } else {
      # 固定サイズなら該当レイヤの size を直接上書き
      for (i in seq_along(p$layers)) {
        lyr <- p$layers[[i]]
        if (inherits(lyr$geom, "GeomPoint") || inherits(lyr$geom, "GeomNodePoint") || inherits(lyr$geom, "GeomCircle")) {
          lyr$aes_params$size  <- node_min
          lyr$geom_params$size <- node_min
          p$layers[[i]] <- lyr
        }
      }
    }
    
    # --- ラベル層：未対応引数を渡さない（警告回避） ---
    # 既存の goplot 出力内の Label/Text レイヤがあればサイズだけ上書き
    for (i in seq_along(p$layers)) {
      lyr <- p$layers[[i]]
      if (inherits(lyr$geom, "GeomLabelRepel") || inherits(lyr$geom, "GeomTextRepel") ||
          inherits(lyr$geom, "GeomLabel")      || inherits(lyr$geom, "GeomText")) {
        lyr$aes_params$size  <- lab_size
        lyr$geom_params$size <- lab_size
        # 未対応のパラメータは設定しない
        p$layers[[i]] <- lyr
      }
    }
  })
  
  
  safe_save_png(out_png, p + ggplot2::ggtitle(title), width=6, height=4, dpi=300)
}



# ---------------------------
# 背景（universe）を作る関数
# ---------------------------
build_universe_from_mode <- function(mode = BG_MODE,
                                     in_xlsx = IN_XLSX,
                                     universe_sheet = UNIVERSE_SHEET_NAME,
                                     universe_symbol_col = UNIVERSE_SYMBOL_COL,
                                     df_sheet_tair = NULL){
  mode <- match.arg(mode, c("sheet","global","excel_universe"))
  
  if (mode == "sheet"){
    if (is.null(df_sheet_tair) || !length(df_sheet_tair))
      stop("[BG] sheet モードですが df_sheet_tair が空です。")
    bg <- unique(toupper(df_sheet_tair))
    all_tair <- AnnotationDbi::keys(org.At.tair.db, keytype="TAIR")
    return(intersect(bg, all_tair))
  }
  
  if (mode == "global"){
    bg <- AnnotationDbi::keys(org.At.tair.db, keytype="TAIR")
    return(unique(toupper(bg)))
  }
  
  if (mode == "excel_universe"){
    df_u <- tryCatch(readxl::read_excel(in_xlsx, sheet = universe_sheet),
                     error = function(e){ stop(sprintf("[BG] universe シート読込失敗: %s", e$message)) })
    if (!(universe_symbol_col %in% names(df_u)))
      stop(sprintf("[BG] universe シートに '%s' 列が見つかりません。", universe_symbol_col))
    v <- clean_chr(df_u[[universe_symbol_col]])
    if (!length(v) || all(is.na(v))) stop("[BG] universe シートの SYMBOL 列が空です。")
    tair_from_agi <- extract_agi_anycase(v)
    need_map <- is.na(tair_from_agi)
    sym_keys <- toupper(trimws(v[need_map]))
    sym2tair  <- map_symbol_to_tair_any(sym_keys)
    tair_from_sym <- rep(NA_character_, length(v))
    if (length(sym2tair) > 0){
      tair_from_sym[need_map] <- unname(sym2tair[toupper(v[need_map])])
    }
    bg <- unique(na.omit(toupper(dplyr::coalesce(tair_from_agi, tair_from_sym))))
    if (!length(bg)) stop("[BG] universe シートから TAIR を1つも作れませんでした。")
    all_tair <- AnnotationDbi::keys(org.At.tair.db, keytype="TAIR")
    bg <- intersect(bg, all_tair)
    if (!length(bg)) stop("[BG] universe シートの TAIR が org.At に存在しません。")
    return(bg)
  }
}

# ---------------------------
# 3) Main
# ---------------------------
for (sh in SHEETS){
  message("\n==============================\nSheet: ", sh, "\n==============================")
  
  df_raw <- tryCatch(readxl::read_excel(IN_XLSX, sheet=sh),
                     error=function(e){ warning(sprintf("[%s] 読み込み失敗: %s", sh, e$message)); NULL })
  if (is.null(df_raw)) next
  
  needed <- c(COL_SYM, COL_LFC)
  miss <- setdiff(needed, names(df_raw))
  if (length(miss) > 0){ warning(sprintf("[%s] 必要列欠如: %s", sh, paste(miss, collapse=", "))); next }
  
  df <- df_raw |>
    dplyr::mutate(!!COL_SYM := clean_chr(.data[[COL_SYM]]),
                  !!COL_LFC := suppressWarnings(as.numeric(.data[[COL_LFC]]))) |>
    dplyr::filter(!is.na(.data[[COL_SYM]]))
  if (!nrow(df)){ warning(sprintf("[%s] SYMBOL 有効行なし", sh)); next }
  
  # --- TAIR マッピング（SYMBOL列にTAIRが入っているケースも対応） ---
  df$TAIR_from_AGI <- extract_agi_anycase(df[[COL_SYM]])
  need_map <- is.na(df$TAIR_from_AGI)
  sym_keys <- toupper(trimws(df[[COL_SYM]][need_map]))
  sym2tair <- map_symbol_to_tair_any(sym_keys)
  df$TAIR_from_SYM <- NA_character_
  if (length(sym2tair) > 0){
    df$TAIR_from_SYM[need_map] <- unname(sym2tair[toupper(df[[COL_SYM]][need_map])])
  }
  df$TAIR <- dplyr::coalesce(df$TAIR_from_AGI, df$TAIR_from_SYM)
  
  map_audit <- df |>
    dplyr::transmute(INPUT_SYMBOL = .data[[COL_SYM]],
                     Log2FC = .data[[COL_LFC]],
                     TAIR_from_AGI, TAIR_from_SYM,
                     TAIR_final = TAIR)
  
  df <- df |>
    dplyr::filter(!is.na(TAIR) & TAIR!="") |>
    dplyr::arrange(dplyr::desc(abs(.data[[COL_LFC]]))) |>
    dplyr::distinct(TAIR, .keep_all = TRUE)
  
  out_dir_sheet <- file.path(OUT_DIR, sh); dir.create(out_dir_sheet, showWarnings=FALSE, recursive=TRUE)
  if (nrow(df) < 10){
    warning(sprintf("[%s] TAIR マップ数が少ない: %d", sh, nrow(df)))
    utils::write.csv(map_audit, file=file.path(out_dir_sheet, "mapping_audit.csv"), row.names=FALSE)
    next
  }
  message(sprintf("[%s] TAIR mapped: %d", sh, nrow(df)))
  
  # 背景（universe）
  bg_universe <- build_universe_from_mode(
    mode = BG_MODE,
    in_xlsx = IN_XLSX,
    universe_sheet = UNIVERSE_SHEET_NAME,
    universe_symbol_col = UNIVERSE_SYMBOL_COL,
    df_sheet_tair = df$TAIR
  )
  message(sprintf("[BG] mode=%s | universe size=%d", BG_MODE, length(bg_universe)))
  
  # Up/Down gene set（背景でフィルタ）
  if (USE_QUANTILE){
    uq <- as.numeric(stats::quantile(df[[COL_LFC]], probs=1 - TOP_PCT/100, na.rm=TRUE))
    dq <- as.numeric(stats::quantile(df[[COL_LFC]], probs=TOP_PCT/100,     na.rm=TRUE))
    up_genes   <- df$TAIR[!is.na(df[[COL_LFC]]) & df[[COL_LFC]] >=  uq]
    down_genes <- df$TAIR[!is.na(df[[COL_LFC]]) & df[[COL_LFC]] <=  dq]
  } else {
    up_genes   <- df$TAIR[!is.na(df[[COL_LFC]]) & df[[COL_LFC]] >=  ABS_CUTOFF]
    down_genes <- df$TAIR[!is.na(df[[COL_LFC]]) & df[[COL_LFC]] <= -ABS_CUTOFF]
  }
  up_genes   <- intersect(unique(toupper(up_genes)),   bg_universe)
  down_genes <- intersect(unique(toupper(down_genes)), bg_universe)
  message(sprintf("[%s] gene sets: Up:%d | Down:%d", sh, length(up_genes), length(down_genes)))
  
  # 監査ログ
  utils::write.csv(map_audit, file=file.path(out_dir_sheet, "mapping_audit.csv"), row.names=FALSE)
  utils::write.csv(data.frame(TAIR=up_genes),   file=file.path(out_dir_sheet, "UP_genes.csv"),   row.names=FALSE)
  utils::write.csv(data.frame(TAIR=down_genes), file=file.path(out_dir_sheet, "DOWN_genes.csv"), row.names=FALSE)
  
  # --- enrichGO 実行 ---
  ego_up   <- run_enrichGO_TAIR(up_genes,   bg_universe, paste0(sh, "_UP"))
  ego_down <- run_enrichGO_TAIR(down_genes, bg_universe, paste0(sh, "_DOWN"))
  
  # ---- UP 出力（★ループ内★）----
  if (!is.null(ego_up)){
    utils::write.csv(as.data.frame(ego_up), file=file.path(out_dir_sheet, paste0("UP_GO_", GO_ONTOLOGY, ".csv")), row.names=FALSE)
    plot_bar_save(ego_up, paste0(sh," (UP) : GO ", GO_ONTOLOGY, " enrichment"),
                  file.path(out_dir_sheet, paste0("UP_GO_", GO_ONTOLOGY, "_bar.png")))
    de_up <- df %>% dplyr::select(ID=TAIR, logFC=!!rlang::sym(COL_LFC)) %>% dplyr::mutate(ID=toupper(ID))
    gp_up <- prep_GOplot_inputs(ego_up, de_up, GO_ONTOLOGY)
    if (!is.null(gp_up)){
      save_GOplot_figs(gp_up$circ, file.path(out_dir_sheet, paste0("UP_", GO_ONTOLOGY)),
                       ego=ego_up, de_df=de_up)
      plot_go_sig_tree(ego_up, paste0(sh," (UP) sig-only tree"),
                       file.path(out_dir_sheet, paste0("UP_GO_", GO_ONTOLOGY, "_tree.png")), GO_ONTOLOGY)
      if (DRAW_DAG){
        plot_goplot_DAG(ego_up, paste0(sh," (UP) goplot DAG"),
                        file.path(out_dir_sheet, paste0("UP_GO_", GO_ONTOLOGY, "_DAG.png")))
      }
    }
  }
  
  # ---- DOWN 出力（★ループ内★）----

  
  # このシートで作られたファイル一覧を必ずログ
  created <- try(list.files(out_dir_sheet, full.names = TRUE), silent=TRUE)
  if (!inherits(created, "try-error")) {
    message("[files] ", out_dir_sheet, " :\n  - ", paste(basename(created), collapse = "\n  - "))
  }
}

cat("完了：", normalizePath(OUT_DIR), "\n")

