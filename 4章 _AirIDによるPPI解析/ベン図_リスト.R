# ============================================================
# 完全版：
# AirID(暖色4条件) × CoIP(寒色4条件) の「交差集合(prot_acc)」を
# 16通り (4×4) で Excel に出力（prot_acc + SYMBOL列を併記）
#
# 仕様：
#  - 入力：2つのxlsx（AirID / CoIP）それぞれ4シート
#  - 出力：Excel
#     * summary_counts   : 交差数の行列
#     * summary_jaccard  : Jaccardの行列（0〜1）
#     * pair_...         : 16シート、prot_acc + SYMBOL_AirID + SYMBOL_CoIP
#  - SYMBOL列名が完全一致しない場合に備え、候補列を自動検出
#    （SYMBOL/symbol/gene っぽい列名を優先）
#  - 空白列名や末尾スペース等の “見えない事故” に耐性あり
# ============================================================

IN_XLSX_AIRID <- "/Users/ryo/Downloads/251220_AirID.xlsx"
SHEETS_AIRID  <- c("DMSO","bikinin","mock","BAK1")

IN_XLSX_COIP  <- "/Users/ryo/Downloads/251108_CoIP.xlsx"
SHEETS_COIP   <- c("DMSO","bikinin","BR-up","BR-down")

OUT_XLSX <- "/Users/ryo/Downloads/AirID_vs_CoIP_intersections_with_SYMBOL.xlsx"

# ---- packages ----
pkgs <- c("readxl", "stringr", "openxlsx", "dplyr")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(readxl)
library(stringr)
library(openxlsx)
library(dplyr)

# ---- helpers ----

# 列名の空白などを除去しつつシートをdata.frameで読む
read_sheet_df <- function(xlsx, sheet) {
  df <- readxl::read_excel(xlsx, sheet = sheet)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  # 列名の空白除去（SYMBOL などの一致事故を防ぐ）
  names(df) <- stringr::str_trim(names(df))
  
  df
}

# ベクトル長を変えない正規化（NAは残す）
norm_keep_len <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x[x == ""] <- NA_character_
  x
}

# prot_acc のユニーク集合
read_prot_acc_set <- function(df, prot_col = "prot_acc") {
  names(df) <- stringr::str_trim(names(df))
  if (!prot_col %in% names(df)) stop(sprintf("Column '%s' not found", prot_col))
  x <- norm_keep_len(df[[prot_col]])
  x <- x[!is.na(x)]
  unique(x)
}

# SYMBOL列名を自動推定
guess_symbol_col <- function(df, prefer = "SYMBOL") {
  nms <- stringr::str_trim(names(df))
  names(df) <- nms
  
  # まずは完全一致優先（大文字小文字も許容）
  if (prefer %in% nms) return(prefer)
  hit <- nms[tolower(nms) == tolower(prefer)]
  if (length(hit) >= 1) return(hit[1])
  
  # 次に候補（symbol/gene/name/locus など）
  cand <- nms[stringr::str_detect(tolower(nms), "symbol|gene|name|locus|tair")]
  if (length(cand) >= 1) return(cand[1])
  
  NA_character_
}

# prot_acc -> SYMBOL 対応表（同一prot_accが複数SYMBOLなら "A;B" で連結）
make_id_symbol_map <- function(df, prot_col = "prot_acc", sym_col = "SYMBOL") {
  names(df) <- stringr::str_trim(names(df))
  if (!prot_col %in% names(df)) stop(sprintf("Column '%s' not found", prot_col))
  
  sym_use <- sym_col
  if (!sym_use %in% names(df)) {
    sym_use <- guess_symbol_col(df, prefer = sym_col)
  }
  
  prot <- norm_keep_len(df[[prot_col]])
  sym  <- if (!is.na(sym_use)) norm_keep_len(df[[sym_use]]) else rep(NA_character_, length(prot))
  
  tmp <- data.frame(
    prot_acc = prot,
    SYMBOL   = sym,
    stringsAsFactors = FALSE
  )
  
  # prot_acc無し行は落とす（SYMBOLはNAでもOK）
  tmp <- tmp[!is.na(tmp$prot_acc), , drop = FALSE]
  
  tmp %>%
    group_by(prot_acc) %>%
    summarise(
      SYMBOL = {
        s <- SYMBOL[!is.na(SYMBOL)]
        if (length(s) == 0) NA_character_ else paste(sort(unique(s)), collapse = ";")
      },
      .groups = "drop"
    )
}

# Excelシート名安全化（禁止文字除去 + 31文字制限）
sanitize_sheet <- function(x) {
  x <- stringr::str_replace_all(x, "[\\[\\]\\:\\*\\?\\/\\\\]", "_")
  if (nchar(x) > 31) x <- substr(x, 1, 31)
  x
}

# Jaccard
jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  uni   <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

# ---- 1) read all sheets ----
airid_df_list <- setNames(
  lapply(SHEETS_AIRID, function(sh) read_sheet_df(IN_XLSX_AIRID, sh)),
  paste0("AirID_", SHEETS_AIRID)
)
coip_df_list <- setNames(
  lapply(SHEETS_COIP, function(sh) read_sheet_df(IN_XLSX_COIP, sh)),
  paste0("CoIP_", SHEETS_COIP)
)

# ---- 2) sets + maps ----
sets_airid <- lapply(airid_df_list, read_prot_acc_set)
sets_coip  <- lapply(coip_df_list,  read_prot_acc_set)

map_airid <- lapply(airid_df_list, make_id_symbol_map)
map_coip  <- lapply(coip_df_list,  make_id_symbol_map)

warm_names <- names(sets_airid)
cool_names <- names(sets_coip)

# ---- 3) summary matrices ----
count_mat <- matrix(0L, nrow = length(warm_names), ncol = length(cool_names),
                    dimnames = list(warm_names, cool_names))
jacc_mat  <- matrix(NA_real_, nrow = length(warm_names), ncol = length(cool_names),
                    dimnames = list(warm_names, cool_names))

# ---- 4) workbook ----
wb <- createWorkbook()
addWorksheet(wb, "summary_counts")
addWorksheet(wb, "summary_jaccard")

# ---- 5) 16 pair sheets ----
for (wn in warm_names) {
  for (cn in cool_names) {
    
    inter_vec <- sort(intersect(sets_airid[[wn]], sets_coip[[cn]]))
    count_mat[wn, cn] <- length(inter_vec)
    jacc_mat[wn, cn]  <- jaccard(sets_airid[[wn]], sets_coip[[cn]])
    
    ma <- map_airid[[wn]]
    mc <- map_coip[[cn]]
    
    # もし念のため SYMBOL 列が無い（ありえないが保険）
    if (!"SYMBOL" %in% names(ma)) ma$SYMBOL <- NA_character_
    if (!"SYMBOL" %in% names(mc)) mc$SYMBOL <- NA_character_
    
    ma <- dplyr::rename(ma, SYMBOL_AirID = SYMBOL)
    mc <- dplyr::rename(mc, SYMBOL_CoIP  = SYMBOL)
    
    df_out <- data.frame(prot_acc = inter_vec, stringsAsFactors = FALSE) %>%
      dplyr::left_join(ma, by = "prot_acc") %>%
      dplyr::left_join(mc, by = "prot_acc")
    
    sh <- sanitize_sheet(paste0("pair_", wn, "_x_", cn))
    addWorksheet(wb, sh)
    writeData(wb, sh, df_out)
    setColWidths(wb, sh, cols = 1:ncol(df_out), widths = "auto")
  }
}

# ---- 6) write summaries ----
writeData(wb, "summary_counts", as.data.frame(count_mat), rowNames = TRUE)
writeData(wb, "summary_jaccard", as.data.frame(round(jacc_mat, 3)), rowNames = TRUE)

# ---- 7) style ----
style_header <- createStyle(textDecoration = "bold", halign = "center")
style_int    <- createStyle(numFmt = "0")
style_dec    <- createStyle(numFmt = "0.000")

# counts
addStyle(wb, "summary_counts", style_header,
         rows = 1, cols = 1:(ncol(count_mat)+1), gridExpand = TRUE, stack = TRUE)
addStyle(wb, "summary_counts", style_int,
         rows = 2:(nrow(count_mat)+1), cols = 2:(ncol(count_mat)+1), gridExpand = TRUE, stack = TRUE)
setColWidths(wb, "summary_counts", cols = 1:(ncol(count_mat)+1), widths = "auto")

# jaccard
addStyle(wb, "summary_jaccard", style_header,
         rows = 1, cols = 1:(ncol(jacc_mat)+1), gridExpand = TRUE, stack = TRUE)
addStyle(wb, "summary_jaccard", style_dec,
         rows = 2:(nrow(jacc_mat)+1), cols = 2:(ncol(jacc_mat)+1), gridExpand = TRUE, stack = TRUE)
setColWidths(wb, "summary_jaccard", cols = 1:(ncol(jacc_mat)+1), widths = "auto")

# ---- 8) save ----
saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
cat("Saved:", OUT_XLSX, "\n")











