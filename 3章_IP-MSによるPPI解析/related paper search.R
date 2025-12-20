# ===========================================
# Safe literature count add-on (+ BR-context counts)
# 入力: Protein_Zscore_and_Log2FC_fromZ.xlsx（TAIR/SYMBOL列あり推奨）
# 出力: GenePriority_with_literature.xlsx
# ===========================================
suppressPackageStartupMessages({
  if (!requireNamespace("rentrez", quietly=TRUE)) install.packages("rentrez")
  library(readxl); library(dplyr); library(stringr); library(tidyr)
  library(writexl); library(rentrez); library(scales)
  library(org.At.tair.db); library(AnnotationDbi)
})

# --- NCBI API key (to avoid rate limits) ---
Sys.setenv(ENTREZ_KEY = "c0533edd325919c40ae8b1a6a3ba48956e08")

IN_XLSX  <- "/Users/ryo/Downloads/251027_CoIP.xlsx"
SHEET_FC <- "Log2FC_from_Z_bySeries"
OUT_XLSX <- "GenePriority_with_literature.xlsx"

# fromZの列名（必要なら調整）
COL_S1 <- "Log2FC_Z_S1_BIL7GFP_DMSO_vs_Col0_DMSO"
COL_S2 <- "Log2FC_Z_S2_BIL7GFP_bikinin_vs_Col0_bikinin"
COL_S3 <- "Log2FC_Z_S3_BIL7GFP_bikinin_vs_BIL7GFP_DMSO"

# 1) 候補表
df <- read_excel(IN_XLSX, sheet = SHEET_FC) %>%
  mutate(across(any_of(c(COL_S1,COL_S2,COL_S3)), as.numeric))

if (!"TAIR"   %in% names(df)) df$TAIR   <- NA_character_
if (!"SYMBOL" %in% names(df)) df$SYMBOL <- NA_character_

# 2) TAIR -> SYMBOL を org.At.tair.db から可能な範囲で補完（ALIASは使わない）
avail_cols <- AnnotationDbi::columns(org.At.tair.db)
if ("SYMBOL" %in% avail_cols) {
  tair_keys <- unique(na.omit(df$TAIR))
  if (length(tair_keys) > 0) {
    sym_tab <- AnnotationDbi::select(org.At.tair.db,
                                     keys = tair_keys,
                                     columns = "SYMBOL",
                                     keytype = "TAIR") |>
      dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") |>
      dplyr::group_by(TAIR) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup()
    df <- dplyr::left_join(df, sym_tab, by = "TAIR", suffix = c("", "_DB"))
    df$SYMBOL <- dplyr::coalesce(df$SYMBOL, df$SYMBOL_DB)
    df$SYMBOL_DB <- NULL
  }
}

# 3) PubMed クエリ作成（一般 & BR文脈）
MAKE_TERM <- function(sym, tair, organism = TRUE, since_year = 2000){
  keys <- unique(na.omit(c(tair, sym)))
  if (length(keys) == 0) return(NA_character_)
  keys <- keys[nchar(keys) >= 3 & nchar(keys) <= 20]
  ors  <- paste(sprintf('("%s"[All Fields])', keys), collapse = " OR ")
  ctx  <- if (organism) '(arabidopsis[Title/Abstract] OR plant[Title/Abstract])' else ''
  years <- if (!is.null(since_year)) sprintf("(%d:3000[dp])", since_year) else ''
  paste(ors, ctx, years)
}

# ★ ご指定の BR 文脈検索式（括弧ごと固定文字列）
BR_CONTEXT <- paste0(
  '(',
  '(',
  '"Brassinosteroids"[MeSH] OR brassinosteroid*[tiab] OR brassinolide[tiab] OR castasterone[tiab] ',
  'OR BZR1[tiab] OR BES1[tiab] OR BIL1[tiab] OR BIN2[tiab] OR BRI1[tiab] OR BAK1[tiab] OR BSU1[tiab] OR PP2A[tiab]',
  ') AND (',
  '"Signal Transduction"[MeSH] OR signaling[tiab] OR pathway[tiab] OR cascade[tiab] OR phosphorylation[tiab] OR "protein kinase"[tiab]',
  ')',
  ') OR (',
  '(',
  '"Brassinosteroids"[MeSH] OR brassinosteroid*[tiab]',
  ') AND (',
  'auxin[tiab] OR cytokinin[tiab] OR gibberellin*[tiab] OR "abscisic acid"[tiab] OR jasmonic acid[tiab] OR "salicylic acid"[tiab] ',
  'OR ethylene[tiab] OR strigolactone*[tiab] OR karrikin*[tiab] OR PIF[tiab] OR "PHYTOCHROME INTERACTING FACTOR"[tiab] ',
  'OR "light signaling"[tiab] OR photomorphogenesis[tiab] OR "abiotic stress"[tiab] OR drought[tiab] OR salt[tiab] OR salinity[tiab] ',
  'OR cold[tiab] OR freezing[tiab] OR heat[tiab] OR "temperature stress"[tiab]',
  ') AND (',
  '"crosstalk"[tiab] OR interaction[tiab] OR synergy[tiab] OR antagonism[tiab]',
  ')',
  ')'
)

MAKE_TERM_BR <- function(sym, tair, organism = TRUE, since_year = 2000){
  base <- MAKE_TERM(sym, tair, organism = organism, since_year = since_year)
  if (is.na(base) || base == "") return(NA_character_)
  paste0('(', base, ') AND (', BR_CONTEXT, ')')
}

get_pubmed_count <- function(term){
  if (is.na(term) || term == "") return(0L)
  tryCatch({
    res <- rentrez::entrez_search(db = "pubmed", term = term, retmax = 0)
    as.integer(res$count)
  }, error = function(e) 0L)
}

# 4) 件数取得（一般 & BR 文脈）
since_year <- 2000   # 必要なら変更可
counts_general <- integer(nrow(df))
counts_br      <- integer(nrow(df))

for (i in seq_len(nrow(df))) {
  term_g  <- MAKE_TERM(sym = df$SYMBOL[i], tair = df$TAIR[i], organism = TRUE, since_year = since_year)
  term_br <- MAKE_TERM_BR(sym = df$SYMBOL[i], tair = df$TAIR[i], organism = TRUE, since_year = since_year)
  counts_general[i] <- get_pubmed_count(term_g)
  counts_br[i]      <- get_pubmed_count(term_br)
  if ((i %% 5) == 0) Sys.sleep(0.4)  # rate limit care
}

df$literature_count      <- counts_general
df$br_literature_count   <- counts_br
df$literature_score      <- rescale(log1p(df$literature_count), to = c(0,1))
df$br_literature_score   <- rescale(log1p(df$br_literature_count), to = c(0,1))

# 5) 効果量/一貫性の簡易スコア（以前の定義を流用）
score_strength <- function(x){
  x <- abs(x); x[is.na(x)] <- 0
  if (all(x == 0)) return(rep(0, length(x)))
  (x - min(x)) / (max(x) - min(x))
}
sign_consistency <- function(s1, s2, s3){
  s1s <- sign(s1); s2s <- sign(s2); s3s <- sign(s3)
  same12 <- (s1s != 0 & s1s == s2s)
  same_with_s3 <- (same12 & s3s == s1s)
  as.numeric(same12) + as.numeric(same_with_s3) # 0–2
}

df <- df %>%
  mutate(
    strength_S1    = score_strength(.data[[COL_S1]]),
    strength_S2    = score_strength(.data[[COL_S2]]),
    strength_S3    = score_strength(.data[[COL_S3]]),
    strength_combo = 0.4*strength_S1 + 0.4*strength_S2 + 0.6*strength_S3,
    consistency    = sign_consistency(.data[[COL_S1]], .data[[COL_S2]], .data[[COL_S3]])
  )

# GOスコアが未結合ならゼロに（後で結合して再計算もOK）
if (!all(c("GO_score_max","GO_score_sum") %in% names(df))){
  df$GO_score_max <- 0; df$GO_score_sum <- 0
}

# 6) 総合スコア（例）— BR 文脈スコアを加味したい場合は W_LITER_BR を追加
W_EFFECT = 0.47; W_CONSIS = 0.18; W_GO_MAX = 0.20; W_GO_SUM = 0.05
W_LITER  = 0.06  # 一般文献
W_LITER_BR = 0.04 # BR文脈（合計で従来の0.10に相当）

df$score <- 
  W_EFFECT*df$strength_combo +
  W_CONSIS*(df$consistency/2) +
  W_GO_MAX*rescale(df$GO_score_max, to=c(0,1), from=range(df$GO_score_max, na.rm=TRUE)) +
  W_GO_SUM*rescale(df$GO_score_sum, to=c(0,1), from=range(df$GO_score_sum, na.rm=TRUE)) +
  W_LITER   *df$literature_score +
  W_LITER_BR*df$br_literature_score

# 7) 出力
df_out <- df %>%
  arrange(desc(score)) %>%
  select(TAIR, SYMBOL, starts_with("Log2FC_Z_"),
         literature_count, literature_score,
         br_literature_count, br_literature_score,
         strength_S1, strength_S2, strength_S3, strength_combo, consistency,
         GO_score_max, GO_score_sum, score)

write_xlsx(list(Priority = df_out), OUT_XLSX)
message("出力: ", normalizePath(OUT_XLSX))

