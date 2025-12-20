#!/usr/bin/env Rscript

# ============================================================
# Safe Literature Count Add-on (+ Context-Specific Counts)
#
# Description:
#   Calculates priority scores based on experimental data (Z-scores)
#   and literature mining (PubMed counts).
#   Includes specific context filtering for Brassinosteroid (BR) research.
#
# Input:  Excel file with Z-scores and Log2FC (Recommended columns: TAIR, SYMBOL)
# Output: Excel file with added literature counts and priority scores
# ============================================================

suppressPackageStartupMessages({
  if (!requireNamespace("rentrez", quietly=TRUE)) install.packages("rentrez")
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(writexl)
  library(rentrez)
  library(scales)
  library(org.At.tair.db)
  library(AnnotationDbi)
})

# ---------------------------
# 1) Configuration
# ---------------------------

# --- NCBI API key ---
# NOTE: Set your API key to avoid strict rate limits (3 requests/sec vs 10/sec)
Sys.setenv(ENTREZ_KEY = "YOUR_NCBI_API_KEY") 

# Paths
IN_XLSX   <- "./input_data/Protein_Zscore_Data.xlsx"  # Generic input path
SHEET_FC  <- "Log2FC_from_Z_bySeries"
OUT_XLSX  <- "GenePriority_with_Literature_Scores.xlsx"

# Columns in the input file (Adjust as necessary)
COL_S1 <- "Log2FC_Z_S1_BIL7GFP_DMSO_vs_Col0_DMSO"
COL_S2 <- "Log2FC_Z_S2_BIL7GFP_bikinin_vs_Col0_bikinin"
COL_S3 <- "Log2FC_Z_S3_BIL7GFP_bikinin_vs_BIL7GFP_DMSO"

# Scoring Weights (Total should ideally be 1.0 or balanced)
W_EFFECT   <- 0.47
W_CONSIS   <- 0.18
W_GO_MAX   <- 0.20
W_GO_SUM   <- 0.05
W_LITER    <- 0.06  # General literature
W_LITER_BR <- 0.04  # Context-specific (BR) literature

# PubMed Search Settings
SINCE_YEAR <- 2000

# ---------------------------
# 2) Helper Functions
# ---------------------------

# Construct PubMed query for general search
build_pubmed_query <- function(sym, tair, organism = TRUE, since_year = 2000){
  keys <- unique(na.omit(c(tair, sym)))
  if (length(keys) == 0) return(NA_character_)
  
  # Filter short/long keywords to avoid noise
  keys <- keys[nchar(keys) >= 3 & nchar(keys) <= 20]
  if (length(keys) == 0) return(NA_character_)
  
  ors  <- paste(sprintf('("%s"[All Fields])', keys), collapse = " OR ")
  ctx  <- if (organism) '(arabidopsis[Title/Abstract] OR plant[Title/Abstract])' else ''
  years <- if (!is.null(since_year)) sprintf("(%d:3000[dp])", since_year) else ''
  
  paste(ors, ctx, years)
}

# Define Brassinosteroid (BR) Context Terms
# (Includes synonyms, signaling components, and crosstalk factors)
BR_CONTEXT_TERMS <- paste0(
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

# Construct PubMed query with BR context
build_br_context_query <- function(sym, tair, organism = TRUE, since_year = 2000){
  base <- build_pubmed_query(sym, tair, organism = organism, since_year = since_year)
  if (is.na(base) || base == "") return(NA_character_)
  paste0('(', base, ') AND (', BR_CONTEXT_TERMS, ')')
}

# Execute search and retrieve count
get_pubmed_count <- function(term){
  if (is.na(term) || term == "") return(0L)
  tryCatch({
    res <- rentrez::entrez_search(db = "pubmed", term = term, retmax = 0)
    as.integer(res$count)
  }, error = function(e) 0L)
}

# Calculate strength score (Normalized absolute value)
calc_strength_score <- function(x){
  x <- abs(x); x[is.na(x)] <- 0
  if (all(x == 0)) return(rep(0, length(x)))
  (x - min(x)) / (max(x) - min(x))
}

# Calculate consistency score (0 to 2)
# Checks if signs match across S1/S2 and S1/S3
calc_consistency_score <- function(s1, s2, s3){
  s1s <- sign(s1); s2s <- sign(s2); s3s <- sign(s3)
  same12 <- (s1s != 0 & s1s == s2s)
  same_with_s3 <- (same12 & s3s == s1s)
  as.numeric(same12) + as.numeric(same_with_s3)
}

# ---------------------------
# 3) Main Processing
# ---------------------------

# -- Load Data --
message("Loading data from: ", IN_XLSX)
df <- read_excel(IN_XLSX, sheet = SHEET_FC) %>%
  mutate(across(any_of(c(COL_S1, COL_S2, COL_S3)), as.numeric))

if (!"TAIR"   %in% names(df)) df$TAIR   <- NA_character_
if (!"SYMBOL" %in% names(df)) df$SYMBOL <- NA_character_

# -- Map TAIR to SYMBOL (fill missing using org.At.tair.db) --
avail_cols <- AnnotationDbi::columns(org.At.tair.db)
if ("SYMBOL" %in% avail_cols) {
  tair_keys <- unique(na.omit(df$TAIR))
  if (length(tair_keys) > 0) {
    sym_tab <- AnnotationDbi::select(org.At.tair.db,
                                     keys = tair_keys,
                                     columns = "SYMBOL",
                                     keytype = "TAIR") %>%
      dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
      dplyr::group_by(TAIR) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup()
    
    df <- dplyr::left_join(df, sym_tab, by = "TAIR", suffix = c("", "_DB"))
    df$SYMBOL <- dplyr::coalesce(df$SYMBOL, df$SYMBOL_DB)
    df$SYMBOL_DB <- NULL
  }
}

# -- Fetch PubMed Counts --
message("Starting PubMed queries...")
counts_general <- integer(nrow(df))
counts_br      <- integer(nrow(df))

for (i in seq_len(nrow(df))) {
  term_g  <- build_pubmed_query(sym = df$SYMBOL[i], tair = df$TAIR[i], organism = TRUE, since_year = SINCE_YEAR)
  term_br <- build_br_context_query(sym = df$SYMBOL[i], tair = df$TAIR[i], organism = TRUE, since_year = SINCE_YEAR)
  
  counts_general[i] <- get_pubmed_count(term_g)
  counts_br[i]      <- get_pubmed_count(term_br)
  
  # Rate limit safety (approx. 2-3 requests per sec with sleep)
  if ((i %% 5) == 0) Sys.sleep(0.4) 
  
  if (i %% 100 == 0) message(sprintf("Processed %d / %d genes...", i, nrow(df)))
}

df$literature_count    <- counts_general
df$br_literature_count <- counts_br

# -- Normalize Counts --
df$literature_score    <- rescale(log1p(df$literature_count), to = c(0,1))
df$br_literature_score <- rescale(log1p(df$br_literature_count), to = c(0,1))

# -- Calculate Signal Scores --
df <- df %>%
  mutate(
    strength_S1    = calc_strength_score(.data[[COL_S1]]),
    strength_S2    = calc_strength_score(.data[[COL_S2]]),
    strength_S3    = calc_strength_score(.data[[COL_S3]]),
    strength_combo = 0.4 * strength_S1 + 0.4 * strength_S2 + 0.6 * strength_S3,
    consistency    = calc_consistency_score(.data[[COL_S1]], .data[[COL_S2]], .data[[COL_S3]])
  )

# Initialize GO scores if missing
if (!all(c("GO_score_max", "GO_score_sum") %in% names(df))){
  df$GO_score_max <- 0; df$GO_score_sum <- 0
}

# -- Calculate Final Priority Score --
df$score <- 
  W_EFFECT * df$strength_combo +
  W_CONSIS * (df$consistency / 2) +
  W_GO_MAX * rescale(df$GO_score_max, to=c(0,1), from=range(df$GO_score_max, na.rm=TRUE)) +
  W_GO_SUM * rescale(df$GO_score_sum, to=c(0,1), from=range(df$GO_score_sum, na.rm=TRUE)) +
  W_LITER      * df$literature_score +
  W_LITER_BR   * df$br_literature_score

# -- Output --
df_out <- df %>%
  arrange(desc(score)) %>%
  select(TAIR, SYMBOL, starts_with("Log2FC_Z_"),
         literature_count, literature_score,
         br_literature_count, br_literature_score,
         strength_S1, strength_S2, strength_S3, strength_combo, consistency,
         GO_score_max, GO_score_sum, score)

write_xlsx(list(Priority = df_out), OUT_XLSX)
message("Output saved to: ", normalizePath(OUT_XLSX))