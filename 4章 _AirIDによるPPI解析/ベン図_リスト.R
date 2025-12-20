#!/usr/bin/env Rscript

# ============================================================
# Pairwise Intersection Analysis: AirID (4 sets) x CoIP (4 sets)
#
# Description:
#   Calculates intersections of 'prot_acc' between two datasets (AirID and CoIP),
#   each containing multiple conditions (sheets).
#   Generates a comprehensive Excel report containing:
#    1. Summary Matrix: Intersection counts
#    2. Summary Matrix: Jaccard indices
#    3. Detailed Lists (16 sheets): Intersected IDs with SYMBOLs from both sources.
#
# Features:
#   - Automatic detection of "SYMBOL" columns (handling case/variations like 'gene', 'name').
#   - Robust string normalization (trimming whitespace).
#   - formatted Excel output with auto-width columns.
# ============================================================

# ----------------------------
# 1) Configuration
# ----------------------------

# Input Files
IN_XLSX_AIRID <- "./input_data/AirID_Dataset.xlsx"
SHEETS_AIRID  <- c("DMSO", "bikinin", "mock", "BAK1")

IN_XLSX_COIP  <- "./input_data/CoIP_Dataset.xlsx"
SHEETS_COIP   <- c("DMSO", "bikinin", "BR-up", "BR-down")

# Output File
OUT_XLSX      <- "./results/AirID_vs_CoIP_Intersections.xlsx"
dir.create(dirname(OUT_XLSX), showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 2) Libraries
# ----------------------------
pkgs <- c("readxl", "stringr", "openxlsx", "dplyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# ----------------------------
# 3) Helper Functions
# ----------------------------

# Read sheet as data.frame, trimming whitespace from column names
read_clean_sheet <- function(xlsx, sheet) {
  df <- tryCatch(
    readxl::read_excel(xlsx, sheet = sheet),
    error = function(e) stop(sprintf("Error reading sheet '%s' in '%s'", sheet, xlsx))
  )
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  names(df) <- stringr::str_trim(names(df)) # Prevent whitespace issues in headers
  df
}

# Normalize IDs: Trim whitespace, empty strings to NA (preserves vector length)
normalize_ids <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x[x == ""] <- NA_character_
  x
}

# Extract unique prot_acc set from a dataframe
get_prot_acc_set <- function(df, prot_col = "prot_acc") {
  names(df) <- stringr::str_trim(names(df))
  if (!prot_col %in% names(df)) stop(sprintf("Column '%s' not found.", prot_col))
  
  x <- normalize_ids(df[[prot_col]])
  x <- x[!is.na(x)]
  unique(x)
}

# Automatically guess the SYMBOL column name if exact match fails
guess_symbol_col <- function(df, prefer = "SYMBOL") {
  nms <- stringr::str_trim(names(df))
  
  # 1. Prefer exact match (case-insensitive)
  if (prefer %in% nms) return(prefer)
  hit <- nms[tolower(nms) == tolower(prefer)]
  if (length(hit) >= 1) return(hit[1])
  
  # 2. Fuzzy match for common alternatives (symbol, gene, name, locus, tair)
  cand <- nms[stringr::str_detect(tolower(nms), "symbol|gene|name|locus|tair")]
  if (length(cand) >= 1) return(cand[1])
  
  NA_character_
}

# Create a mapping table: prot_acc -> SYMBOL (collapses multiple symbols with ';')
make_id_symbol_map <- function(df, prot_col = "prot_acc", sym_col = "SYMBOL") {
  names(df) <- stringr::str_trim(names(df))
  if (!prot_col %in% names(df)) stop(sprintf("Column '%s' not found.", prot_col))
  
  sym_use <- sym_col
  if (!sym_use %in% names(df)) {
    sym_use <- guess_symbol_col(df, prefer = sym_col)
  }
  
  prot <- normalize_ids(df[[prot_col]])
  sym  <- if (!is.na(sym_use)) normalize_ids(df[[sym_use]]) else rep(NA_character_, length(prot))
  
  tmp <- data.frame(
    prot_acc = prot,
    SYMBOL   = sym,
    stringsAsFactors = FALSE
  )
  
  # Remove rows without prot_acc
  tmp <- tmp[!is.na(tmp$prot_acc), , drop = FALSE]
  
  # Aggregate symbols per prot_acc
  tmp %>%
    dplyr::group_by(prot_acc) %>%
    dplyr::summarise(
      SYMBOL = {
        s <- SYMBOL[!is.na(SYMBOL)]
        if (length(s) == 0) NA_character_ else paste(sort(unique(s)), collapse = ";")
      },
      .groups = "drop"
    )
}

# Sanitize Excel sheet names (remove invalid chars + 31 char limit)
sanitize_sheet_name <- function(x) {
  x <- stringr::str_replace_all(x, "[\\[\\]\\:\\*\\?\\/\\\\]", "_")
  if (nchar(x) > 31) x <- substr(x, 1, 31)
  x
}

# Calculate Jaccard Index
calc_jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  uni   <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

# ----------------------------
# 4) Main Processing
# ----------------------------

message("Loading datasets...")

# Load raw dataframes
airid_dfs <- setNames(
  lapply(SHEETS_AIRID, function(sh) read_clean_sheet(IN_XLSX_AIRID, sh)),
  paste0("AirID_", SHEETS_AIRID)
)
coip_dfs <- setNames(
  lapply(SHEETS_COIP, function(sh) read_clean_sheet(IN_XLSX_COIP, sh)),
  paste0("CoIP_", SHEETS_COIP)
)

# Extract sets and create symbol maps
sets_airid <- lapply(airid_dfs, get_prot_acc_set)
sets_coip  <- lapply(coip_dfs,  get_prot_acc_set)

map_airid <- lapply(airid_dfs, make_id_symbol_map)
map_coip  <- lapply(coip_dfs,  make_id_symbol_map)

warm_names <- names(sets_airid)
cool_names <- names(sets_coip)

# Initialize summary matrices
count_mat <- matrix(0L, nrow = length(warm_names), ncol = length(cool_names),
                    dimnames = list(warm_names, cool_names))
jacc_mat  <- matrix(NA_real_, nrow = length(warm_names), ncol = length(cool_names),
                    dimnames = list(warm_names, cool_names))

# Initialize Workbook
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "summary_counts")
openxlsx::addWorksheet(wb, "summary_jaccard")

message("Processing intersections (4x4)...")

# Iterate through all 16 combinations
for (wn in warm_names) {
  for (cn in cool_names) {
    
    # Calculate intersection and stats
    inter_vec <- sort(intersect(sets_airid[[wn]], sets_coip[[cn]]))
    count_mat[wn, cn] <- length(inter_vec)
    jacc_mat[wn, cn]  <- calc_jaccard(sets_airid[[wn]], sets_coip[[cn]])
    
    # Retrieve maps
    ma <- map_airid[[wn]]
    mc <- map_coip[[cn]]
    
    # Safety check for missing SYMBOL columns
    if (!"SYMBOL" %in% names(ma)) ma$SYMBOL <- NA_character_
    if (!"SYMBOL" %in% names(mc)) mc$SYMBOL <- NA_character_
    
    # Rename columns for clarity in merged output
    ma <- dplyr::rename(ma, SYMBOL_AirID = SYMBOL)
    mc <- dplyr::rename(mc, SYMBOL_CoIP  = SYMBOL)
    
    # Create detailed dataframe
    df_out <- data.frame(prot_acc = inter_vec, stringsAsFactors = FALSE) %>%
      dplyr::left_join(ma, by = "prot_acc") %>%
      dplyr::left_join(mc, by = "prot_acc")
    
    # Write to sheet
    sh_name <- sanitize_sheet_name(paste0("pair_", wn, "_x_", cn))
    openxlsx::addWorksheet(wb, sh_name)
    openxlsx::writeData(wb, sh_name, df_out)
    openxlsx::setColWidths(wb, sh_name, cols = 1:ncol(df_out), widths = "auto")
  }
}

# Write Summaries
openxlsx::writeData(wb, "summary_counts",  as.data.frame(count_mat), rowNames = TRUE)
openxlsx::writeData(wb, "summary_jaccard", as.data.frame(round(jacc_mat, 3)), rowNames = TRUE)

# ----------------------------
# 5) Styling & Saving
# ----------------------------
style_header <- openxlsx::createStyle(textDecoration = "bold", halign = "center")
style_int    <- openxlsx::createStyle(numFmt = "0")
style_dec    <- openxlsx::createStyle(numFmt = "0.000")

# Apply style: Counts
openxlsx::addStyle(wb, "summary_counts", style_header,
                   rows = 1, cols = 1:(ncol(count_mat)+1), gridExpand = TRUE, stack = TRUE)
openxlsx::addStyle(wb, "summary_counts", style_int,
                   rows = 2:(nrow(count_mat)+1), cols = 2:(ncol(count_mat)+1), gridExpand = TRUE, stack = TRUE)
openxlsx::setColWidths(wb, "summary_counts", cols = 1:(ncol(count_mat)+1), widths = "auto")

# Apply style: Jaccard
openxlsx::addStyle(wb, "summary_jaccard", style_header,
                   rows = 1, cols = 1:(ncol(jacc_mat)+1), gridExpand = TRUE, stack = TRUE)
openxlsx::addStyle(wb, "summary_jaccard", style_dec,
                   rows = 2:(nrow(jacc_mat)+1), cols = 2:(ncol(jacc_mat)+1), gridExpand = TRUE, stack = TRUE)
openxlsx::setColWidths(wb, "summary_jaccard", cols = 1:(ncol(jacc_mat)+1), widths = "auto")

# Save file
openxlsx::saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
message("Completed. File saved to: ", normalizePath(OUT_XLSX))