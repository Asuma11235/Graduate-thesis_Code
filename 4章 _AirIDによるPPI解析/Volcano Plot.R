#!/usr/bin/env Rscript

# ============================================================
# Volcano Plot & Differential Abundance Analysis
#
# Description:
#   Generates volcano plots for all (Abundance Ratio) x (P-Value) pairs found in the input.
#   - Automatically pairs columns by matching suffixes.
#   - Classifies proteins as Up/Down/Not Significant based on fold-change and p-value.
#   - Exports:
#       1. PNG plots (saved to a directory).
#       2. PPTX presentation (one slide per comparison).
#       3. Excel file with UP/DOWN protein lists and summary statistics.
#
# Input:
#   - Excel file containing "Abundance Ratio" and "P-Value" columns.
#
# Output:
#   - PPTX, XLSX, and PNG files.
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(stringr)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
  library(officer)
})

# ==============================================================================
# Configuration
# ==============================================================================
# Input/Output paths
IN_XLSX     <- "path/to/your/AmPR2_data.xlsx"
OUT_PPTX    <- "path/to/output/Volcano_AllComparisons.pptx"
OUT_XLSX    <- "path/to/output/Volcano_UP_DOWN_lists.xlsx"
OUT_PNG_DIR <- "path/to/output/volcano_pngs" # Directory for PNGs

# Excel Settings
SHEET_IDX   <- 1  # Sheet index or name (e.g., "Sheet1")

# Thresholds for Significance
ALPHA       <- 0.05
FC_UP       <- 2   # Up-regulated if Ratio >= FC_UP (and p < ALPHA)
# Down-regulated if Ratio <= 1/FC_UP (and p < ALPHA)

# Labeling Settings
LABEL_TOP_N <- 15  # Label top N proteins (smallest p-value) for Up/Down

# ==============================================================================
# Helper Functions
# ==============================================================================

# Numeric coercion helper (handles commas and empty strings)
to_num <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x[x == ""] <- NA_character_
  x <- gsub(",", "", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

# Create safe Excel sheet name (Max 31 chars, no illegal chars)
make_sheet_name <- function(prefix, suf, idx) {
  s <- gsub("[:\\\\/\\?\\*\\[\\]]", "_", suf) # Remove illegal chars
  s <- paste0(prefix, "_", idx, "_", s)
  s <- substr(s, 1, 31)
  s
}

# ==============================================================================
# 1. Load Data & Identify Columns
# ==============================================================================
stopifnot(file.exists(IN_XLSX))
if (!dir.exists(OUT_PNG_DIR)) dir.create(OUT_PNG_DIR, recursive = TRUE)

df <- read_excel(IN_XLSX, sheet = SHEET_IDX) %>%
  as.data.frame(stringsAsFactors = FALSE)
names(df) <- str_trim(names(df))

# Annotation Columns
anno_cols_wanted <- c(
  "Accession", "Description", "Biological Process", "Cellular Component",
  "Molecular Function", "Gene Symbol", "Ensembl Gene ID"
)
anno_cols_exist <- intersect(anno_cols_wanted, names(df))
anno_cols_missing <- setdiff(anno_cols_wanted, anno_cols_exist)

if (length(anno_cols_missing) > 0) {
  message("NOTE: The following annotation columns are missing and will be skipped:")
  message(paste(anno_cols_missing, collapse = ", "))
}

# Identify Abundance Columns
scaled_cols <- grep("^Abundances \\(Scaled\\):", names(df), value = TRUE)

# Identify Ratio and P-value Columns
ratio_prefix <- "^Abundance Ratio:\\s*"
pval_prefix  <- "^Abundance Ratio P-Value:\\s*"

ratio_cols <- grep(ratio_prefix, names(df), value = TRUE)
pval_cols  <- grep(pval_prefix,  names(df), value = TRUE)

if (length(ratio_cols) == 0 || length(pval_cols) == 0) {
  stop("No ratio or p-value columns found.\n",
       "Ratio cols found: ", length(ratio_cols), "\n",
       "P-value cols found: ", length(pval_cols))
}

# Match Pairs by Suffix
ratio_suffix <- sub(ratio_prefix, "", ratio_cols)
pval_suffix  <- sub(pval_prefix,  "", pval_cols)

common_suffix <- intersect(ratio_suffix, pval_suffix)
if (length(common_suffix) == 0) {
  stop("No matching suffixes found between ratio and p-value columns.\n",
       "Example Ratio Suffix: ", head(ratio_suffix, 3), "\n",
       "Example Pval Suffix: ", head(pval_suffix, 3))
}

pairs <- lapply(common_suffix, function(suf) {
  list(
    suffix    = suf,
    ratio_col = ratio_cols[match(suf, ratio_suffix)],
    pval_col  = pval_cols[match(suf, pval_suffix)]
  )
})

message("=== Matched Comparisons ===")
for (p in pairs) message(" - ", p$suffix)

# Determine Label Column
label_col <- if ("Gene Symbol" %in% names(df)) "Gene Symbol" else if ("Accession" %in% names(df)) "Accession" else NULL

# Color Definition
COLS_VOLCANO <- c(
  Up   = "#D32F2F",  # Red
  Down = "#1976D2",  # Blue
  NS   = "#9E9E9E"   # Grey
)

# ==============================================================================
# 2. Process Each Comparison (Loop)
# ==============================================================================
wb <- createWorkbook()
addWorksheet(wb, "SUMMARY")
summary_rows <- list()

ppt <- read_pptx() # Initialize empty slide deck

for (i in seq_along(pairs)) {
  suf       <- pairs[[i]]$suffix
  ratio_col <- pairs[[i]]$ratio_col
  pval_col  <- pairs[[i]]$pval_col
  
  # Prepare Data
  dat <- df %>%
    mutate(
      .ratio = to_num(.data[[ratio_col]]),
      .pval  = to_num(.data[[pval_col]])
    )
  
  # Sanitize p-values (Replace 0 with min positive double to allow log transform)
  dat$.pval[is.finite(dat$.pval) & dat$.pval <= 0] <- .Machine$double.xmin
  
  # Calculate Plot Axes
  dat <- dat %>%
    mutate(
      log2_ratio = ifelse(is.finite(.ratio) & .ratio > 0, log2(.ratio), NA_real_),
      neglog10_p = ifelse(is.finite(.pval)  & .pval  > 0, -log10(.pval), NA_real_)
    )
  
  # Classify Proteins
  dat <- dat %>%
    mutate(
      category = case_when(
        !is.na(.pval) & .pval < ALPHA & !is.na(.ratio) & .ratio >= FC_UP       ~ "Up",
        !is.na(.pval) & .pval < ALPHA & !is.na(.ratio) & .ratio <= 1/FC_UP     ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  n_up <- sum(dat$category == "Up", na.rm = TRUE)
  n_dn <- sum(dat$category == "Down", na.rm = TRUE)
  
  # Select Labels for Top Hits
  label_df <- bind_rows(
    dat %>% filter(category == "Up")   %>% arrange(.pval) %>% slice_head(n = LABEL_TOP_N),
    dat %>% filter(category == "Down") %>% arrange(.pval) %>% slice_head(n = LABEL_TOP_N)
  )
  
  if (!is.null(label_col)) {
    label_df <- label_df %>%
      mutate(.label = as.character(.data[[label_col]])) %>%
      mutate(.label = ifelse(is.na(.label) | .label == "", NA_character_, .label))
  } else {
    label_df$.label <- NA_character_
  }
  
  # Generate Volcano Plot
  vline <- log2(FC_UP)
  hline <- -log10(ALPHA)
  
  p <- ggplot(dat, aes(x = log2_ratio, y = neglog10_p)) +
    geom_point(aes(color = category), alpha = 0.8, size = 1.6, na.rm = TRUE) +
    scale_color_manual(values = COLS_VOLCANO) +
    geom_vline(xintercept = c(-vline, vline), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = hline, linetype = "dashed", color = "grey40") +
    labs(
      title = paste0("Volcano: ", suf),
      subtitle = sprintf("Up=%d, Down=%d | Threshold: p<%.3g, Ratio>=%.2f (Up), <=%.3f (Down)",
                         n_up, n_dn, ALPHA, FC_UP, 1/FC_UP),
      x = "log2(Abundance Ratio)",
      y = "-log10(P-value)",
      color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), legend.position = "top")
  
  if (!is.null(label_col) && any(!is.na(label_df$.label))) {
    p <- p + geom_text(
      data = label_df,
      aes(label = .label),
      size = 3,
      vjust = -0.6,
      check_overlap = TRUE,
      na.rm = TRUE
    )
  }
  
  # Save PNG
  png_path <- file.path(OUT_PNG_DIR, paste0("volcano_", i, ".png"))
  ggsave(png_path, p, width = 10, height = 7, dpi = 200)
  
  # Add Slide to PPTX
  ppt <- ppt %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = paste0("Volcano: ", suf), location = ph_location_type(type = "title")) %>%
    ph_with(value = external_img(png_path, width = 9.5, height = 5.5),
            location = ph_location(left = 0.5, top = 1.3, width = 9.5, height = 5.5))
  
  # Prepare Export Tables
  base_cols <- c(anno_cols_exist, scaled_cols)
  base_cols <- base_cols[base_cols %in% names(dat)] # Ensure columns exist
  
  out_common <- dat %>%
    mutate(
      Comparison = suf,
      Abundance_Ratio = .ratio,
      P_Value = .pval
    )
  
  up_tbl <- out_common %>%
    filter(category == "Up") %>%
    select(Comparison, Abundance_Ratio, P_Value, all_of(base_cols))
  
  dn_tbl <- out_common %>%
    filter(category == "Down") %>%
    select(Comparison, Abundance_Ratio, P_Value, all_of(base_cols))
  
  # Add to Excel Workbook
  sh_up <- make_sheet_name("UP", suf, i)
  sh_dn <- make_sheet_name("DOWN", suf, i)
  
  addWorksheet(wb, sh_up)
  addWorksheet(wb, sh_dn)
  
  writeData(wb, sh_up, up_tbl)
  writeData(wb, sh_dn, dn_tbl)
  
  # Record Summary Info
  summary_rows[[i]] <- data.frame(
    idx = i,
    suffix = suf,
    ratio_col = ratio_col,
    pval_col = pval_col,
    n_up = n_up,
    n_down = n_dn,
    sheet_up = sh_up,
    sheet_down = sh_dn,
    plot_png = png_path,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# 3. Finalize Outputs
# ==============================================================================
# Write Summary Sheet
summary_df <- bind_rows(summary_rows)
writeData(wb, "SUMMARY", summary_df)

# Save Files
saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
print(ppt, target = OUT_PPTX)

message("\nAnalysis Complete.")
message("Saved PPTX: ", normalizePath(OUT_PPTX))
message("Saved XLSX: ", normalizePath(OUT_XLSX))
message("Saved PNGs: ", normalizePath(OUT_PNG_DIR))