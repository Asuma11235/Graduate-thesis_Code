#!/usr/bin/env Rscript

# ============================================================
# Euler Diagram Visualization for 8 Sets (AirID + CoIP)
#
# Description:
#   Reads protein/gene identifiers from 8 Excel sheets.
#   Generates an approximate Euler diagram to visualize overlaps.
#   - AirID sets: Warm colors
#   - CoIP sets:  Cool colors
#
# Input:
#   - AirID Excel file (4 sheets)
#   - CoIP Excel file (4 sheets)
#
# Output:
#   - Euler diagram (PDF & PNG)
# ============================================================

# ----------------------------
# 1) Configuration
# ----------------------------

# Input Files
IN_XLSX_AIRID <- "./input_data/AirID_Dataset.xlsx"
SHEETS_AIRID  <- c("DMSO", "bikinin", "mock", "BAK1")

IN_XLSX_COIP  <- "./input_data/CoIP_Dataset.xlsx"
SHEETS_COIP   <- c("DMSO", "bikinin", "BR-up", "BR-down")

# Target Column (e.g., "SYMBOL" or "prot_acc")
ID_COLUMN     <- "SYMBOL"

# Output Files
OUT_DIR       <- "./results_euler"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_PDF       <- file.path(OUT_DIR, "AirID_CoIP_Euler8sets_Symbol.pdf")
OUT_PNG       <- file.path(OUT_DIR, "AirID_CoIP_Euler8sets_Symbol.png")

# ----------------------------
# 2) Libraries
# ----------------------------
pkgs <- c("readxl", "dplyr", "stringr", "eulerr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# ----------------------------
# 3) Helper Functions
# ----------------------------

# Read unique identifiers from a specific sheet
read_id_set <- function(xlsx, sheet, col = ID_COLUMN) {
  df <- tryCatch(
    readxl::read_excel(xlsx, sheet = sheet),
    error = function(e) stop(sprintf("Failed to read sheet '%s' from '%s'", sheet, xlsx))
  )
  
  if (!col %in% colnames(df)) {
    stop(sprintf("Column '%s' not found in sheet '%s' of '%s'", col, sheet, xlsx))
  }
  
  x <- df[[col]] |>
    as.character() |>
    stringr::str_trim() |>
    stats::na.omit()
  
  # Remove empty strings
  x <- x[nzchar(x)]
  
  # Optional: Remove additional annotations if necessary (e.g., "ID description")
  # x <- sub("\\s+.*$", "", x)
  
  unique(x)
}

# ----------------------------
# 4) Main Processing
# ----------------------------

message("Loading datasets...")

# Load AirID sets
sets_airid <- setNames(
  lapply(SHEETS_AIRID, function(sh) read_id_set(IN_XLSX_AIRID, sh)),
  paste0("AirID_", SHEETS_AIRID)
)

# Load CoIP sets
sets_coip <- setNames(
  lapply(SHEETS_COIP, function(sh) read_id_set(IN_XLSX_COIP, sh)),
  paste0("CoIP_", SHEETS_COIP)
)

# Combine all sets
sets_all <- c(sets_airid, sets_coip)

# --- Sanity Check ---
message("=== Set Sizes ===")
print(sapply(sets_all, length))

# --- Fit Euler Diagram ---
# shape="ellipse" is generally more stable for many sets
message("Fitting Euler diagram...")
fit <- eulerr::euler(sets_all, shape = "ellipse")

# --- Define Colors ---
# AirID (Warm colors) vs CoIP (Cool colors)
fill_cols <- c(
  # AirID
  "AirID_DMSO"    = "#F4A261",
  "AirID_bikinin" = "#E76F51",
  "AirID_mock"    = "#F6BD60",
  "AirID_BAK1"    = "#E9C46A",
  # CoIP
  "CoIP_DMSO"     = "#4D96FF",
  "CoIP_bikinin"  = "#1F77B4",
  "CoIP_BR-up"    = "#3FB6B2",
  "CoIP_BR-down"  = "#2A9D8F"
)

# Verify color mapping
missing_cols <- setdiff(names(sets_all), names(fill_cols))
if (length(missing_cols) > 0) {
  stop("Missing color definitions for: ", paste(missing_cols, collapse = ", "))
}

# --- Plotting Function ---
plot_euler <- function() {
  print(plot(
    fit,
    fills = list(fill = fill_cols, alpha = 0.45),
    edges = list(col = "grey30", lwd = 1),
    labels = list(col = "grey10", cex = 0.9),
    quantities = list(col = "grey10", cex = 0.75),
    legend = list(cex = 0.8),
    main = "AirID (Warm) vs CoIP (Cool): Overlap Analysis"
  ))
}

# --- Save Outputs ---
pdf(OUT_PDF, width = 10, height = 7)
plot_euler()
invisible(dev.off())

png(OUT_PNG, width = 1800, height = 1200, res = 200)
plot_euler()
invisible(dev.off())

message("\nSaved Euler diagrams:")
message("- ", normalizePath(OUT_PDF))
message("- ", normalizePath(OUT_PNG))