#!/usr/bin/env Rscript

# ============================================================
# Venn/Euler Region Exporter
#
# Description:
#   Exports data corresponding to specific Venn/Euler diagram regions (intersections)
#   into a multi-sheet Excel file.
#   - Identifies which conditions (sheets) each protein belongs to.
#   - Groups proteins by their specific intersection pattern (e.g., "bikinin & DMSO").
#   - Creates a separate Excel sheet for each region.
#   - Generates a summary sheet with counts per region.
#
# Input:
#   - Excel file containing multiple sheets (one per condition).
#
# Output:
#   - Excel file where each sheet represents a specific Venn region.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
  library(tidyr)
  library(readxl)
  library(stringr)
})

# ==============================================================================
# Configuration
# ==============================================================================
# Input/Output paths
IN_XLSX          <- "path/to/your/CoIP_data.xlsx"
OUT_XLSX_REGIONS <- "path/to/output/Venn_BR_regions_details.xlsx"

# Data Structure
SHEETS   <- c("bikinin", "DMSO", "BR-up", "BR-down")
PROT_COL <- "prot_acc"  # Column containing unique Protein IDs

# ==============================================================================
# Main Process
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Read and Merge Data
# ------------------------------------------------------------------------------
# Read each sheet, clean column names, and append a source identifier
df_list <- lapply(SHEETS, function(sh) {
  # Check if sheet exists
  if (!sh %in% excel_sheets(IN_XLSX)) {
    warning(sprintf("Sheet '%s' not found. Skipping.", sh))
    return(NULL)
  }
  
  df <- read_excel(IN_XLSX, sheet = sh)
  names(df) <- str_trim(names(df)) # Remove whitespace from column names
  
  # Validate and clean Protein ID column
  if (PROT_COL %in% names(df)) {
    df[[PROT_COL]] <- str_trim(as.character(df[[PROT_COL]]))
    df <- df[!is.na(df[[PROT_COL]]) & nzchar(df[[PROT_COL]]), ]
  } else {
    warning(paste("Sheet", sh, "does not have column", PROT_COL))
    return(NULL)
  }
  
  # Add source column to track origin
  df$`_Source_` <- sh
  return(df)
})

# Filter out NULLs (failed reads) and bind into a single master dataframe
df_list <- df_list[!sapply(df_list, is.null)]
df_master <- bind_rows(df_list)

# ------------------------------------------------------------------------------
# 2. Determine Membership (Intersections)
# ------------------------------------------------------------------------------
# Group by Protein ID and determine the exact combination of conditions
# e.g., ProtA -> "bikinin & DMSO"
membership_df <- df_master %>%
  group_by(!!sym(PROT_COL)) %>%
  summarise(
    # Sort source names alphabetically to ensure consistency
    Belongs_To = paste(sort(unique(`_Source_`)), collapse = " & "),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 3. Export Each Region to Excel
# ------------------------------------------------------------------------------
# Identify all unique intersection patterns
regions <- sort(unique(membership_df$Belongs_To))

message("=== Detected Regions ===")
print(regions)

# Initialize new Workbook
wb <- createWorkbook()

for (region_name in regions) {
  
  # Identify proteins belonging to this specific region
  target_prots <- membership_df %>%
    filter(Belongs_To == region_name) %>%
    pull(!!sym(PROT_COL))
  
  # Extract original rows for these proteins
  # (Includes rows from all sources where the protein was present)
  region_data <- df_master %>%
    filter(!!sym(PROT_COL) %in% target_prots) %>%
    arrange(!!sym(PROT_COL), `_Source_`) # Sort by Protein ID, then Source
  
  # Create a valid Excel sheet name (Max 31 chars, no special symbols)
  sheet_name <- region_name %>%
    gsub(" & ", "_", .) %>%  # Replace separator
    gsub("-", "", .)         # Remove hyphens (optional style choice)
  
  # Truncate to 31 characters if necessary
  if (nchar(sheet_name) > 31) {
    sheet_name <- substr(sheet_name, 1, 31)
  }
  
  # Add sheet and write data
  if (nrow(region_data) > 0) {
    # Check if sheet name already exists (handling truncation collisions)
    if (sheet_name %in% names(wb)) {
      warning(sprintf("Sheet name collision for '%s'. Appending suffix.", sheet_name))
      sheet_name <- substr(paste0(sheet_name, "_2"), 1, 31)
    }
    
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, region_data)
    message(sprintf("Sheet added: %-25s (%d proteins, %d rows)", 
                    sheet_name, length(target_prots), nrow(region_data)))
  }
}

# ------------------------------------------------------------------------------
# 4. Create Summary Sheet
# ------------------------------------------------------------------------------
summary_df <- membership_df %>%
  group_by(Belongs_To) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))

addWorksheet(wb, "Summary_Counts")
writeData(wb, "Summary_Counts", summary_df)

# ------------------------------------------------------------------------------
# Save Output
# ------------------------------------------------------------------------------
# Ensure output directory exists
out_dir <- dirname(OUT_XLSX_REGIONS)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

saveWorkbook(wb, OUT_XLSX_REGIONS, overwrite = TRUE)

message("\nSuccessfully saved Excel file to:")
message(normalizePath(OUT_XLSX_REGIONS))