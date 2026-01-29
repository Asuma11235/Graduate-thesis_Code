#!/usr/bin/env Rscript

# ============================================================
# Fisher's Exact Test for Set Overlap
#
# Description:
#   Evaluates the statistical significance of the overlap between two sets (A and B)
#   within a defined universe (Total population).
#   - Constructs a 2x2 contingency table.
#   - Performs a two-sided Fisher's Exact Test.
#   - Reports Odds Ratio, P-value, and biological interpretation.
#

# ============================================================

# ==============================================================================
# Function Definition
# ==============================================================================
calculate_fisher_overlap <- function(n_a, n_b, n_overlap, n_total) {
  
  # 1. Calculate 2x2 Contingency Table Elements
  # a: In both A and B (Overlap)
  # b: In A, but not B
  # c: In B, but not A
  # d: In neither (Background)
  a <- n_overlap
  b <- n_a - n_overlap
  c <- n_b - n_overlap
  d <- n_total - n_a - n_b + n_overlap
  
  # Error Check: Ensure universe size is sufficient
  if (d < 0) {
    message("Error: n_total is smaller than the union of Set A and Set B.")
    return(invisible(NULL))
  }
  
  # 2. Construct Matrix
  # Rows: Set A status (+/-)
  # Cols: Set B status (+/-)
  tab <- matrix(
    c(a, b, c, d),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("In Set A", "Not in Set A"),
      c("In Set B", "Not in Set B")
    )
  )
  
  # 3. Perform Fisher's Exact Test (Two-sided)
  ft <- fisher.test(tab, alternative = "two.sided")
  
  # Extract statistics
  odds_ratio <- unname(ft$estimate)
  p_value    <- ft$p.value
  
  # 4. Print Results
  cat("\n--- Contingency Table ---\n")
  print(as.data.frame.matrix(tab))
  cat(paste0(strrep("-", 35), "\n"))
  cat(sprintf("Odds Ratio : %.4f\n", odds_ratio))
  cat(sprintf("P-value    : %.4e\n", p_value))
  cat(paste0(strrep("-", 35), "\n"))
  
  # 5. Interpretation
  if (p_value < 0.05) {
    cat("Result: Significant association (p < 0.05)\n")
    if (odds_ratio > 1) {
      cat("Trend : Positive association (Overlap is greater than expected)\n")
    } else {
      cat("Trend : Negative association (Overlap is less than expected)\n")
    }
  } else {
    cat("Result: No significant association (p >= 0.05)\n")
  }
  cat(paste0(strrep("=", 35), "\n\n"))
  
  # Return results invisibly for further use
  invisible(list(table = tab, fisher = ft))
}

# ==============================================================================
# Main Execution Example
# ==============================================================================

# Input Parameters
# Change these values according to your data
num_A       <- 22   # Size of Set A
num_B       <- 10   # Size of Set B
num_overlap <- 6    # Overlap count
num_total   <- 216  # Total Universe Size (Background population)

# Execute Function
calculate_fisher_overlap(num_A, num_B, num_overlap, num_total)