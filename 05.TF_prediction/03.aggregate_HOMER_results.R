#!/usr/bin/env Rscript
# Aggregate HOMER knownResults.txt from all subfamilies into a single table

library(tidyverse)

# Set working directory
setwd("/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_prediction_add_bg_v3")

# Get all subfamily directories
result_dirs <- list.dirs("HOMER_results_per_subfamily", recursive = FALSE)
subfamily_names <- basename(result_dirs)

# Read and combine all knownResults.txt files with standardized column names
all_results <- list()

for (i in seq_along(result_dirs)) {
  dir <- result_dirs[i]
  sf <- subfamily_names[i]
  result_file <- file.path(dir, "knownResults.txt")
  
  if (file.exists(result_file)) {
    df <- read.delim(result_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    
    # Standardize column names
    colnames(df) <- c("Motif_Name", "Consensus", "P_value", "Log_P_value", "q_value_Benjamini",
                      "Target_with_Motif_count", "Target_with_Motif_pct",
                      "Background_with_Motif_count", "Background_with_Motif_pct")
    
    # Convert percentage columns to numeric (remove % sign)
    df$Target_with_Motif_pct <- as.numeric(gsub("%", "", df$Target_with_Motif_pct))
    df$Background_with_Motif_pct <- as.numeric(gsub("%", "", df$Background_with_Motif_pct))
    
    df$Subfamily <- sf
    all_results[[sf]] <- df
    cat("Read:", sf, "- ", nrow(df), "motifs\n")
  } else {
    cat("Missing:", sf, "\n")
  }
}

# Combine all results
combined_df <- bind_rows(all_results)

# Extract TF name from Motif Name (remove the part in parentheses for cleaner names)
combined_df$TF_Name <- gsub("\\(.*\\)", "", combined_df$Motif_Name)
combined_df$TF_Name <- trimws(combined_df$TF_Name)

# Calculate fold enrichment
# Fold = (% Target) / (% Background)
combined_df$Fold_Enrichment <- combined_df$Target_with_Motif_pct / combined_df$Background_with_Motif_pct
combined_df$Fold_Enrichment[is.infinite(combined_df$Fold_Enrichment)] <- NA
combined_df$Fold_Enrichment[is.nan(combined_df$Fold_Enrichment)] <- NA

# Save combined raw results
write.table(combined_df, "combined_HOMER_results_all_subfamilies.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nSaved combined results to: combined_HOMER_results_all_subfamilies.txt\n")
cat("Total rows:", nrow(combined_df), "\n")

# Filter significant results (P-value < 1e-5)
sig_df <- combined_df %>%
  filter(P_value < 1e-5)
# Filter 3-fold enrichment 
sig_df <- sig_df %>%
  filter(Fold_Enrichment >= 3)

cat("Significant results (P < 0.05):", nrow(sig_df), "\n")

write.table(sig_df, "combined_HOMER_results_significant_p1e-5_fc3.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Filter significant results (P-value < 1e-5)
sig_df <- combined_df %>%
  filter(P_value < 1e-5)
# Filter 3-fold enrichment 
sig_df <- sig_df %>%
  filter(Fold_Enrichment >= 1)

cat("Significant results (P < 0.05):", nrow(sig_df), "\n")

write.table(sig_df, "combined_HOMER_results_significant_p1e-5_fc1.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
