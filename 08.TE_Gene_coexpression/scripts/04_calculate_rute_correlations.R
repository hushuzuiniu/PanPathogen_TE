#!/usr/bin/env Rscript

# Step 4: calculate expression correlations between ruTEs and nearby genes.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
project_root <- if (length(script_arg) > 0L) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."))
} else {
  normalizePath(getwd())
}
source(file.path(project_root, "R", "project_setup.R"))

assert_packages(c("data.table", "doParallel", "dplyr", "foreach"))

library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

assert_files(c(TE_GENE_MAP_FILE, RUTE_LIST_FILE))
ensure_parent_dirs(c(CORRELATION_RDS, CORRELATION_TSV, CORRELATION_SUMMARY_TSV))

te_gene_map_file <- TE_GENE_MAP_FILE
ruTE_list_file <- RUTE_LIST_FILE
expr_dir <- NORMALIZED_EXPRESSION_DIR
output_file1 <- CORRELATION_RDS
output_file2 <- CORRELATION_TSV
output_file3 <- CORRELATION_SUMMARY_TSV

min_timepoints <- MIN_TIMEPOINTS
n_cores <- N_CORES

datasets <- DATASETS

# Average biological replicates within each time point.
aggregate_replicates <- function(expr_matrix, timepoint_pattern = "_rep") {
  ID_col <- colnames(expr_matrix)[1]
  col_names <- colnames(expr_matrix)[-1]
  
  timepoints <- gsub(paste0(timepoint_pattern, ".*"), "", col_names)
  unique_tps <- unique(timepoints)
  unique_tps <- unique_tps[order(as.numeric(gsub("t", "", unique_tps)))]
  
  result <- data.table(ID = expr_matrix[[ID_col]])
  
  for(tp in unique_tps) {
    cols <- which(timepoints == tp)
    if(length(cols) == 1) {
      result[[tp]] <- expr_matrix[[cols + 1]]
    } else {
      result[[tp]] <- apply(expr_matrix[, cols + 1, with = FALSE], 1, mean, na.rm = TRUE)
    }
  }
  
  return(result)
}

cat("Phase 1: Reading the ruTE list and TE-gene map...\n")

ruTE_ids <- fread(ruTE_list_file, header = FALSE)[[1]]
ruTE_ids <- gsub('"', '', ruTE_ids)
ruTE_ids <- trimws(ruTE_ids)
cat(sprintf("Number of ruTEs: %d\n", length(ruTE_ids)))

if(!file.exists(te_gene_map_file)) {
  stop(paste("TE-gene map file not found:", te_gene_map_file))
}
te_gene_map <- fread(te_gene_map_file, header = TRUE)
cat(sprintf("Loaded TE-gene map: %d rows\n", nrow(te_gene_map)))
cat("Columns:", paste(colnames(te_gene_map), collapse = ", "), "\n")

ruTE_gene_map <- te_gene_map[te_id %in% ruTE_ids & !is.na(closest_gene_in_50kb), 
                             .(te_id, gene_id = closest_gene_in_50kb)]
cat(sprintf("ruTEs with nearest gene: %d\n", nrow(ruTE_gene_map)))

if(nrow(ruTE_gene_map) == 0) {
  stop("No ruTEs with nearest gene found")
}

cat("\nPhase 2: Creating all ruTE-by-dataset tasks...\n")

task_list <- list()
task_id <- 1

for(i in 1:nrow(ruTE_gene_map)) {
  te <- ruTE_gene_map$te_id[i]
  gene <- ruTE_gene_map$gene_id[i]
  
  for(ds in datasets) {
    task_list[[task_id]] <- data.table(
      te_id = te,
      gene_id = gene,
      dataset = ds,
      stringsAsFactors = FALSE
    )
    task_id <- task_id + 1
  }
}

tasks_df <- rbindlist(task_list)
cat(sprintf("Total tasks: %d (%d ruTEs × %d datasets)\n", 
            nrow(tasks_df), nrow(ruTE_gene_map), length(datasets)))

cat("\nPhase 3: Loading and aggregating expression matrices...\n")

expr_list <- list()
for(ds in datasets) {
  te_file <- file.path(expr_dir, paste0(ds, "_TE_log2TPM.txt"))
  gene_file <- file.path(expr_dir, paste0(ds, "_Gene_log2TPM.txt"))
  
  if(!file.exists(te_file)) {
    stop(paste("TE file not found:", te_file))
  }
  if(!file.exists(gene_file)) {
    stop(paste("Gene file not found:", gene_file))
  }
  
  te_raw <- fread(te_file, header = TRUE)
  gene_raw <- fread(gene_file, header = TRUE)
  
  cat(sprintf("  %s: raw TE: %d rows, %d columns; Gene: %d rows, %d columns\n", 
              ds, nrow(te_raw), ncol(te_raw), nrow(gene_raw), ncol(gene_raw)))
  
  te_expr <- aggregate_replicates(te_raw)
  gene_expr <- aggregate_replicates(gene_raw)
  
  setkey(te_expr, ID)
  setkey(gene_expr, ID)
  
  n_timepoints <- ncol(te_expr) - 1
  
  expr_list[[ds]] <- list(
    te = te_expr,
    gene = gene_expr,
    timepoints = n_timepoints,
    timepoint_names = colnames(te_expr)[-1]
  )
  
  cat(sprintf("    Aggregated: %d TEs, %d genes, %d timepoints\n", 
              nrow(te_expr), nrow(gene_expr), n_timepoints))
}

calc_correlation_fast <- function(te_id, gene_id, dataset, expr_list, min_timepoints = 3) {
  
  ds_data <- expr_list[[dataset]]
  
  if(ds_data$timepoints < min_timepoints) {
    return(list(rho = NA, reason = paste0("insufficient_timepoints_", ds_data$timepoints)))
  }
  
  te_row <- ds_data$te[.(te_id)]
  gene_row <- ds_data$gene[.(gene_id)]
  
  if(nrow(te_row) == 0) {
    return(list(rho = NA, reason = "TE_not_found"))
  }
  if(nrow(gene_row) == 0) {
    return(list(rho = NA, reason = "Gene_not_found"))
  }
  
  te_vec <- as.numeric(te_row[1, -1])
  gene_vec <- as.numeric(gene_row[1, -1])
  
  if(any(is.na(te_vec)) || any(is.na(gene_vec))) {
    return(list(rho = NA, reason = "contains_NA"))
  }
  
  if(sd(te_vec) == 0 || sd(gene_vec) == 0) {
    return(list(rho = NA, reason = "zero_variance"))
  }
  
  rho <- cor(te_vec, gene_vec, method = "spearman")
  
  return(list(rho = rho, reason = NA))
}

cat("\nPhase 4: Calculating correlations in parallel...\n")

if (nrow(tasks_df) > 0L) {
  n_cores <- min(n_cores, nrow(tasks_df))
}
cat(sprintf("Using %d cores\n", n_cores))
cl <- makeCluster(n_cores)
registerDoParallel(cl)

clusterExport(cl, c("calc_correlation_fast", "expr_list", "min_timepoints"))

cat("Starting parallel computation...\n")
start_time <- Sys.time()

results_list <- foreach(i = 1:nrow(tasks_df), 
                        .combine = rbind,
                        .packages = c("data.table")) %dopar% {
                          
                          task <- tasks_df[i]
                          result <- calc_correlation_fast(
                            te_id = task$te_id,
                            gene_id = task$gene_id,
                            dataset = task$dataset,
                            expr_list = expr_list,
                            min_timepoints = min_timepoints
                          )
                          
                          data.table(
                            te_id = task$te_id,
                            gene_id = task$gene_id,
                            dataset = task$dataset,
                            rho = result$rho,
                            reason = ifelse(is.na(result$reason), NA_character_, result$reason)
                          )
                        }

stopCluster(cl)

end_time <- Sys.time()
cat(sprintf("Computation completed in %.2f minutes\n", 
            difftime(end_time, start_time, units = "mins")))

cat("\nPhase 5: Adding pathogen metadata...\n")

dataset_to_pathogen <- data.table(
  dataset = datasets,
  pathogen = c("Virus_EBV", "Virus_EBV", 
               "Fungi_Af", "Fungi_Af",
               "Virus_IAV", "Virus_RSV", "Virus_HCV")
)

results_list <- results_list %>%
  left_join(dataset_to_pathogen, by = "dataset") %>%
  select(te_id, gene_id, pathogen, dataset, rho, reason)
results_list <- as.data.table(results_list)

cat("\nPhase 6: Summarizing correlation results...\n")

cat(sprintf("Total results: %d\n", nrow(results_list)))

successful <- results_list[!is.na(rho)]
failed <- results_list[is.na(rho)]

cat(sprintf("Successful calculations: %d (%.1f%%)\n", 
            nrow(successful), 100 * nrow(successful) / nrow(results_list)))
cat(sprintf("Failed calculations: %d (%.1f%%)\n", 
            nrow(failed), 100 * nrow(failed) / nrow(results_list)))

if(nrow(failed) > 0) {
  cat("\nFailure reasons:\n")
  reason_counts <- failed[, .N, by = reason][order(-N)]
  print(reason_counts)
}

if(nrow(successful) > 0) {
  cat("\nCorrelation distribution:\n")
  cat(sprintf("  Min: %.3f\n", min(successful$rho, na.rm = TRUE)))
  cat(sprintf("  Q1: %.3f\n", quantile(successful$rho, 0.25, na.rm = TRUE)))
  cat(sprintf("  Median: %.3f\n", median(successful$rho, na.rm = TRUE)))
  cat(sprintf("  Mean: %.3f\n", mean(successful$rho, na.rm = TRUE)))
  cat(sprintf("  Q3: %.3f\n", quantile(successful$rho, 0.75, na.rm = TRUE)))
  cat(sprintf("  Max: %.3f\n", max(successful$rho, na.rm = TRUE)))
}

cat("\nPhase 7: Saving results...\n")

saveRDS(results_list, output_file1)
cat(sprintf("Saved RDS to: %s\n", output_file1))

fwrite(results_list, output_file2, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("Saved TSV to: %s\n", output_file2))

if(nrow(successful) > 0) {
  # Summarize each TE by its median correlation across valid datasets.
  te_summary <- successful[, .(median_rho = median(rho, na.rm = TRUE)), by = te_id]
  
  summary_stats <- data.table(
    metric = c("total_tasks", "successful", "failed", 
               "median_rho", "mean_rho", "min_rho", "max_rho",
               "n_ruTEs", "ruTE_median_rho"),
    value = c(nrow(results_list), nrow(successful), nrow(failed),
              median(successful$rho, na.rm = TRUE),
              mean(successful$rho, na.rm = TRUE),
              min(successful$rho, na.rm = TRUE),
              max(successful$rho, na.rm = TRUE),
              nrow(te_summary),
              median(te_summary$median_rho, na.rm = TRUE))
  )
  fwrite(summary_stats, output_file3, sep = "\t", quote = FALSE)
  cat(sprintf("Saved summary to: %s\n", output_file3))
}

cat("\nTime points per dataset after aggregation:\n")
for(ds in datasets) {
  cat(sprintf("  %s: %d timepoints (%s)\n", 
              ds, expr_list[[ds]]$timepoints,
              paste(expr_list[[ds]]$timepoint_names, collapse = ", ")))
}

cat("\nStep 4 completed.\n")
