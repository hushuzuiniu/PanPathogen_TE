#!/usr/bin/env Rscript

# Step 1: standardize sample names and convert raw counts to log2(TPM + 1).

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
project_root <- if (length(script_arg) > 0L) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."))
} else {
  normalizePath(getwd())
}
source(file.path(project_root, "R", "project_setup.R"))

assert_packages(c("data.table", "dplyr", "stringr"))
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

assert_files(c(SAMPLE_METADATA_FILE, GENE_BED_FILE))
ensure_parent_dirs(FINAL_METADATA_FILE)
dir.create(NORMALIZED_EXPRESSION_DIR, recursive = TRUE, showWarnings = FALSE)

message("Reading sample metadata: ", SAMPLE_METADATA_FILE)
metadata <- fread(SAMPLE_METADATA_FILE, data.table = FALSE, encoding = "UTF-8")
required_metadata_columns <- c("Dataset", "Sample_Name", "Time_Point")
missing_metadata_columns <- setdiff(required_metadata_columns, names(metadata))
if (length(missing_metadata_columns) > 0L) {
  stop(
    "Metadata is missing columns: ",
    paste(missing_metadata_columns, collapse = ", "),
    call. = FALSE
  )
}

metadata <- metadata %>%
  mutate(Time_num = as.numeric(gsub("[^0-9]", "", Time_Point))) %>%
  group_by(Dataset) %>%
  mutate(Time_rank = dense_rank(Time_num) - 1L) %>%
  ungroup() %>%
  group_by(Dataset, Time_Point) %>%
  arrange(Sample_Name, .by_group = TRUE) %>%
  mutate(rep_id = row_number()) %>%
  ungroup() %>%
  mutate(Sample_Name_New = paste0("t", Time_rank, "_rep", rep_id))

# HCV samples require an explicit biological ordering that cannot be inferred
# from the sample names alone.
hcv_names <- c(
  "Virus_HCV_1_Pre_73dm1" = "t0_rep1",
  "Virus_HCV_1_Pre_73dm2" = "t0_rep2",
  "Virus_HCV_1_Pre_73dm3" = "t0_rep3",
  "Virus_HCV_1_Post_73di1" = "t1_rep1",
  "Virus_HCV_1_Post_73di2" = "t1_rep2",
  "Virus_HCV_1_Post_73di3" = "t1_rep3",
  "Virus_HCV_1_Post_75di1" = "t2_rep1",
  "Virus_HCV_1_Post_75di2" = "t2_rep2",
  "Virus_HCV_1_Post_75di3" = "t2_rep3"
)
hcv_rows <- metadata$Dataset == "Virus_HCV_1" & metadata$Sample_Name %in% names(hcv_names)
metadata$Sample_Name_New[hcv_rows] <- unname(hcv_names[metadata$Sample_Name[hcv_rows]])
metadata <- metadata %>% select(-Time_num, -Time_rank, -rep_id)

duplicate_names <- metadata %>%
  count(Dataset, Sample_Name_New) %>%
  filter(n > 1L)
if (nrow(duplicate_names) > 0L) {
  stop("Sample renaming produced duplicated names within a dataset.", call. = FALSE)
}

fwrite(metadata, FINAL_METADATA_FILE)
message("Wrote standardized metadata: ", FINAL_METADATA_FILE)

gene_bed <- fread(GENE_BED_FILE, header = FALSE, data.table = FALSE)
if (ncol(gene_bed) < 4L) {
  stop("Gene BED must contain at least four columns.", call. = FALSE)
}
names(gene_bed)[1:4] <- c("chr", "start", "end", "info")

gene_length <- gene_bed %>%
  mutate(
    gene_id = str_extract(info, 'gene_id "[^"]+"'),
    gene_id = str_replace_all(gene_id, 'gene_id "|"', ""),
    length = end - start
  ) %>%
  filter(!is.na(gene_id), length > 0) %>%
  select(gene_id, length) %>%
  group_by(gene_id) %>%
  summarise(length = sum(length), .groups = "drop")

gene_length_vector <- gene_length$length
names(gene_length_vector) <- gene_length$gene_id

counts_to_tpm <- function(counts, feature_length) {
  rate <- sweep(counts, 1L, feature_length, "/")
  library_sizes <- colSums(rate, na.rm = TRUE)
  zero_columns <- library_sizes <= 0
  if (any(zero_columns)) {
    warning("One or more samples have a zero effective library size.", call. = FALSE)
    library_sizes[zero_columns] <- NA_real_
  }
  sweep(rate, 2L, library_sizes, "/") * 1e6
}

sample_map_for <- function(dataset_name) {
  map <- metadata %>%
    filter(Dataset == dataset_name) %>%
    select(Sample_Name, Sample_Name_New)
  if (nrow(map) == 0L) {
    stop("No metadata rows found for dataset: ", dataset_name, call. = FALSE)
  }
  map
}

process_gene <- function(file, dataset_name) {
  message("Processing gene counts: ", basename(file))
  mat <- fread(file, data.table = FALSE)
  gene_id <- mat[[1L]]
  expression <- mat[, -1L, drop = FALSE]
  rownames(expression) <- gene_id

  map <- sample_map_for(dataset_name)
  common_samples <- intersect(names(expression), map$Sample_Name)
  if (length(common_samples) == 0L) {
    stop("No gene-count columns match metadata for ", dataset_name, call. = FALSE)
  }
  expression <- expression[, common_samples, drop = FALSE]
  map <- map[match(common_samples, map$Sample_Name), , drop = FALSE]
  names(expression) <- map$Sample_Name_New

  gene_lengths <- gene_length_vector[rownames(expression)]
  keep <- !is.na(gene_lengths) & gene_lengths > 0
  message("  matched gene lengths: ", sum(keep), "/", length(keep))
  expression <- expression[keep, , drop = FALSE]
  gene_lengths <- gene_lengths[keep]

  log2(counts_to_tpm(as.matrix(expression), gene_lengths) + 1)
}

process_te <- function(file, dataset_name) {
  message("Processing TE counts: ", basename(file))
  mat <- fread(file, data.table = FALSE)
  if (ncol(mat) < 7L) {
    stop("TE count matrix must have six annotation columns plus samples.", call. = FALSE)
  }

  te_id <- mat[[1L]]
  te_length <- as.numeric(mat[[6L]])
  expression <- mat[, -(1:6), drop = FALSE]
  rownames(expression) <- te_id

  map <- sample_map_for(dataset_name)
  common_samples <- intersect(names(expression), map$Sample_Name)
  if (length(common_samples) == 0L) {
    stop("No TE-count columns match metadata for ", dataset_name, call. = FALSE)
  }
  expression <- expression[, common_samples, drop = FALSE]
  map <- map[match(common_samples, map$Sample_Name), , drop = FALSE]
  names(expression) <- map$Sample_Name_New

  keep <- !is.na(te_length) & te_length > 0
  expression <- expression[keep, , drop = FALSE]
  te_length <- te_length[keep]

  log2(counts_to_tpm(as.matrix(expression), te_length) + 1)
}

write_expression <- function(matrix, file) {
  output <- data.table(ID = rownames(matrix))
  output <- cbind(output, as.data.table(matrix))
  fwrite(output, file, sep = "\t", quote = FALSE)
}

processed_datasets <- character()
for (dataset in DATASETS) {
  gene_file <- file.path(
    RAW_EXPRESSION_DIR,
    paste0(dataset, "_readscounts_matrix_Gene.txt")
  )
  te_file <- file.path(
    RAW_EXPRESSION_DIR,
    paste0(dataset, "_readscounts_matrix_TE_Loci.txt")
  )

  if (!file.exists(gene_file) || !file.exists(te_file)) {
    warning("Skipping ", dataset, ": one or both count files are missing.", call. = FALSE)
    next
  }

  gene_matrix <- process_gene(gene_file, dataset)
  write_expression(
    gene_matrix,
    file.path(NORMALIZED_EXPRESSION_DIR, paste0(dataset, "_Gene_log2TPM.txt"))
  )

  te_matrix <- process_te(te_file, dataset)
  write_expression(
    te_matrix,
    file.path(NORMALIZED_EXPRESSION_DIR, paste0(dataset, "_TE_log2TPM.txt"))
  )
  processed_datasets <- c(processed_datasets, dataset)
}

if (length(processed_datasets) == 0L) {
  stop(
    "No datasets were normalized. Check RAW_EXPRESSION_DIR and data/README.md.",
    call. = FALSE
  )
}

message("Completed normalization for: ", paste(processed_datasets, collapse = ", "))
