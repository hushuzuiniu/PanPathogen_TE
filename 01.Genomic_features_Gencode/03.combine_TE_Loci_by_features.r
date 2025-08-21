library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

setwd("/data2t_2/hushu/01.Genomic_features_Gencode_v2/")
result <- fread("hg38_TE_anno_custom_v20240110_0based_with_feature_summary_v2.bed")
TE_Loci_info <- result[, c("gene_id", "transcript_id", "class_id", "family_id", "prioritized_feature")]
TE_Loci_info$region <- TE_Loci_info$prioritized_feature
TE_Loci_info$region[TE_Loci_info$prioritized_feature %in% 
                      c("3UTR", "5UTR", "CDS", "noncoding_exon")] <- "exon"

input_dir <- "/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data/"
output_dir <- "/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data_TE_subfamily_feature_v2/"

files <- list.files(path = input_dir, pattern = ".*readscounts_matrix_TE_Loci\\.txt$", full.names = TRUE)

transcript_mapping <- TE_Loci_info %>%
  select(gene_id, transcript_id, region)


for(file in files) {
  file_name <- basename(file)
  output_file_name <- gsub("matrix_TE_Loci", "matrix_TE_subfamily_feature", file_name)
  
  prefix <- str_extract(file_name, "^[^_]+_[^_]+_[^_]+")
  if(is.na(prefix)) {
    prefix <- str_extract(file_name, "^[^_]+")
  }
  raw_counts <- read.csv(file, sep = "\t")
  merged_data <- raw_counts %>%
    rename(transcript_id = Geneid) %>%
    inner_join(transcript_mapping, by = "transcript_id")
  
  expression_cols <- names(raw_counts)[grep(paste0("^", prefix), names(raw_counts))]
  
  if(length(expression_cols) == 0) {
    meta_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
    expression_cols <- setdiff(names(raw_counts), meta_cols)
  }
  
  result <- merged_data %>%
    group_by(gene_id, region) %>%
    summarize(across(all_of(expression_cols), sum),
              .groups = "drop")
  
  result <- result %>%
    mutate(gene_id_region = paste0(gene_id, "_", region))
  
  result <- result %>%
    select(gene_id_region, all_of(expression_cols))
  
  names(result)[1] <- "GeneID"
  output_file <- file.path(output_dir, output_file_name)
  write.table(result, output_file, row.names = FALSE, sep = "\t", quote = FALSE)
}

##################################################################################################
input_dir <- "/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data/"
output_dir <- "/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data_TE_subfamily/"

files <- list.files(path = input_dir, pattern = ".*readscounts_matrix_TE_Loci\\.txt$", full.names = TRUE)

transcript_mapping <- TE_Loci_info %>%
  select(gene_id, transcript_id)

for(file in files) {
  file_name <- basename(file)
  output_file_name <- gsub("matrix_TE_Loci", "matrix_TE_subfamily", file_name)
  
  prefix <- str_extract(file_name, "^[^_]+_[^_]+_[^_]+")
  if(is.na(prefix)) {
    prefix <- str_extract(file_name, "^[^_]+")
  }
  
  raw_counts <- read.csv(file, sep = "\t")
  
  merged_data <- raw_counts %>%
    rename(transcript_id = Geneid) %>%
    inner_join(transcript_mapping, by = "transcript_id")
  
  expression_cols <- names(raw_counts)[grep(paste0("^", prefix), names(raw_counts))]
  
  if(length(expression_cols) == 0) {
    meta_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
    expression_cols <- setdiff(names(raw_counts), meta_cols)
  }
  
  result <- merged_data %>%
    group_by(gene_id) %>%
    summarize(across(all_of(expression_cols), sum),
              .groups = "drop")
  
  result <- result %>%
    select(gene_id, all_of(expression_cols))
  
  names(result)[1] <- "GeneID"
  
  output_file <- file.path(output_dir, output_file_name)
  write.table(result, output_file, row.names = FALSE, sep = "\t", quote = FALSE)
}