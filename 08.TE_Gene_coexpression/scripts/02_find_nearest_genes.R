#!/usr/bin/env Rscript

# Step 2: find the nearest gene within 50 kb of each intergenic TE.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
project_root <- if (length(script_arg) > 0L) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."))
} else {
  normalizePath(getwd())
}
source(file.path(project_root, "R", "project_setup.R"))

assert_packages(c("data.table", "dplyr"))

library(dplyr)
library(data.table)

assert_files(c(TE_BED_FILE, GENE_BED_FILE, RUTE_LIST_FILE))
ensure_parent_dirs(c(NEAREST_GENE_FILE, TE_GENE_MAP_FILE))
if (!nzchar(Sys.which("bedtools"))) {
  stop("bedtools is required for step 2. See README.md for installation notes.", call. = FALSE)
}

te_file <- TE_BED_FILE
gene_file <- GENE_BED_FILE
ruTE_list_file <- RUTE_LIST_FILE
output_file <- NEAREST_GENE_FILE
simple_output_file <- TE_GENE_MAP_FILE
max_distance <- MAX_DISTANCE

# Restrict the analysis to standard human chromosomes.
standard_chr <- c(paste0("chr", 1:22), "chrX", "chrY")

get_chr_order <- function(chr) {
  chr <- factor(chr, levels = standard_chr)
  return(chr)
}

sort_bed <- function(df, chr_col = "chr", start_col = "start") {
  df %>%
    filter(!!sym(chr_col) %in% standard_chr) %>%
    mutate(chr_ordered = get_chr_order(!!sym(chr_col))) %>%
    arrange(chr_ordered, !!sym(start_col)) %>%
    select(-chr_ordered)
}

# Extract a TE identifier from a GTF-style transcript_id attribute.
parse_te_id <- function(info_string) {
  pattern <- 'transcript_id\\s+([^;]+);'
  match <- regmatches(info_string, regexec(pattern, info_string))
  if(length(match[[1]]) > 1) {
    result <- trimws(match[[1]][2])
    result <- gsub('"', '', result)
    return(result)
  }
  return(NA)
}

# Extract a gene identifier from a GTF-style gene_id attribute.
parse_gene_id <- function(info_string) {
  pattern <- 'gene_id\\s+([^;]+);'
  match <- regmatches(info_string, regexec(pattern, info_string))
  if(length(match[[1]]) > 1) {
    result <- trimws(match[[1]][2])
    result <- gsub('"', '', result)
    return(result)
  }
  return(NA)
}

cat("Reading ruTE list...\n")
ruTE_ids <- fread(ruTE_list_file, header = FALSE)[[1]]
ruTE_ids <- gsub('"', '', ruTE_ids)
ruTE_ids <- trimws(ruTE_ids)
cat(sprintf("Number of ruTEs: %d\n", length(ruTE_ids)))

cat("\nReading TE file...\n")

te_raw <- fread(te_file, header = FALSE, sep = "\t", fill = TRUE)

cat(sprintf("Raw TE file has %d columns\n", ncol(te_raw)))

# Detect the strand column because the supported TE formats place it in
# either column 6 or column 9.
strand_candidates <- intersect(c(6L, 9L), seq_len(ncol(te_raw)))
strand_scores <- vapply(
  strand_candidates,
  function(index) mean(as.character(te_raw[[index]]) %in% c("+", "-", ".")),
  numeric(1)
)
if (length(strand_scores) == 0L || max(strand_scores) < 0.95) {
  stop("Could not identify a valid TE strand column (expected column 6 or 9).")
}
te_strand_column <- strand_candidates[[which.max(strand_scores)]]
cat(sprintf("Using TE strand column: %d\n", te_strand_column))

tes <- data.frame(
  chr = as.character(te_raw[[1]]),
  start = as.numeric(te_raw[[2]]),
  end = as.numeric(te_raw[[3]]),
  te_id_raw = as.character(te_raw[[4]]),
  strand = as.character(te_raw[[te_strand_column]]),
  stringsAsFactors = FALSE
)

tes$te_id <- sapply(tes$te_id_raw, parse_te_id)

tes <- tes %>% filter(!is.na(te_id), te_id != "")

tes <- tes %>% 
  filter(!is.na(chr), !is.na(start), !is.na(end), chr != "", start >= 0, end > start)

tes <- tes %>% filter(chr %in% standard_chr)

tes$te_group <- ifelse(tes$te_id %in% ruTE_ids, "ruTE", "background")

cat(sprintf("Loaded %d TEs (after filtering to standard chromosomes)\n", nrow(tes)))
cat(sprintf("  ruTEs: %d\n", sum(tes$te_group == "ruTE")))
cat(sprintf("  background TEs: %d\n", sum(tes$te_group == "background")))

cat("\nReading Gene file...\n")

genes_raw <- fread(gene_file, header = FALSE, sep = "\t", fill = TRUE)

if(ncol(genes_raw) >= 5) {
  colnames(genes_raw)[1:5] <- c("chr", "start", "end", "info", "strand")
} else {
  stop("Gene file has unexpected number of columns")
}

parse_gene_info <- function(info_string) {
  gene_id <- parse_gene_id(info_string)
  
  pattern_name <- 'gene_name\\s+([^;]+);'
  match_name <- regmatches(info_string, regexec(pattern_name, info_string))
  if(length(match_name[[1]]) > 1) {
    gene_name <- trimws(match_name[[1]][2])
    gene_name <- gsub('"', '', gene_name)
  } else {
    gene_name <- gene_id
  }
  
  if(is.na(gene_id)) gene_id <- "unknown"
  if(is.na(gene_name)) gene_name <- gene_id
  
  return(data.frame(gene_id = gene_id, gene_name = gene_name, stringsAsFactors = FALSE))
}

gene_info <- do.call(rbind, lapply(genes_raw$info, parse_gene_info))

genes <- data.frame(
  chr = as.character(genes_raw$chr),
  start = as.numeric(genes_raw$start),
  end = as.numeric(genes_raw$end),
  gene_id = gene_info$gene_id,
  gene_name = gene_info$gene_name,
  strand = as.character(genes_raw$strand),
  stringsAsFactors = FALSE
)

genes <- genes %>%
  filter(!is.na(chr), !is.na(start), !is.na(end), chr != "", start >= 0, end > start)

genes <- genes %>% filter(chr %in% standard_chr)

# Represent each transcription start site as a one-base BED interval.
genes <- genes %>%
  mutate(
    tss_start = ifelse(strand == "+", start, end - 1),
    tss_end = tss_start + 1
  ) %>%
  filter(!is.na(chr), !is.na(tss_start))

cat(sprintf("Loaded %d genes (after filtering to standard chromosomes)\n", nrow(genes)))

cat("\nSorting BED files...\n")

te_bed <- tes %>%
  select(chr, start, end, te_id, strand) %>%
  rename(name = te_id)

te_bed_sorted <- sort_bed(te_bed, "chr", "start")

gene_bed <- genes %>%
  select(chr, tss_start, tss_end, gene_id, gene_name, strand) %>%
  rename(name = gene_id)

gene_bed_sorted <- sort_bed(gene_bed, "chr", "tss_start")

cat(sprintf("TE file sorted: %d rows\n", nrow(te_bed_sorted)))
cat(sprintf("Gene file sorted: %d rows\n", nrow(gene_bed_sorted)))

temp_te <- tempfile(fileext = ".bed")
temp_gene <- tempfile(fileext = ".bed")

fwrite(te_bed_sorted, temp_te, sep = "\t", col.names = FALSE, quote = FALSE)
fwrite(gene_bed_sorted, temp_gene, sep = "\t", col.names = FALSE, quote = FALSE)

cat("\nRunning bedtools closest...\n")

cmd_sorted <- sprintf("bedtools closest -a %s -b %s -D b -t first -sorted", 
                      temp_te, temp_gene)
cat("Trying with -sorted parameter...\n")

closest <- tryCatch({
  fread(cmd = cmd_sorted, sep = "\t", header = FALSE)
}, error = function(e) {
  cat("Error with -sorted:", e$message, "\n")
  cat("Trying without -sorted parameter...\n")
  cmd_unsorted <- sprintf("bedtools closest -a %s -b %s -D b -t first", 
                          temp_te, temp_gene)
  fread(cmd = cmd_unsorted, sep = "\t", header = FALSE)
})

if(is.null(closest) || nrow(closest) == 0) {
  stop("bedtools closest returned no results. Check input files.")
}

unlink(c(temp_te, temp_gene))

cat(sprintf("bedtools returned %d rows\n", nrow(closest)))

expected_cols <- 12
if(ncol(closest) != expected_cols) {
  cat(sprintf("Warning: Expected %d columns, got %d\n", expected_cols, ncol(closest)))
}

colnames(closest) <- c("te_chr", "te_start", "te_end", "te_id", "te_strand",
                       "gene_chr", "gene_tss_start", "gene_tss_end", "gene_id", "gene_name", "gene_strand",
                       "distance")

# Retain gene annotations only when the closest TSS is within 50 kb.
closest <- closest %>%
  mutate(
    distance = as.numeric(distance),
    abs_distance = abs(distance),
    distance_kb = abs_distance / 1000,
    has_gene_within_50kb = abs_distance <= max_distance
  ) %>%
  mutate(
    gene_id = ifelse(has_gene_within_50kb, gene_id, NA),
    gene_name = ifelse(has_gene_within_50kb, gene_name, NA),
    gene_chr = ifelse(has_gene_within_50kb, gene_chr, NA),
    gene_strand = ifelse(has_gene_within_50kb, gene_strand, NA),
    distance = ifelse(has_gene_within_50kb, distance, NA),
    distance_kb = ifelse(has_gene_within_50kb, distance_kb, NA)
  )

closest <- closest %>%
  left_join(tes %>% select(te_id, te_group), by = "te_id")

output_full <- closest %>%
  select(
    te_chr, te_start, te_end, te_id, te_strand, te_group,
    gene_id, gene_name, gene_chr, gene_strand,
    distance, distance_kb, has_gene_within_50kb
  )

cat(sprintf("\nSaving full results to %s...\n", output_file))
fwrite(output_full, output_file, sep = "\t", na = "NA", quote = FALSE)

cat("\nSaving simple version (te_id, closest_gene_in_50kb, te_group)...\n")

output_simple <- output_full %>%
  select(te_id, closest_gene_in_50kb = gene_id, te_group) %>%
  distinct()

fwrite(output_simple, simple_output_file, sep = "\t", na = "NA", quote = FALSE)
cat(sprintf("Saved simple version to: %s\n", simple_output_file))

cat("\n========================================\n")
cat("Summary Statistics\n")
cat("========================================\n")

cat(sprintf("\nTotal TEs: %d\n", nrow(tes)))
cat(sprintf("  ruTEs: %d\n", sum(tes$te_group == "ruTE")))
cat(sprintf("  background: %d\n", sum(tes$te_group == "background")))

cat(sprintf("\nTEs with nearest gene within 50kb: %d (%.1f%%)\n", 
            sum(output_full$has_gene_within_50kb, na.rm = TRUE),
            100 * sum(output_full$has_gene_within_50kb, na.rm = TRUE) / nrow(output_full)))

cat("\nBy group:\n")
for(g in c("ruTE", "background")) {
  n_te <- sum(output_full$te_group == g, na.rm = TRUE)
  n_with_gene <- sum(output_full$te_group == g & output_full$has_gene_within_50kb, na.rm = TRUE)
  if(n_te > 0) {
    cat(sprintf("  %s: %d/%d (%.1f%%)\n", g, n_with_gene, n_te, 100 * n_with_gene / n_te))
  }
}

cat("\nSimple version summary:\n")
cat(sprintf("  Total rows: %d\n", nrow(output_simple)))
cat(sprintf("  ruTEs with gene: %d\n", 
            sum(output_simple$te_group == "ruTE" & !is.na(output_simple$closest_gene_in_50kb))))
cat(sprintf("  background TEs with gene: %d\n", 
            sum(output_simple$te_group == "background" & !is.na(output_simple$closest_gene_in_50kb))))

cat("\nDone!\n")
cat("  Full output:", output_file, "\n")
cat("  Simple output:", simple_output_file, "\n")
