# Shared configuration for the TE-gene co-expression workflow.
#
# Run scripts from the repository root. Paths can be overridden with the
# environment variables documented in README.md.

env_integer <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(as.integer(default))
  }

  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1L) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  parsed
}

env_path <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

DATASETS <- c(
  "Virus_EBV_2",
  "Virus_EBV_5",
  "Fungi_Af_4",
  "Fungi_Af_6",
  "Virus_IAV_11",
  "Virus_RSV_1",
  "Virus_HCV_1"
)

DATA_DIR <- file.path(PROJECT_ROOT, "data")
RESULTS_DIR <- file.path(PROJECT_ROOT, "results")

SAMPLE_METADATA_FILE <- env_path(
  "TE_COEXPR_METADATA",
  file.path(DATA_DIR, "reference", "sample_info_7_datasets_with_timecourse.csv")
)
GENE_BED_FILE <- env_path(
  "TE_COEXPR_GENE_BED",
  file.path(DATA_DIR, "reference", "Genes_0based.bed")
)
TE_BED_FILE <- env_path(
  "TE_COEXPR_TE_BED",
  file.path(DATA_DIR, "reference", "TEs_all_other_intergenic_0based.bed")
)
RUTE_LIST_FILE <- env_path(
  "TE_COEXPR_RUTE_LIST",
  file.path(DATA_DIR, "reference", "ruTE_other_intergenic_n595.txt")
)
RAW_EXPRESSION_DIR <- env_path(
  "TE_COEXPR_RAW_EXPRESSION",
  file.path(DATA_DIR, "raw_expression")
)
NORMALIZED_EXPRESSION_DIR <- env_path(
  "TE_COEXPR_NORMALIZED_EXPRESSION",
  file.path(DATA_DIR, "normalized_expression")
)

FINAL_METADATA_FILE <- file.path(
  DATA_DIR,
  "processed",
  "sample_info_7_datasets_with_timecourse_final.csv"
)
NEAREST_GENE_FILE <- file.path(
  RESULTS_DIR,
  "nearest_gene",
  "TEs_with_nearest_gene.tsv"
)
TE_GENE_MAP_FILE <- file.path(
  RESULTS_DIR,
  "nearest_gene",
  "TE_gene_map_simple.tsv"
)

CORRELATION_RDS <- file.path(
  RESULTS_DIR,
  "correlation",
  "rute_correlations_raw.rds"
)
CORRELATION_TSV <- file.path(
  RESULTS_DIR,
  "correlation",
  "rute_correlations_raw.tsv"
)
CORRELATION_SUMMARY_TSV <- file.path(
  RESULTS_DIR,
  "correlation",
  "rute_correlation_summary.tsv"
)

FIGURE4_SAMPLING_RDS <- file.path(
  RESULTS_DIR,
  "figure4",
  "figure4_1000resamples_sampling_results.rds"
)
FIGURE4_NULL_TSV <- file.path(
  RESULTS_DIR,
  "figure4",
  "figure4_1000resamples_null_distribution.tsv"
)
FIGURE4_SUMMARY_RDS <- file.path(
  RESULTS_DIR,
  "figure4",
  "figure4_1000resamples_summary_stats.rds"
)
FIGURE4_SUMMARY_TSV <- file.path(
  RESULTS_DIR,
  "figure4",
  "figure4_1000resamples_summary_stats.tsv"
)
FIGURE4A_PDF <- file.path(RESULTS_DIR, "figure4", "Figure4A_1000resamples.pdf")
FIGURE4B_PDF <- file.path(RESULTS_DIR, "figure4", "Figure4B_1000resamples.pdf")
FIGURE4AB_PDF <- file.path(RESULTS_DIR, "figure4", "Figure4AB_1000resamples.pdf")

MAX_DISTANCE <- 50000L
MIN_TIMEPOINTS <- env_integer("TE_COEXPR_MIN_TIMEPOINTS", 3L)
detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores) || detected_cores < 1L) {
  detected_cores <- 1L
}
N_CORES <- env_integer(
  "TE_COEXPR_N_CORES",
  max(1L, min(8L, detected_cores))
)
N_PERMUTATIONS <- env_integer("TE_COEXPR_N_PERMUTATIONS", 1000L)
N_SAMPLE <- env_integer("TE_COEXPR_N_SAMPLE", 595L)
RANDOM_SEED <- env_integer("TE_COEXPR_RANDOM_SEED", 123L)
THRESHOLDS <- c(0.3, 0.4, 0.5, 0.6, 0.7)
