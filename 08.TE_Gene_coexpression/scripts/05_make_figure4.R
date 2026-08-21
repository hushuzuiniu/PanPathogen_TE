#!/usr/bin/env Rscript

# Step 5: run background resampling and generate Figure 4A-B.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
project_root <- if (length(script_arg) > 0L) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."))
} else {
  normalizePath(getwd())
}
source(file.path(project_root, "R", "project_setup.R"))

assert_packages(c(
  "data.table",
  "doParallel",
  "dplyr",
  "foreach",
  "ggplot2",
  "patchwork"
))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(foreach)
  library(doParallel)
  library(ggplot2)
  library(patchwork)
})

# Input and output paths
te_gene_map_file <- TE_GENE_MAP_FILE
ruTE_results_file <- CORRELATION_RDS
expr_dir <- NORMALIZED_EXPRESSION_DIR

sampling_rds <- FIGURE4_SAMPLING_RDS
null_csv <- FIGURE4_NULL_TSV
summary_rds <- FIGURE4_SUMMARY_RDS
summary_csv <- FIGURE4_SUMMARY_TSV

figureA_pdf <- FIGURE4A_PDF
figureB_pdf <- FIGURE4B_PDF
figureAB_pdf <- FIGURE4AB_PDF

# Analysis parameters
n_permutations <- N_PERMUTATIONS
n_sample <- N_SAMPLE
min_timepoints <- MIN_TIMEPOINTS
n_cores <- N_CORES
random_seed <- RANDOM_SEED

datasets <- DATASETS

thresholds <- THRESHOLDS

ruTE_color <- "#E15B3F"
otherTE_color <- "#6D8FBF"

assert_files(c(te_gene_map_file, ruTE_results_file))
ensure_parent_dirs(c(
  sampling_rds,
  null_csv,
  summary_rds,
  summary_csv,
  figureA_pdf,
  figureB_pdf,
  figureAB_pdf
))

cat("============================================================\n")
cat("Figure 4A-B: 1,000-resample TE-gene co-expression analysis\n")
cat("============================================================\n")
cat(sprintf("n_permutations = %d\n", n_permutations))
cat(sprintf("n_sample = %d\n", n_sample))
cat(sprintf("n_cores = %d\n", n_cores))
cat(sprintf("random_seed = %d\n\n", random_seed))

# Average biological replicates within each time point.
aggregate_replicates <- function(expr_matrix, timepoint_pattern = "_rep") {

  ID_col <- colnames(expr_matrix)[1]
  col_names <- colnames(expr_matrix)[-1]

  timepoints <- gsub(
    paste0(timepoint_pattern, ".*"),
    "",
    col_names
  )

  unique_tps <- unique(timepoints)

  tp_num <- suppressWarnings(
    as.numeric(gsub("^t", "", unique_tps))
  )

  if (all(!is.na(tp_num))) {
    unique_tps <- unique_tps[order(tp_num)]
  }

  result <- data.table(ID = expr_matrix[[ID_col]])

  for (tp in unique_tps) {

    cols <- which(timepoints == tp)

    if (length(cols) == 1) {

      result[[tp]] <- expr_matrix[[cols + 1]]

    } else {

      result[[tp]] <- apply(
        expr_matrix[, cols + 1, with = FALSE],
        1,
        mean,
        na.rm = TRUE
      )
    }
  }

  return(result)
}

cat("Phase 1: Reading the TE-gene map...\n")

if (!file.exists(te_gene_map_file)) {
  stop("TE-gene mapping file not found: ", te_gene_map_file)
}

te_gene_map <- fread(
  te_gene_map_file,
  header = TRUE
)

cat(sprintf(
  "  Loaded %d TE-gene pairs\n",
  nrow(te_gene_map)
))

if (!all(c("te_id", "closest_gene_in_50kb", "te_group") %in%
         colnames(te_gene_map))) {

  stop(
    "Required columns are missing from TE_gene_map_simple.txt. ",
    "Expected: te_id, closest_gene_in_50kb, te_group"
  )
}

cat("\nPhase 2: Building the background pool...\n")

background_pool <- te_gene_map[
  te_group == "background" &
    !is.na(closest_gene_in_50kb) &
    closest_gene_in_50kb != "",
  .(
    te_id,
    gene_id = closest_gene_in_50kb
  )
]

background_pool <- unique(background_pool)

cat(sprintf(
  "  Other intergenic TE loci with a nearest gene within 50 kb: %d\n",
  nrow(background_pool)
))

if (nrow(background_pool) < n_sample) {

  warning(
    "Background pool is smaller than n_sample. ",
    "Sampling will be performed with replacement."
  )

  replace_sample <- TRUE

} else {

  replace_sample <- FALSE
}

if (nrow(background_pool) == 0L) {
  stop("No background TE-gene pairs are available for resampling.")
}

cat(sprintf(
  "  Sampling with replacement: %s\n",
  replace_sample
))

gene_lookup <- setNames(
  background_pool$gene_id,
  background_pool$te_id
)

cat("\nPhase 3: Reading observed ruTE correlations...\n")

if (!file.exists(ruTE_results_file)) {
  stop("ruTE correlation file not found: ", ruTE_results_file)
}

ruTE_results <- readRDS(ruTE_results_file)
setDT(ruTE_results)

if (!all(c("te_id", "rho") %in% colnames(ruTE_results))) {
  stop("ruTE results must contain columns 'te_id' and 'rho'.")
}

ruTE_obs <- ruTE_results[
  !is.na(rho),
  .(
    median_rho = median(rho, na.rm = TRUE)
  ),
  by = te_id
]

ruTE_obs <- ruTE_obs[
  is.finite(median_rho)
]

cat(sprintf(
  "  ruTE loci with valid median correlations: %d\n",
  nrow(ruTE_obs)
))

ruTE_obs_median <- median(
  ruTE_obs$median_rho,
  na.rm = TRUE
)

ruTE_obs_props <- sapply(
  thresholds,
  function(x) mean(
    ruTE_obs$median_rho > x,
    na.rm = TRUE
  )
)

names(ruTE_obs_props) <- paste0(
  "prop_gt_",
  thresholds
)

cat(sprintf(
  "  Observed ruTE median rho: %.3f\n",
  ruTE_obs_median
))

for (i in seq_along(thresholds)) {
  cat(sprintf(
    "  Observed proportion rho > %.1f: %.3f\n",
    thresholds[i],
    ruTE_obs_props[i]
  ))
}

cat("\nPhase 4: Loading expression matrices...\n")

expr_list <- list()

for (ds in datasets) {

  te_file <- file.path(
    expr_dir,
    paste0(ds, "_TE_log2TPM.txt")
  )

  gene_file <- file.path(
    expr_dir,
    paste0(ds, "_Gene_log2TPM.txt")
  )

  if (!file.exists(te_file)) {
    warning("Missing TE expression file: ", te_file)
    next
  }

  if (!file.exists(gene_file)) {
    warning("Missing gene expression file: ", gene_file)
    next
  }

  cat(sprintf("  Loading %s...\n", ds))

  te_raw <- fread(te_file, header = TRUE)
  gene_raw <- fread(gene_file, header = TRUE)

  te_expr <- aggregate_replicates(te_raw)
  gene_expr <- aggregate_replicates(gene_raw)

  setkey(te_expr, ID)
  setkey(gene_expr, ID)

  expr_list[[ds]] <- list(
    te = te_expr,
    gene = gene_expr,
    timepoints = ncol(te_expr) - 1
  )

  cat(sprintf(
    "    %d time points\n",
    expr_list[[ds]]$timepoints
  ))
}

datasets_available <- datasets[
  datasets %in% names(expr_list)
]

if (length(datasets_available) == 0) {
  stop("No expression datasets were successfully loaded.")
}

cat(sprintf(
  "  Successfully loaded %d/%d datasets\n",
  length(datasets_available),
  length(datasets)
))

calc_rho <- function(
    te_id,
    gene_id,
    dataset,
    expr_list,
    min_timepoints = 3) {

  ds_data <- expr_list[[dataset]]

  if (is.null(ds_data)) {
    return(NA_real_)
  }

  if (ds_data$timepoints < min_timepoints) {
    return(NA_real_)
  }

  te_row <- ds_data$te[.(te_id)]
  gene_row <- ds_data$gene[.(gene_id)]

  if (nrow(te_row) == 0 ||
      nrow(gene_row) == 0) {
    return(NA_real_)
  }

  te_vec <- as.numeric(
    te_row[1, -1]
  )

  gene_vec <- as.numeric(
    gene_row[1, -1]
  )

  if (length(te_vec) < min_timepoints ||
      length(gene_vec) < min_timepoints) {
    return(NA_real_)
  }

  if (anyNA(te_vec) ||
      anyNA(gene_vec)) {
    return(NA_real_)
  }

  if (sd(te_vec) == 0 ||
      sd(gene_vec) == 0) {
    return(NA_real_)
  }

  suppressWarnings(
    cor(
      te_vec,
      gene_vec,
      method = "spearman"
    )
  )
}

cat("\nPhase 5: Running background resampling...\n")
cat("This step may take some time.\n")

n_cores <- min(n_cores, n_permutations)
cl <- makeCluster(n_cores)

registerDoParallel(cl)

# Reproducible independent random streams for parallel workers
parallel::clusterSetRNGStream(
  cl,
  iseed = random_seed
)

clusterExport(
  cl,
  varlist = c(
    "background_pool",
    "gene_lookup",
    "expr_list",
    "datasets_available",
    "n_sample",
    "replace_sample",
    "calc_rho",
    "min_timepoints"
  ),
  envir = environment()
)

start_time <- Sys.time()

perm_results <- foreach(
  iter = 1:n_permutations,
  .combine = rbind,
  .packages = "data.table"
) %dopar% {

  sampled_indices <- sample(
    seq_len(nrow(background_pool)),
    size = n_sample,
    replace = replace_sample
  )

  sampled_tes <- background_pool$te_id[
    sampled_indices
  ]

  sampled_genes <- gene_lookup[
    sampled_tes
  ]

  te_rhos <- numeric(0)

  for (i in seq_along(sampled_tes)) {

    te <- sampled_tes[i]
    gene <- sampled_genes[i]

    ds_rhos <- numeric(0)

    for (ds in datasets_available) {

      rho <- calc_rho(
        te,
        gene,
        ds,
        expr_list,
        min_timepoints
      )

      if (!is.na(rho) &&
          is.finite(rho)) {

        ds_rhos <- c(
          ds_rhos,
          rho
        )
      }
    }

    if (length(ds_rhos) > 0) {

      te_rhos <- c(
        te_rhos,
        median(
          ds_rhos,
          na.rm = TRUE
        )
      )
    }
  }

  if (length(te_rhos) > 0) {

    data.table(
      iteration = iter,
      n_valid = length(te_rhos),
      median_rho = median(
        te_rhos,
        na.rm = TRUE
      ),
      mean_rho = mean(
        te_rhos,
        na.rm = TRUE
      ),
      prop_gt_0.3 = mean(
        te_rhos > 0.3,
        na.rm = TRUE
      ),
      prop_gt_0.4 = mean(
        te_rhos > 0.4,
        na.rm = TRUE
      ),
      prop_gt_0.5 = mean(
        te_rhos > 0.5,
        na.rm = TRUE
      ),
      prop_gt_0.6 = mean(
        te_rhos > 0.6,
        na.rm = TRUE
      ),
      prop_gt_0.7 = mean(
        te_rhos > 0.7,
        na.rm = TRUE
      )
    )

  } else {

    data.table(
      iteration = iter,
      n_valid = 0,
      median_rho = NA_real_,
      mean_rho = NA_real_,
      prop_gt_0.3 = NA_real_,
      prop_gt_0.4 = NA_real_,
      prop_gt_0.5 = NA_real_,
      prop_gt_0.6 = NA_real_,
      prop_gt_0.7 = NA_real_
    )
  }
}

stopCluster(cl)

elapsed_minutes <- as.numeric(
  difftime(
    Sys.time(),
    start_time,
    units = "mins"
  )
)

cat(sprintf(
  "\n  Resampling completed in %.2f minutes\n",
  elapsed_minutes
))

n_valid_iters <- sum(
  is.finite(perm_results$median_rho)
)

cat(sprintf(
  "  Valid resamples: %d/%d\n",
  n_valid_iters,
  n_permutations
))

cat("\nPhase 6: Calculating null distributions and P values...\n")

valid_perm <- perm_results[
  is.finite(median_rho)
]

expected_stats <- list(
  mean_median_rho = mean(
    valid_perm$median_rho
  ),
  sd_median_rho = sd(
    valid_perm$median_rho
  ),
  mean_prop_gt_0.3 = mean(
    valid_perm$prop_gt_0.3,
    na.rm = TRUE
  ),
  sd_prop_gt_0.3 = sd(
    valid_perm$prop_gt_0.3,
    na.rm = TRUE
  ),
  mean_prop_gt_0.4 = mean(
    valid_perm$prop_gt_0.4,
    na.rm = TRUE
  ),
  sd_prop_gt_0.4 = sd(
    valid_perm$prop_gt_0.4,
    na.rm = TRUE
  ),
  mean_prop_gt_0.5 = mean(
    valid_perm$prop_gt_0.5,
    na.rm = TRUE
  ),
  sd_prop_gt_0.5 = sd(
    valid_perm$prop_gt_0.5,
    na.rm = TRUE
  ),
  mean_prop_gt_0.6 = mean(
    valid_perm$prop_gt_0.6,
    na.rm = TRUE
  ),
  sd_prop_gt_0.6 = sd(
    valid_perm$prop_gt_0.6,
    na.rm = TRUE
  ),
  mean_prop_gt_0.7 = mean(
    valid_perm$prop_gt_0.7,
    na.rm = TRUE
  ),
  sd_prop_gt_0.7 = sd(
    valid_perm$prop_gt_0.7,
    na.rm = TRUE
  )
)

# Plus-one corrected one-sided empirical P value
empirical_p <- function(
    null_values,
    observed_value) {

  null_values <- null_values[
    is.finite(null_values)
  ]

  (
    sum(
      null_values >= observed_value
    ) + 1
  ) / (
    length(null_values) + 1
  )
}

p_values <- data.table(
  metric = c(
    "median_rho",
    "prop_gt_0.3",
    "prop_gt_0.4",
    "prop_gt_0.5",
    "prop_gt_0.6",
    "prop_gt_0.7"
  ),
  observed = c(
    ruTE_obs_median,
    ruTE_obs_props[1],
    ruTE_obs_props[2],
    ruTE_obs_props[3],
    ruTE_obs_props[4],
    ruTE_obs_props[5]
  ),
  p_value = c(
    empirical_p(
      valid_perm$median_rho,
      ruTE_obs_median
    ),
    empirical_p(
      valid_perm$prop_gt_0.3,
      ruTE_obs_props[1]
    ),
    empirical_p(
      valid_perm$prop_gt_0.4,
      ruTE_obs_props[2]
    ),
    empirical_p(
      valid_perm$prop_gt_0.5,
      ruTE_obs_props[3]
    ),
    empirical_p(
      valid_perm$prop_gt_0.6,
      ruTE_obs_props[4]
    ),
    empirical_p(
      valid_perm$prop_gt_0.7,
      ruTE_obs_props[5]
    )
  )
)

cat("\nObserved vs background:\n")
cat(sprintf(
  "  Median rho: %.3f vs %.3f ± %.3f; p = %.4g\n",
  ruTE_obs_median,
  expected_stats$mean_median_rho,
  expected_stats$sd_median_rho,
  p_values[metric == "median_rho", p_value]
))

for (i in seq_along(thresholds)) {

  metric_name <- paste0(
    "prop_gt_",
    thresholds[i]
  )

  cat(sprintf(
    "  rho > %.1f: %.1f%% vs %.1f%% ± %.1f%%; p = %.4g\n",
    thresholds[i],
    100 * ruTE_obs_props[i],
    100 * expected_stats[[paste0("mean_", metric_name)]],
    100 * expected_stats[[paste0("sd_", metric_name)]],
    p_values[metric == metric_name, p_value]
  ))
}

cat("\nPhase 7: Saving statistical results...\n")

null_dist <- valid_perm[
  ,
  .(
    iteration,
    n_valid,
    median_rho,
    prop_gt_0.3,
    prop_gt_0.4,
    prop_gt_0.5,
    prop_gt_0.6,
    prop_gt_0.7
  )
]

saveRDS(
  perm_results,
  sampling_rds
)

fwrite(
  null_dist,
  null_csv,
  sep = "\t",
  quote = FALSE
)

summary_output <- list(
  n_permutations = n_permutations,
  n_sample = n_sample,
  n_valid_iterations = n_valid_iters,
  background_pool_size = nrow(background_pool),
  ruTE_n_valid = nrow(ruTE_obs),
  ruTE_observed = data.table(
    metric = c(
      "median_rho",
      "prop_gt_0.3",
      "prop_gt_0.4",
      "prop_gt_0.5",
      "prop_gt_0.6",
      "prop_gt_0.7"
    ),
    value = c(
      ruTE_obs_median,
      ruTE_obs_props
    )
  ),
  expected_stats = expected_stats,
  p_values = p_values
)

saveRDS(
  summary_output,
  summary_rds
)

fwrite(
  p_values,
  summary_csv,
  sep = "\t",
  quote = FALSE
)

cat("\nPhase 8: Creating Figure 4A...\n")

background_medians <- valid_perm[
  ,
  .(
    iteration,
    median_rho
  )
]

perm_p <- p_values[
  metric == "median_rho",
  p_value
]

sig_label <- ifelse(
  perm_p < 0.001,
  "***",
  ifelse(
    perm_p < 0.01,
    "**",
    ifelse(
      perm_p < 0.05,
      "*",
      "ns"
    )
  )
)

background_x <- 1
ruTE_x <- 2

all_y <- c(
  background_medians$median_rho,
  ruTE_obs_median
)

y_min <- min(
  all_y,
  na.rm = TRUE
)

y_max_data <- max(
  all_y,
  na.rm = TRUE
)

y_range <- max(
  y_max_data - y_min,
  0.1
)

y_upper <- y_max_data + 0.28 * y_range
y_lower <- y_min - 0.05 * y_range

bracket_y <- y_max_data + 0.13 * y_range
bracket_height <- 0.018 * y_range

pA <- ggplot(
  background_medians,
  aes(
    x = background_x,
    y = median_rho
  )
) +

  geom_boxplot(
    width = 0.42,
    fill = otherTE_color,
    color = otherTE_color,
    alpha = 0.35,
    outlier.shape = NA,
    linewidth = 0.65
  ) +

  geom_jitter(
    width = 0.12,
    height = 0,
    size = 1.25,
    alpha = 0.55,
    color = otherTE_color
  ) +

  geom_point(
    data = data.frame(
      x = ruTE_x,
      median_rho = ruTE_obs_median
    ),
    aes(
      x = x,
      y = median_rho
    ),
    shape = 18,
    size = 5.2,
    color = ruTE_color,
    inherit.aes = FALSE
  ) +

  geom_segment(
    aes(
      x = background_x,
      xend = background_x,
      y = bracket_y,
      yend = bracket_y - bracket_height
    ),
    linewidth = 0.45,
    inherit.aes = FALSE
  ) +

  geom_segment(
    aes(
      x = background_x,
      xend = ruTE_x,
      y = bracket_y,
      yend = bracket_y
    ),
    linewidth = 0.45,
    inherit.aes = FALSE
  ) +

  geom_segment(
    aes(
      x = ruTE_x,
      xend = ruTE_x,
      y = bracket_y,
      yend = bracket_y - bracket_height
    ),
    linewidth = 0.45,
    inherit.aes = FALSE
  ) +

  annotate(
    "text",
    x = 1.5,
    y = bracket_y + 0.045 * y_range,
    label = sig_label,
    size = 4.2,
    fontface = "bold"
  ) +

  annotate(
    "text",
    x = 1.5,
    y = bracket_y + 0.005 * y_range,
    label = paste0(
      "p=",
      format(
        perm_p,
        digits = 3,
        scientific = FALSE
      )
    ),
    size = 3.2,
    fontface = "italic"
  ) +

  annotate(
    "text",
    x = ruTE_x,
    y = ruTE_obs_median - 0.055 * y_range,
    label = paste0(
      "rho==",
      sprintf(
        "%.3f",
        ruTE_obs_median
      )
    ),
    parse = TRUE,
    color = ruTE_color,
    size = 3.3,
    fontface = "bold"
  ) +

  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "(a)",
    size = 6.5,
    fontface = "bold",
    hjust = -0.15,
    vjust = 1.25
  ) +

  scale_x_continuous(
    breaks = c(
      background_x,
      ruTE_x
    ),
    labels = c(
      paste0(
        "Other intergenic TEs\n(",
        n_permutations,
        " resamples)"
      ),
      paste0(
        "ruTE\n(",
        nrow(ruTE_obs),
        " loci)"
      )
    ),
    limits = c(
      0.55,
      2.45
    ),
    expand = c(0, 0)
  ) +

  scale_y_continuous(
    limits = c(
      y_lower,
      y_upper
    ),
    expand = c(0, 0)
  ) +

  labs(
    x = "",
    y = "Median Spearman Correlation Coefficient (ρ)"
  ) +

  theme_bw(
    base_size = 11
  ) +

  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.x = element_text(
      face = "bold",
      size = 10
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 10
    ),
    axis.text.x = element_text(
      face = "bold",
      size = 9.5,
      lineheight = 0.95
    ),
    axis.text.y = element_text(
      size = 9
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(
      color = "gray90",
      linewidth = 0.3
    ),
    panel.border = element_blank(),
    axis.line = element_line(
      color = "black",
      linewidth = 0.45
    ),
    plot.margin = margin(
      8, 8, 5, 8
    )
  )

ggsave(
  figureA_pdf,
  pA,
  width = 3.45,
  height = 4.1,
  units = "in"
)

cat("\nPhase 9: Creating Figure 4B...\n")

# Observed ruTE proportions above each correlation threshold.
threshold_data <- data.table(
  threshold = thresholds,
  observed = as.numeric(ruTE_obs_props)
)

# Mean and standard deviation from the background resamples.
bg_props <- rbindlist(
  lapply(thresholds, function(x) {
    
    column_name <- paste0(
      "prop_gt_",
      format(x, nsmall = 1, trim = TRUE)
    )
    
    vals <- valid_perm[[column_name]]
    vals <- vals[is.finite(vals)]
    
    data.table(
      threshold = x,
      observed = mean(vals),
      sd = sd(vals),
      n_valid = length(vals)
    )
  })
)

print(threshold_data)
print(bg_props)
threshold_pvals <- p_values[
  metric %in% paste0(
    "prop_gt_",
    thresholds
  )
]

threshold_pvals[, threshold := as.numeric(
  sub(
    "prop_gt_",
    "",
    metric
  )
)]

threshold_pvals[, sig_label := fifelse(
  p_value < 0.001,
  "***",
  fifelse(
    p_value < 0.01,
    "**",
    fifelse(
      p_value < 0.05,
      "*",
      "ns"
    )
  )
)]

enrichment <- merge(
  threshold_data,
  bg_props,
  by = "threshold",
  suffixes = c(
    "_ruTE",
    "_background"
  )
)

enrichment[, fold_enrichment :=
            observed_ruTE /
            observed_background]

bar_data <- rbind(
  threshold_data[
    ,
    .(
      threshold,
      proportion = observed,
      group = "ruTE"
    )
  ],
  bg_props[
    ,
    .(
      threshold,
      proportion = observed,
      group = "background"
    )
  ]
)

bar_data[, x_base := match(
  threshold,
  thresholds
)]

bar_data[
  group == "ruTE",
  x_pos := x_base - 0.20
]

bar_data[
  group == "background",
  x_pos := x_base + 0.20
]

y_max_bar_data <- max(
  c(
    bar_data$proportion +
      bar_data$proportion * 0.05,
    bg_props$observed +
      bg_props$sd
  ),
  na.rm = TRUE
)

y_max_bar <- max(
  0.65,
  y_max_bar_data * 1.18
)

sig_y <- y_max_bar * 0.91

pB <- ggplot() +

  geom_col(
    data = bar_data,
    aes(
      x = x_pos,
      y = proportion,
      fill = group
    ),
    width = 0.32,
    alpha = 0.9
  ) +

  geom_errorbar(
    data = bg_props,
    aes(
      x = match(
        threshold,
        thresholds
      ) + 0.20,
      ymin = pmax(
        0,
        observed - sd
      ),
      ymax = observed + sd
    ),
    width = 0.10,
    linewidth = 0.4
  ) +

  geom_text(
    data = bar_data,
    aes(
      x = x_pos,
      y = proportion + 0.018,
      label = paste0(
        round(
          proportion * 100,
          0
        ),
        "%"
      )
    ),
    size = 3.0
  ) +

  geom_text(
    data = enrichment,
    aes(
      x = match(
        threshold,
        thresholds
      ) - 0.20,
      y = observed_ruTE + 0.065,
      label = paste0(
        round(
          fold_enrichment,
          1
        ),
        "×"
      )
    ),
    size = 3.0,
    fontface = "bold",
    color = ruTE_color
  ) +

  geom_segment(
    data = threshold_pvals,
    aes(
      x = match(
        threshold,
        thresholds
      ) - 0.28,
      xend = match(
        threshold,
        thresholds
      ) + 0.28,
      y = sig_y,
      yend = sig_y
    ),
    linewidth = 0.4
  ) +

  geom_segment(
    data = threshold_pvals,
    aes(
      x = match(
        threshold,
        thresholds
      ) - 0.28,
      xend = match(
        threshold,
        thresholds
      ) - 0.28,
      y = sig_y - 0.018,
      yend = sig_y
    ),
    linewidth = 0.4
  ) +

  geom_segment(
    data = threshold_pvals,
    aes(
      x = match(
        threshold,
        thresholds
      ) + 0.28,
      xend = match(
        threshold,
        thresholds
      ) + 0.28,
      y = sig_y - 0.018,
      yend = sig_y
    ),
    linewidth = 0.4
  ) +

  geom_text(
    data = threshold_pvals,
    aes(
      x = match(
        threshold,
        thresholds
      ),
      y = sig_y + 0.025,
      label = sig_label
    ),
    size = 3.7,
    fontface = "bold"
  ) +

  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "(b)",
    size = 6.5,
    fontface = "bold",
    hjust = -0.15,
    vjust = 1.25
  ) +

  scale_fill_manual(
    values = c(
      "ruTE" = ruTE_color,
      "background" = otherTE_color
    ),
    labels = c(
      "ruTE" = "ruTEs",
      "background" = "Other TEs"
    )
  ) +

  scale_x_continuous(
    breaks = seq_along(
      thresholds
    ),
    labels = paste0(
      ">",
      thresholds
    ),
    limits = c(
      0.45,
      length(thresholds) + 0.55
    ),
    expand = c(0, 0)
  ) +

  scale_y_continuous(
    labels = scales::percent,
    limits = c(
      0,
      y_max_bar
    ),
    expand = c(0, 0)
  ) +

  labs(
    x = "Median Spearman Correlation Coefficient (ρ)",
    y = "Proportion of TEs",
    fill = ""
  ) +

  theme_bw(
    base_size = 10.5
  ) +

  theme(
    legend.position = c(
      0.82,
      0.91
    ),
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(
      size = 8.5
    ),
    legend.key.height = grid::unit(
      0.35,
      "cm"
    ),
    legend.background = element_blank(),
    axis.title.x = element_text(
      face = "bold",
      size = 9.5
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 9.5
    ),
    axis.text.x = element_text(
      size = 9
    ),
    axis.text.y = element_text(
      size = 9
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(
      color = "gray90",
      linewidth = 0.3
    ),
    panel.border = element_blank(),
    axis.line = element_line(
      color = "black",
      linewidth = 0.45
    ),
    plot.margin = margin(
      8, 8, 5, 8
    )
  )

ggsave(
  figureB_pdf,
  pB,
  width = 6.2,
  height = 4.1,
  units = "in"
)

cat("\nPhase 10: Creating combined Figure 4A-B...\n")

pAB <- pA +
  pB +
  plot_layout(
    widths = c(
      0.95,
      1.75
    )
  )

ggsave(
  figureAB_pdf,
  pAB,
  width = 9.2,
  height = 4.1,
  units = "in"
)

cat("\n============================================================\n")
cat("Completed successfully.\n")
cat("============================================================\n")

cat("\nObserved ruTE statistics:\n")
cat(sprintf(
  "  n = %d loci\n",
  nrow(ruTE_obs)
))
cat(sprintf(
  "  Median rho = %.3f\n",
  ruTE_obs_median
))

for (i in seq_along(thresholds)) {

  cat(sprintf(
    "  rho > %.1f: %.1f%%\n",
    thresholds[i],
    100 * ruTE_obs_props[i]
  ))
}

cat("\nEmpirical P values:\n")

for (i in seq_len(nrow(p_values))) {

  cat(sprintf(
    "  %-15s p = %.6g\n",
    p_values$metric[i],
    p_values$p_value[i]
  ))
}

cat("\nOutput files:\n")
cat("  ", sampling_rds, "\n", sep = "")
cat("  ", null_csv, "\n", sep = "")
cat("  ", summary_rds, "\n", sep = "")
cat("  ", summary_csv, "\n", sep = "")
cat("  ", figureA_pdf, "\n", sep = "")
cat("  ", figureB_pdf, "\n", sep = "")
cat("  ", figureAB_pdf, "\n", sep = "")

cat("\nDone.\n")
