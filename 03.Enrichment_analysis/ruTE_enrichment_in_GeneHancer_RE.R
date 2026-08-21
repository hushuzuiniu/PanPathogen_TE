#!/usr/bin/env Rscript

# Test whether ruTE loci are enriched in elite GeneHancer regulatory elements.
#
# Usage:
#   Rscript ruTE_enrichment_in_GeneHancer_RE_analysis_v2.R \
#     [te_table] [genehancer_table] [output_directory]
#
# Both input tables are expected to use zero-based, half-open coordinates.

required_packages <- c("data.table", "GenomicRanges", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
})

get_script_directory <- function() {
  script_argument <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )
  if (length(script_argument) == 0L) {
    return(normalizePath(getwd(), mustWork = TRUE))
  }
  
  script_path <- sub("^--file=", "", script_argument[[1L]])
  dirname(normalizePath(script_path, mustWork = TRUE))
}

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) > 3L) {
  stop(
    paste(
      "Usage: Rscript ruTE_enrichment_in_GeneHancer_RE_analysis_v2.R",
      "[te_table] [genehancer_table] [output_directory]"
    ),
    call. = FALSE
  )
}

script_directory <- get_script_directory()
te_file <- if (length(arguments) >= 1L) {
  normalizePath(arguments[[1L]], mustWork = TRUE)
} else {
  file.path(script_directory, "TEs_with_nearest_gene_with_n838_ruTEs.txt")
}
genehancer_file <- if (length(arguments) >= 2L) {
  normalizePath(arguments[[2L]], mustWork = TRUE)
} else {
  file.path(script_directory, "GeneHancer_AnnotSV_elements_v5.24.txt")
}
output_directory <- if (length(arguments) >= 3L) {
  arguments[[3L]]
} else {
  file.path(script_directory, "output")
}

missing_files <- c(te_file, genehancer_file)[
  !file.exists(c(te_file, genehancer_file))
]
if (length(missing_files) > 0L) {
  stop(
    "Missing input file(s):\n- ",
    paste(missing_files, collapse = "\n- "),
    call. = FALSE
  )
}
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
output_directory <- normalizePath(output_directory, mustWork = TRUE)

message("Reading TE and GeneHancer tables...")
te <- fread(te_file, na.strings = c("NA", ""))
regulatory_elements <- fread(genehancer_file, na.strings = c("NA", ""))

required_te_columns <- c("te_chr", "te_start", "te_end", "te_group")
required_regulatory_columns <- c(
  "chr",
  "element_start",
  "element_end",
  "is_elite"
)
missing_te_columns <- setdiff(required_te_columns, names(te))
missing_regulatory_columns <- setdiff(
  required_regulatory_columns,
  names(regulatory_elements)
)
if (length(missing_te_columns) > 0L) {
  stop(
    "The TE table is missing columns: ",
    paste(missing_te_columns, collapse = ", "),
    call. = FALSE
  )
}
if (length(missing_regulatory_columns) > 0L) {
  stop(
    "The GeneHancer table is missing columns: ",
    paste(missing_regulatory_columns, collapse = ", "),
    call. = FALSE
  )
}

te[, te_group := trimws(te_group)]
te <- te[te_group %in% c("ruTE", "background")]
te[, `:=`(
  te_start = as.numeric(te_start),
  te_end = as.numeric(te_end)
)]
te <- te[
  !is.na(te_chr) & !is.na(te_start) & !is.na(te_end) &
    te_start >= 0 & te_end > te_start
]
te[, te_length := te_end - te_start]
te[, te_index := .I]
if (nrow(te) == 0L) {
  stop("No valid ruTE or background TE records remain after filtering.")
}
if (!all(c("ruTE", "background") %in% te$te_group)) {
  stop("The TE table must contain both ruTE and background groups.")
}

regulatory_elements <- regulatory_elements[is_elite == 1]
regulatory_elements[, `:=`(
  element_start = as.numeric(element_start),
  element_end = as.numeric(element_end)
)]
regulatory_elements <- regulatory_elements[
  !is.na(chr) & !is.na(element_start) & !is.na(element_end) &
    element_start >= 0 & element_end > element_start
]
if (nrow(regulatory_elements) == 0L) {
  stop("No valid elite GeneHancer elements remain after filtering.")
}

normalize_chromosome <- function(chromosome) {
  chromosome <- as.character(chromosome)
  ifelse(grepl("^chr", chromosome), chromosome, paste0("chr", chromosome))
}

te_ranges <- GRanges(
  seqnames = normalize_chromosome(te$te_chr),
  ranges = IRanges(start = te$te_start + 1, end = te$te_end)
)
regulatory_ranges <- GRanges(
  seqnames = normalize_chromosome(regulatory_elements$chr),
  ranges = IRanges(
    start = regulatory_elements$element_start + 1,
    end = regulatory_elements$element_end
  )
)

message("Calculating TE-regulatory element overlaps...")
hits <- findOverlaps(te_ranges, regulatory_ranges, ignore.strand = TRUE)

if (length(hits) > 0L) {
  overlap_width <- width(
    pintersect(
      te_ranges[queryHits(hits)],
      regulatory_ranges[subjectHits(hits)]
    )
  )
  overlap_table <- data.table(
    te_index = queryHits(hits),
    regulatory_index = subjectHits(hits),
    overlap_bp = overlap_width
  )
  
  te_metadata_columns <- setdiff(names(te), "te_index")
  overlap_table <- cbind(
    overlap_table,
    te[overlap_table$te_index, ..te_metadata_columns],
    regulatory_elements[overlap_table$regulatory_index]
  )
  overlap_table[, overlap_fraction := overlap_bp / te_length]
  
  overlap_summary <- overlap_table[, .(
    max_overlap_bp = max(overlap_bp),
    max_overlap_fraction = max(overlap_fraction)
  ), by = te_index]
} else {
  overlap_table <- data.table(
    te_index = integer(),
    regulatory_index = integer(),
    overlap_bp = integer(),
    overlap_fraction = numeric()
  )
  overlap_summary <- data.table(
    te_index = integer(),
    max_overlap_bp = numeric(),
    max_overlap_fraction = numeric()
  )
}

te_summary <- merge(te, overlap_summary, by = "te_index", all.x = TRUE, sort = FALSE)
setorder(te_summary, te_index)
te_summary[is.na(max_overlap_bp), `:=`(
  max_overlap_bp = 0,
  max_overlap_fraction = 0
)]
te_summary[, `:=`(
  overlap_any = max_overlap_bp > 0,
  overlap_50 = max_overlap_fraction >= 0.5,
  overlap_100 = max_overlap_fraction >= 1
)]

fwrite(
  overlap_table,
  file.path(output_directory, "TE_regulatory_overlap_all.tsv"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

run_fisher_test <- function(overlap_flag, label) {
  contingency_table <- table(
    factor(te_summary$te_group, levels = c("ruTE", "background")),
    factor(as.logical(overlap_flag), levels = c(TRUE, FALSE))
  )
  colnames(contingency_table) <- c("Overlap", "No_overlap")
  rownames(contingency_table) <- c("ruTE", "background")
  
  fisher_result <- fisher.test(contingency_table)
  data.table(
    threshold = label,
    odds_ratio = unname(fisher_result$estimate),
    ci_lower = fisher_result$conf.int[[1L]],
    ci_upper = fisher_result$conf.int[[2L]],
    p_value = fisher_result$p.value,
    rute_overlap = contingency_table["ruTE", "Overlap"],
    rute_no_overlap = contingency_table["ruTE", "No_overlap"],
    background_overlap = contingency_table["background", "Overlap"],
    background_no_overlap = contingency_table["background", "No_overlap"]
  )
}

message("Running Fisher's exact tests...")
statistics <- rbindlist(list(
  run_fisher_test(te_summary$overlap_any, "Any"),
  run_fisher_test(te_summary$overlap_50, "50%"),
  run_fisher_test(te_summary$overlap_100, "100%")
))
statistics[, threshold := factor(threshold, levels = c("Any", "50%", "100%"))]

fwrite(
  statistics,
  file.path(output_directory, "TE_enrichment_statistics.tsv"),
  sep = "\t",
  quote = FALSE
)

contingency_tables <- rbindlist(lapply(seq_len(nrow(statistics)), function(index) {
  data.table(
    threshold = statistics$threshold[index],
    group = c("ruTE", "background"),
    overlap = c(
      statistics$rute_overlap[index],
      statistics$background_overlap[index]
    ),
    no_overlap = c(
      statistics$rute_no_overlap[index],
      statistics$background_no_overlap[index]
    )
  )
}))
fwrite(
  contingency_tables,
  file.path(output_directory, "TE_2x2_contingency_tables.tsv"),
  sep = "\t",
  quote = FALSE
)

finite_upper_bounds <- statistics$ci_upper[is.finite(statistics$ci_upper)]
if (length(finite_upper_bounds) == 0L) {
  stop("All confidence interval upper bounds are non-finite; cannot draw plot.")
}
y_max <- max(c(1, finite_upper_bounds)) * 1.15
statistics[, p_label := paste0("P = ", format.pval(p_value, digits = 2, eps = 1e-300))]
statistics[, label_y := pmax(0.15 * y_max, ci_lower - 0.08 * y_max)]

plot_colors <- c(
  "Any" = "#F2C94C",
  "50%" = "#F19F00",
  "100%" = "#E36B00"
)

enrichment_plot <- ggplot(
  statistics,
  aes(x = threshold, y = odds_ratio, color = threshold)
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.8
  ) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.15,
    linewidth = 1
  ) +
  geom_point(size = 4) +
  geom_text(
    aes(y = label_y, label = p_label),
    color = "black",
    size = 4
  ) +
  scale_color_manual(values = plot_colors) +
  scale_y_continuous(
    limits = c(0, y_max),
    breaks = pretty(c(0, y_max))
  ) +
  labs(
    x = "Minimum fraction of TE overlapping\na regulatory element",
    y = "Enrichment ratio\n(ruTEs / background intergenic TEs)",
    title = "Enrichment of ruTEs\nin elite GeneHancer regulatory elements"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 1),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    legend.position = "none"
  )

ggsave(
  filename = file.path(output_directory, "ruTE_enrichment_forest_plot.pdf"),
  plot = enrichment_plot,
  width = 5,
  height = 4.5,
  units = "in",
  device = "pdf",
  useDingbats = FALSE
)

message("Completed GeneHancer enrichment analysis.")
message("  TE-regulatory overlaps: ", file.path(output_directory, "TE_regulatory_overlap_all.tsv"))
message("  Statistics: ", file.path(output_directory, "TE_enrichment_statistics.tsv"))
message("  Contingency tables: ", file.path(output_directory, "TE_2x2_contingency_tables.tsv"))
message("  Figure: ", file.path(output_directory, "ruTE_enrichment_forest_plot.pdf"))
