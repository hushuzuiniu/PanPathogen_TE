#!/usr/bin/env Rscript

# Run Gene Ontology Biological Process enrichment for genes near ruTE loci.
#
# Usage:
#   Rscript GO_BP_enrichment.R [te_gene_table] [output_directory]

required_packages <- c(
  "clusterProfiler",
  "data.table",
  "ggplot2",
  "org.Hs.eg.db"
)
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
  library(clusterProfiler)
  library(data.table)
  library(ggplot2)
  library(org.Hs.eg.db)
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
if (length(arguments) > 2L) {
  stop(
    "Usage: Rscript GO_BP_enrichment.R [te_gene_table] [output_directory]",
    call. = FALSE
  )
}

script_directory <- get_script_directory()
input_file <- if (length(arguments) >= 1L) {
  normalizePath(arguments[[1L]], mustWork = TRUE)
} else {
  file.path(script_directory, "TEs_with_nearest_gene_with_n838_ruTEs.txt")
}
output_directory <- if (length(arguments) >= 2L) {
  arguments[[2L]]
} else {
  file.path(script_directory, "output")
}

if (!file.exists(input_file)) {
  stop("TE-gene table not found: ", input_file, call. = FALSE)
}
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
output_directory <- normalizePath(output_directory, mustWork = TRUE)

message("Reading TE-gene associations: ", input_file)
te_gene_table <- fread(input_file, na.strings = c("NA", ""))
required_columns <- c("te_group", "has_gene_within_50kb", "gene_id")
missing_columns <- setdiff(required_columns, names(te_gene_table))
if (length(missing_columns) > 0L) {
  stop(
    "The TE-gene table is missing columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

within_50kb <- if (is.logical(te_gene_table$has_gene_within_50kb)) {
  te_gene_table$has_gene_within_50kb
} else {
  toupper(trimws(as.character(te_gene_table$has_gene_within_50kb))) %in%
    c("TRUE", "T", "1")
}
te_gene_table[, te_group := trimws(te_group)]

foreground_genes <- unique(te_gene_table[
  te_group == "ruTE" & within_50kb & !is.na(gene_id),
  sub("\\..*$", "", gene_id)
])
background_genes <- unique(te_gene_table[
  te_group == "background" & within_50kb & !is.na(gene_id),
  sub("\\..*$", "", gene_id)
])

if (length(foreground_genes) == 0L || length(background_genes) == 0L) {
  stop("Foreground and background gene sets must both be non-empty.")
}

fwrite(
  data.table(ENSEMBL = foreground_genes),
  file.path(output_directory, "rute_nearest_genes.tsv"),
  sep = "\t"
)
fwrite(
  data.table(ENSEMBL = background_genes),
  file.path(output_directory, "background_te_nearest_genes.tsv"),
  sep = "\t"
)

map_to_entrez <- function(ensembl_ids, label) {
  mapped <- bitr(
    ensembl_ids,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  entrez_ids <- unique(mapped$ENTREZID)
  message(
    "  ",
    label,
    ": ",
    length(entrez_ids),
    "/",
    length(ensembl_ids),
    " genes mapped to Entrez IDs"
  )
  entrez_ids
}

message("Mapping Ensembl IDs to Entrez IDs...")
foreground_entrez <- map_to_entrez(foreground_genes, "ruTE foreground")
background_entrez <- map_to_entrez(background_genes, "TE background")
if (length(foreground_entrez) == 0L || length(background_entrez) == 0L) {
  stop("Entrez mapping produced an empty foreground or background gene set.")
}

message("Running GO Biological Process enrichment...")
go_bp <- enrichGO(
  gene = foreground_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

go_bp_table <- as.data.table(go_bp@result)
if (nrow(go_bp_table) == 0L) {
  stop("GO BP enrichment returned no terms.")
}
fwrite(
  go_bp_table,
  file.path(output_directory, "GO_BP_enrichment.tsv"),
  sep = "\t",
  quote = FALSE
)

plot_data <- go_bp_table[
  is.finite(p.adjust) & p.adjust <= 0.05
][order(p.adjust, pvalue)]

if (nrow(plot_data) == 0L) {
  message("No significant GO BP terms were available for plotting.")
} else {
  plot_data <- head(plot_data, 20L)
  plot_data[, enrichment_score := -log10(pmax(p.adjust, .Machine$double.xmin))]
  plot_data <- plot_data[order(enrichment_score)]
  plot_data[, Description := factor(Description, levels = Description)]
  
  score_max <- max(plot_data$enrichment_score, na.rm = TRUE)
  count_max <- max(plot_data$Count, na.rm = TRUE)
  score_limit <- max(pretty(c(0, score_max), n = 5))
  count_breaks <- pretty(c(0, count_max), n = 5)
  count_limit <- max(count_breaks)
  
  go_bp_plot <- ggplot(plot_data, aes(y = Description)) +
    geom_col(
      aes(x = enrichment_score),
      fill = "#7D578C",
      alpha = 0.65,
      width = 0.8
    ) +
    geom_path(
      aes(x = Count * score_limit / count_limit, group = 1),
      color = "black",
      linewidth = 0.5
    ) +
    geom_point(
      aes(x = Count * score_limit / count_limit),
      color = "black",
      size = 1.2
    ) +
    scale_x_continuous(
      name = expression(-log[10](adjusted~italic(P))),
      position = "top",
      limits = c(0, score_limit),
      expand = expansion(mult = c(0, 0.05)),
      sec.axis = sec_axis(
        transform = ~ . * count_limit / score_limit,
        name = "Gene count",
        breaks = count_breaks
      )
    ) +
    labs(
      title = "GO biological process enrichment",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black"),
      axis.title.x = element_text(color = "black"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5)
    )
  
  ggsave(
    filename = file.path(output_directory, "GO_BP_enrichment.pdf"),
    plot = go_bp_plot,
    width = 6,
    height = 4,
    units = "in",
    device = "pdf",
    useDingbats = FALSE
  )
}

message("Completed GO Biological Process enrichment analysis.")
message("  Output directory: ", output_directory)
