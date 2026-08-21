#!/usr/bin/env Rscript

# Step 3: create a validated TE-gene map with ruTE/background labels.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
project_root <- if (length(script_arg) > 0L) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."))
} else {
  normalizePath(getwd())
}
source(file.path(project_root, "R", "project_setup.R"))

assert_packages(c("data.table"))
suppressPackageStartupMessages(library(data.table))
assert_files(c(NEAREST_GENE_FILE, RUTE_LIST_FILE))
ensure_parent_dirs(TE_GENE_MAP_FILE)

message("Reading nearest-gene results: ", NEAREST_GENE_FILE)
nearest <- fread(NEAREST_GENE_FILE, na.strings = c("NA", ""))
required_columns <- c("te_id", "gene_id")
missing_columns <- setdiff(required_columns, names(nearest))
if (length(missing_columns) > 0L) {
  stop(
    "Nearest-gene file is missing columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

rute_ids <- trimws(gsub('"', "", fread(RUTE_LIST_FILE, header = FALSE)[[1L]]))

te_gene_map <- unique(nearest[, .(
  te_id,
  closest_gene_in_50kb = gene_id
)])
te_gene_map[, te_group := ifelse(te_id %chin% rute_ids, "ruTE", "background")]

missing_rute <- setdiff(rute_ids, te_gene_map$te_id)
if (length(missing_rute) > 0L) {
  warning(
    length(missing_rute),
    " ruTE loci are absent from the nearest-gene results; adding them with NA genes.",
    call. = FALSE
  )
  te_gene_map <- rbind(
    te_gene_map,
    data.table(
      te_id = missing_rute,
      closest_gene_in_50kb = NA_character_,
      te_group = "ruTE"
    )
  )
}

setorder(te_gene_map, te_group, te_id)
fwrite(te_gene_map, TE_GENE_MAP_FILE, sep = "\t", quote = FALSE, na = "NA")

message("Wrote TE-gene map: ", TE_GENE_MAP_FILE)
message("  total loci: ", nrow(te_gene_map))
message("  ruTE loci: ", sum(te_gene_map$te_group == "ruTE"))
message(
  "  ruTE loci with a gene within ",
  MAX_DISTANCE,
  " bp: ",
  sum(te_gene_map$te_group == "ruTE" & !is.na(te_gene_map$closest_gene_in_50kb))
)
