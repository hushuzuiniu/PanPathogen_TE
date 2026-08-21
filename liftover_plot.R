#!/usr/bin/env Rscript

# Plot ruTE conservation across 23 representative mammalian species.
#
# Usage:
#   Rscript liftover_plot.R [input_directory] [output_directory]
#
# If omitted, input_directory defaults to the script directory and
# output_directory defaults to <script_directory>/output.

required_packages <- c("ComplexHeatmap", "data.table", "openxlsx")
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
  library(ComplexHeatmap)
  library(data.table)
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

assert_files_exist <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop(
      "Missing input file(s):\n- ",
      paste(missing, collapse = "\n- "),
      call. = FALSE
    )
  }
}

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) > 2L) {
  stop(
    "Usage: Rscript liftover_plot.R [input_directory] [output_directory]",
    call. = FALSE
  )
}

script_directory <- get_script_directory()
input_directory <- if (length(arguments) >= 1L) {
  normalizePath(arguments[[1L]], mustWork = TRUE)
} else {
  script_directory
}
output_directory <- if (length(arguments) >= 2L) {
  arguments[[2L]]
} else {
  file.path(script_directory, "output")
}
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
output_directory <- normalizePath(output_directory, mustWork = TRUE)

presence_matrix_file <- file.path(
  input_directory,
  "input1.ruTE_presence_matrix_liftover.csv"
)
species_annotation_file <- file.path(
  input_directory,
  "input2.mammals_uniq_n95_name_list_NEW.csv"
)
te_annotation_file <- file.path(
  input_directory,
  "input3.rmsk_TE_CopyNumber_Subfamily_Family_Class_n1073_v20240110.xlsx"
)
selected_species_file <- file.path(
  input_directory,
  "input4.selected_23_mammal_species_ordered_by_tree.txt"
)

selected_annotation_output <- file.path(
  output_directory,
  "selected_23_mammal_species_annotations.csv"
)
selected_matrix_output <- file.path(
  output_directory,
  "selected_23_mammal_liftover_matrix.csv"
)
figure_output <- file.path(
  output_directory,
  "Figure3d.ruTE_species_conservation.pdf"
)

assert_files_exist(c(
  presence_matrix_file,
  species_annotation_file,
  te_annotation_file,
  selected_species_file
))

message("Reading the ruTE presence matrix...")
presence_table <- fread(
  presence_matrix_file,
  encoding = "UTF-8",
  check.names = FALSE
)
if (ncol(presence_table) < 2L) {
  stop("The presence matrix must contain a ruTE ID column and species columns.")
}

rute_ids <- presence_table[[1L]]
presence_matrix <- as.matrix(presence_table[, -1L])
storage.mode(presence_matrix) <- "numeric"
rownames(presence_matrix) <- rute_ids
species_by_rute <- t(presence_matrix)

# The human reference loci are present by definition.
human_row <- rep(1, ncol(species_by_rute))
names(human_row) <- colnames(species_by_rute)
species_by_rute <- rbind(Hg38 = human_row, species_by_rute)

message("Selecting representative mammalian species...")
selected_latin_names <- fread(
  selected_species_file,
  header = FALSE,
  encoding = "UTF-8"
)[[1L]]
species_annotation <- fread(species_annotation_file, encoding = "UTF-8")

required_species_columns <- c(
  "SpeciesAbrreviation",
  "LatinName",
  "CommonName"
)
missing_species_columns <- setdiff(
  required_species_columns,
  names(species_annotation)
)
if (length(missing_species_columns) > 0L) {
  stop(
    "The species annotation file is missing columns: ",
    paste(missing_species_columns, collapse = ", "),
    call. = FALSE
  )
}

annotation_indices <- match(
  selected_latin_names,
  species_annotation$LatinName
)
missing_latin_names <- selected_latin_names[is.na(annotation_indices)]
if (length(missing_latin_names) > 0L) {
  stop(
    "Selected species are missing from the annotation table: ",
    paste(missing_latin_names, collapse = ", "),
    call. = FALSE
  )
}
selected_annotation <- species_annotation[annotation_indices]

species_order <- selected_annotation$SpeciesAbrreviation
missing_matrix_species <- setdiff(species_order, rownames(species_by_rute))
if (length(missing_matrix_species) > 0L) {
  warning(
    "Species absent from the presence matrix will be omitted: ",
    paste(missing_matrix_species, collapse = ", "),
    call. = FALSE
  )
}
species_order <- species_order[species_order %in% rownames(species_by_rute)]
if (length(species_order) == 0L) {
  stop("None of the selected species are present in the liftover matrix.")
}
selected_matrix <- species_by_rute[species_order, , drop = FALSE]

fwrite(
  selected_annotation,
  selected_annotation_output,
  bom = TRUE
)
fwrite(
  data.table(
    SpeciesAbrreviation = rownames(selected_matrix),
    as.data.table(selected_matrix)
  ),
  selected_matrix_output,
  bom = TRUE
)

message("Preparing TE class and family annotations...")
te_annotation <- openxlsx::read.xlsx(te_annotation_file)
required_te_columns <- c("TE_Subfamily", "TE_Class", "TE_Family")
missing_te_columns <- setdiff(required_te_columns, names(te_annotation))
if (length(missing_te_columns) > 0L) {
  stop(
    "The TE annotation file is missing columns: ",
    paste(missing_te_columns, collapse = ", "),
    call. = FALSE
  )
}

column_annotation <- data.table(
  ruTE_locus = colnames(selected_matrix)
)
column_annotation[, TE_Subfamily := sub("_c.*$", "", ruTE_locus)]
te_indices <- match(
  column_annotation$TE_Subfamily,
  te_annotation$TE_Subfamily
)
column_annotation[, TE_Class := te_annotation$TE_Class[te_indices]]
column_annotation[, TE_Family := te_annotation$TE_Family[te_indices]]

missing_te_annotations <- column_annotation[
  is.na(TE_Class) | is.na(TE_Family),
  unique(TE_Subfamily)
]
if (length(missing_te_annotations) > 0L) {
  warning(
    length(missing_te_annotations),
    " TE subfamilies have incomplete class or family annotations.",
    call. = FALSE
  )
}

class_colors <- c(
  "DNA" = "#F8D196",
  "LINE" = "#D7A9CB",
  "LTR" = "#6FA4AF",
  "SINE" = "#8290BB",
  "Retroposon" = "grey50"
)
family_colors <- c(
  "Alu" = "#619C60",
  "CR1" = "#FFD58F",
  "ERV1" = "#A2D2E7",
  "ERVK" = "#BF360C",
  "ERVL" = "#DC4AA8",
  "ERVL-MaLR" = "#FF9D9F",
  "Gypsy" = "#E99B78",
  "hAT" = "#DBA9A8",
  "hAT-Blackjack" = "#4D6D7A",
  "hAT-Charlie" = "#CA9600",
  "hAT-Tip100" = "#9A7FBD",
  "L1" = "#FF8831",
  "L2" = "#7CAA8C",
  "MIR" = "#706D54",
  "RTE-X" = "#4DD0E1",
  "SVA" = "#A9B8E0",
  "TcMar-Mariner" = "#F7F2C6",
  "TcMar-Tigger" = "#F4C2C6",
  "tRNA" = "#3B6C6B"
)

unmapped_classes <- setdiff(
  unique(na.omit(column_annotation$TE_Class)),
  names(class_colors)
)
unmapped_families <- setdiff(
  unique(na.omit(column_annotation$TE_Family)),
  names(family_colors)
)
if (length(unmapped_classes) > 0L || length(unmapped_families) > 0L) {
  stop(
    "Missing colors for TE classes/families: ",
    paste(c(unmapped_classes, unmapped_families), collapse = ", "),
    call. = FALSE
  )
}

top_annotation <- HeatmapAnnotation(
  TE_Class = column_annotation$TE_Class,
  TE_Family = column_annotation$TE_Family,
  col = list(
    TE_Class = class_colors,
    TE_Family = family_colors
  ),
  show_legend = TRUE
)

heatmap <- Heatmap(
  selected_matrix,
  name = "Presence",
  col = c("0" = "white", "1" = "#2C6FA5"),
  na_col = "grey90",
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  top_annotation = top_annotation,
  row_names_side = "right",
  border = TRUE,
  heatmap_legend_param = list(
    title = "Presence",
    at = c(0, 1),
    labels = c("Absence", "Presence"),
    color_bar = "discrete"
  )
)

message("Writing heatmap: ", figure_output)
grDevices::pdf(
  figure_output,
  width = 10,
  height = 6,
  useDingbats = FALSE
)
tryCatch(
  ComplexHeatmap::draw(heatmap),
  finally = grDevices::dev.off()
)

message("Completed successfully.")
message("  Species annotations: ", selected_annotation_output)
message("  Liftover matrix: ", selected_matrix_output)
message("  Figure: ", figure_output)
