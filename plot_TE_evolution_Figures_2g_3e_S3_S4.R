#!/usr/bin/env Rscript

# Reproduce the TE evolutionary-age analyses shown in Figures 3e, 5E, and 6.
#
# By default, input files are read relative to this script and PDFs are written
# to the same directory. The locations can be overridden when needed:
#
#   PANPATHOGEN_TE_DATA_DIR=/path/to/data \
#   PANPATHOGEN_TE_OUTPUT_DIR=/path/to/output \
#   Rscript plot_TE_evolution_Figures_3e_5E_6.R
#
# Required packages: tidyverse, readxl, and ggrepel.

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ggrepel)
})

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()
data_dir <- Sys.getenv("PANPATHOGEN_TE_DATA_DIR", unset = script_dir)
output_dir <- Sys.getenv("PANPATHOGEN_TE_OUTPUT_DIR", unset = script_dir)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- function(...) file.path(data_dir, ...)
output_path <- function(filename) file.path(output_dir, filename)

required_files <- c(
  input_path("ruTE_loci_n838_with_milliDiv_and_MYA.bed"),
  input_path("input3.rmsk_TE_CopyNumber_Subfamily_Family_Class_n1073_v20240110.xlsx"),
  input_path("all_activated_TE_loci", "TE_subfamily_mean_age.tsv"),
  input_path(
    "all_activated_TE_loci",
    "Recurrent_Up_TE_Loci_intergenic_intron_exon_UniqPathogens_gt1_subfamily_summary.csv"
  ),
  input_path("upTE_loci_n295979_with_milliDiv_and_MYA.bed")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Required input file(s) not found:\n", paste0("  - ", missing_files, collapse = "\n"),
    "\nSet PANPATHOGEN_TE_DATA_DIR to the directory containing the input data."
  )
}

# =============================================================================
# Shared ruTE data
# =============================================================================

rute <- read_tsv(
  input_path("ruTE_loci_n838_with_milliDiv_and_MYA.bed"),
  col_names = FALSE,
  show_col_types = FALSE
)

if (ncol(rute) != 19L) {
  stop("The ruTE BED file must contain 19 columns; found ", ncol(rute), ".")
}

colnames(rute) <- c(
  "Chr", "Start", "End", "Locus_ID", "TE_strand", "UniqPatho_Freq",
  "Locus_Name", "TE_Class", "TE_Family", "Genomic_Region", "Gene_Chr",
  "Gene_Start", "Gene_End", "Gene_Length", "Gene_Strand", "Gene_ID",
  "TE_Gene_Distance", "milliDiv", "Age_Mya"
)

rute <- rute %>%
  mutate(Age_Mya = as.numeric(Age_Mya)) %>%
  filter(!is.na(Age_Mya))

# =============================================================================
# Figure3e: ruTE age distribution with a reversed x-axis
# Output: Figure3e_ruTE_age_distribution_reversed_xaxis.pdf
# =============================================================================

median_age <- median(rute$Age_Mya, na.rm = TRUE)
mean_age <- mean(rute$Age_Mya, na.rm = TRUE)
xmax <- ceiling(max(rute$Age_Mya, na.rm = TRUE) / 25) * 25

p3e <- ggplot(rute, aes(x = Age_Mya)) +
  geom_histogram(
    aes(y = after_stat(density)),
    binwidth = 5,
    boundary = 0,
    fill = "#9EC5C4",
    color = "white",
    alpha = 0.9
  ) +
  geom_density(color = "#1F4E5F", linewidth = 1.2) +
  geom_vline(
    xintercept = median_age,
    linetype = "dashed",
    color = "#C0392B",
    linewidth = 0.9
  ) +
  geom_vline(
    xintercept = mean_age,
    linetype = "dotted",
    color = "#7F8C8D",
    linewidth = 0.9
  ) +
  annotate(
    "text",
    x = median_age,
    y = Inf,
    label = paste0("Median = ", round(median_age, 1), " MYA"),
    vjust = 1.5,
    hjust = 1.1,
    color = "#C0392B",
    size = 4
  ) +
  annotate(
    "text",
    x = mean_age,
    y = Inf,
    label = paste0("Mean = ", round(mean_age, 1), " MYA"),
    vjust = 3.0,
    hjust = 1.1,
    color = "#7F8C8D",
    size = 4
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "A",
    hjust = 1.3,
    vjust = 1.3,
    size = 6,
    fontface = "bold"
  ) +
  scale_x_reverse(
    limits = c(xmax, 0),
    breaks = seq(0, xmax, by = 25),
    minor_breaks = seq(0, xmax, by = 5),
    expand = c(0, 0)
  ) +
  labs(x = "Evolutionary Age (MYA)", y = "Density") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.ticks.length = grid::unit(0.22, "cm"),
    axis.ticks = element_line(linewidth = 0.5),
    plot.margin = margin(10, 15, 10, 15)
  )

print(p3e)
ggsave(
  output_path("Figure3e_ruTE_age_distribution_reversed_xaxis.pdf"),
  p3e,
  width = 8,
  height = 5.5
)

# =============================================================================
# FigureS3: ruTE subfamily age versus activated-copy fraction
# Output: FigureS3_ruTE_age_vs_activated_frequency.pdf
# =============================================================================

rute_by_subfamily <- rute %>%
  mutate(TE_Subfamily = sub("_chr.*", "", Locus_Name))

activated_count <- rute_by_subfamily %>%
  count(TE_Subfamily, name = "ActivatedCopies")

age_info <- rute_by_subfamily %>%
  group_by(TE_Subfamily) %>%
  summarise(Age_Mya = median(Age_Mya, na.rm = TRUE), .groups = "drop")

rmsk <- read_excel(
  input_path("input3.rmsk_TE_CopyNumber_Subfamily_Family_Class_n1073_v20240110.xlsx")
)

figureS3_data <- activated_count %>%
  left_join(rmsk, by = "TE_Subfamily") %>%
  left_join(age_info, by = "TE_Subfamily") %>%
  mutate(ActivationRate = ActivatedCopies / CopyNumber) %>%
  filter(is.finite(Age_Mya), is.finite(ActivationRate))

figureS3_cor <- cor.test(
  figureS3_data$Age_Mya,
  figureS3_data$ActivationRate,
  method = "spearman",
  exact = FALSE
)

print(figureS3_cor)

pS3 <- ggplot(figureS3_data, aes(x = Age_Mya, y = ActivationRate)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_bw(base_size = 14) +
  labs(
    x = "TE Subfamily Age (MYA)",
    y = "Activated Copy Fraction",
    title = "Relationship Between TE Age and Activation Frequency"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.1,
    vjust = 1.5,
    label = paste0(
      "Spearman rho = ", round(figureS3_cor$estimate, 3),
      "\nP = ", signif(figureS3_cor$p.value, 3)
    )
  )

print(pS3)
ggsave(
  output_path("FigureS3_ruTE_age_vs_activated_frequency.pdf"),
  pS3,
  width = 6,
  height = 4
)

# =============================================================================
# Figure2g: all activated TE subfamilies
# Output: Figure2g_upTE_age_vs_activation_percentage_correlation.pdf
# =============================================================================

all_te_age <- read_tsv(
  input_path("all_activated_TE_loci", "TE_subfamily_mean_age.tsv"),
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
) %>%
  # The source file contains one repeated header row; exclude it explicitly.
  filter(TE_Subfamily != "subfamily", Mean_Age_MYA != "mean_age_mya") %>%
  mutate(Mean_Age_MYA = as.numeric(Mean_Age_MYA))

upte_summary <- read_csv(
  input_path(
    "all_activated_TE_loci",
    "Recurrent_Up_TE_Loci_intergenic_intron_exon_UniqPathogens_gt1_subfamily_summary.csv"
  ),
  show_col_types = FALSE
)

figure2g_data <- upte_summary %>%
  left_join(all_te_age, by = c("subfamily" = "TE_Subfamily")) %>%
  transmute(
    subfamily,
    family_id,
    class_id,
    Mean_Age_MYA = as.numeric(Mean_Age_MYA),
    up_pct_in_subfamily = as.numeric(up_pct_in_subfamily)
  ) %>%
  drop_na()

figure2g_cor <- cor.test(
  figure2g_data$Mean_Age_MYA,
  figure2g_data$up_pct_in_subfamily,
  method = "spearman",
  exact = FALSE
)

rho <- round(figure2g_cor$estimate, 4)
p_value <- figure2g_cor$p.value
p_label <- ifelse(
  p_value < 2.2e-16,
  "< 2.2e-16",
  paste0("= ", format(p_value, scientific = TRUE, digits = 3))
)

x_max <- ceiling(max(figure2g_data$Mean_Age_MYA) / 25) * 25
x_breaks <- seq(0, x_max, by = 25)
y_max <- ceiling(max(figure2g_data$up_pct_in_subfamily) / 0.1) * 0.1
y_breaks <- seq(0, y_max, by = 0.1)

top10 <- figure2g_data %>%
  arrange(desc(up_pct_in_subfamily)) %>%
  slice_head(n = 10) %>%
  mutate(label_text = paste0(subfamily, "(", family_id, ")"))

p2g <- ggplot(
  figure2g_data,
  aes(x = Mean_Age_MYA, y = up_pct_in_subfamily, color = class_id)
) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(
    aes(group = 1),
    method = "loess",
    se = TRUE,
    color = "black",
    linetype = "dashed",
    linewidth = 0.8,
    alpha = 0.2
  ) +
  scale_x_continuous(breaks = x_breaks, limits = c(0, x_max)) +
  scale_y_continuous(breaks = y_breaks, limits = c(0, y_max)) +
  scale_color_manual(
    values = c(
      "DNA" = "#66C2A5",
      "LINE" = "#FC8D62",
      "LTR" = "#8DA0CB",
      "SINE" = "#A6D854",
      "Retroposon" = "grey70"
    ),
    name = "TE Class"
  ) +
  geom_text_repel(
    data = top10,
    aes(label = label_text),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey50",
    segment.size = 0.3,
    show.legend = FALSE,
    seed = 123
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = grid::unit(0.15, "cm"),
    legend.position = c(0.88, 0.58),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    legend.key = element_blank()
  ) +
  labs(
    x = "Mean Age of TE Subfamily (MYA)",
    y = "Proportion of Activated TE Loci"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0("Spearman's rho = ", rho, "\np ", p_label),
    hjust = 1.1,
    vjust = 1.5,
    size = 4,
    fontface = "italic"
  )

print(p2g)
ggsave(
  output_path("Figure2g_upTE_age_vs_activation_percentage_correlation.pdf"),
  p2g,
  width = 6,
  height = 4,
  dpi = 300
)

# =============================================================================
# FigureS4: all activated TE loci, age versus unique-pathogen count
# Output: FigureS4_all_upTE_loci_age_vs_unique_pathogen_count.pdf
# =============================================================================

upte_loci <- read_tsv(
  input_path("upTE_loci_n295979_with_milliDiv_and_MYA.bed"),
  col_names = FALSE,
  show_col_types = FALSE
)

if (ncol(upte_loci) < 10L) {
  stop("The upTE BED file must contain at least 10 columns; found ", ncol(upte_loci), ".")
}

colnames(upte_loci)[6] <- "Unique_pathogen_count"
colnames(upte_loci)[10] <- "Region"
colnames(upte_loci)[ncol(upte_loci)] <- "Age_Mya"

upte_loci <- upte_loci %>%
  mutate(
    Unique_pathogen_count = as.numeric(Unique_pathogen_count),
    Age_Mya = as.numeric(Age_Mya),
    LocusType = if_else(
      Region == "intergenic" & Unique_pathogen_count >= 5,
      "ruTE loci",
      "Other TE loci"
    )
  ) %>%
  filter(!is.na(Unique_pathogen_count), !is.na(Age_Mya))

FigureS4_cor <- cor.test(
  upte_loci$Age_Mya,
  upte_loci$Unique_pathogen_count,
  method = "spearman",
  exact = FALSE
)

cor_text <- paste0(
  "Spearman rho = ", round(FigureS4_cor$estimate, 3),
  "\nP = ", signif(FigureS4_cor$p.value, 3)
)

xmax <- ceiling(max(upte_loci$Age_Mya, na.rm = TRUE) / 25) * 25
ymax <- ceiling(max(upte_loci$Unique_pathogen_count, na.rm = TRUE))

set.seed(123)
p6 <- ggplot(upte_loci, aes(x = Age_Mya, y = Unique_pathogen_count)) +
  geom_jitter(
    aes(color = LocusType),
    width = 1.2,
    height = 0.15,
    size = 2,
    alpha = 0.7
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#E69F00",
    fill = "#F6D58B",
    linewidth = 1.2
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "D",
    hjust = -0.4,
    vjust = 1.3,
    size = 6,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = cor_text,
    hjust = 1.05,
    vjust = 1.3,
    size = 4.2
  ) +
  scale_x_continuous(
    breaks = seq(0, xmax, by = 25),
    minor_breaks = seq(0, xmax, by = 5)
  ) +
  scale_y_continuous(
    breaks = seq(0, ymax, by = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_color_manual(
    values = c("Other TE loci" = "#4C78A8", "ruTE loci" = "#D62728"),
    name = NULL
  ) +
  coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax)) +
  labs(
    x = "Evolutionary Age (MYA)",
    y = "Number of Unique Pathogens"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.ticks.length = grid::unit(0.22, "cm"),
    axis.ticks = element_line(linewidth = 0.5),
    plot.margin = margin(10, 15, 10, 15)
  )

print(p6)
ggsave(
  output_path("FigureS4_all_upTE_loci_age_vs_unique_pathogen_count.pdf"),
  p6,
  width = 12,
  height = 6
)

message("Finished. Figure PDFs were written to: ", normalizePath(output_dir))
