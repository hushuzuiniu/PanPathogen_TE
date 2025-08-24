library(data.table)
library(readxl)
library(openxlsx)
library(ggplot2)
library(dplyr)

setwd("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_subfamily_feature/plots/")
P1 <- read.csv("All_species_split_DE-TEs_feature_padj0.05log2FC1.csv", header = TRUE)
# P1 <- P1[!grepl("exon", P1$GeneID), ]
TE_Loci_info <- data.frame(transcript_id = P1$GeneID,stringsAsFactors = FALSE)

TE_Loci_info$region <- sapply(strsplit(as.character(P1$GeneID), "_"), function(x) {
  x[length(x)]})

P1_long <- reshape2::melt(P1, id.vars = "GeneID", 
                          variable.name = "Species", 
                          value.name = "Direction")

P1_filtered <- P1_long %>%
  filter(Direction %in% c("Up", "Down"))

merged_data <- merge(P1_filtered, TE_Loci_info, by.x = "GeneID", by.y = "transcript_id")

merged_data$region <- factor(merged_data$region, levels = c("exon", "intron", "intergenic"))

all_species <- unique(P1_long$Species)
all_directions <- c("Up", "Down")
all_regions <- c("exon","intron", "intergenic") 

complete_combinations <- expand.grid(
  Species = all_species,
  Direction = all_directions,
  region = all_regions,
  stringsAsFactors = FALSE
)

region_count_data <- merged_data %>%
  group_by(Species, Direction, region) %>%
  summarise(Count = n(), .groups = "drop")

region_count_data_complete <- merge(complete_combinations, region_count_data, 
                                    by = c("Species", "Direction", "region"), 
                                    all.x = TRUE)

region_count_data_complete$Count[is.na(region_count_data_complete$Count)] <- 0

region_count_data_complete$Count <- ifelse(region_count_data_complete$Direction == "Down", 
                                           -region_count_data_complete$Count, 
                                           region_count_data_complete$Count)

region_count_data_complete$fill_group <- paste(region_count_data_complete$Direction, 
                                               region_count_data_complete$region, sep="_")

region_count_data_complete$fill_group <- factor(region_count_data_complete$fill_group, 
                                                levels = c("Up_exon", "Up_intron", "Up_intergenic",
                                                           "Down_exon", "Down_intron", "Down_intergenic"))
fill_colors <- c(
  "Up_exon" = "#B82132",      
  "Up_intron" = "#D2665A",    
  "Up_intergenic" = "#F2B28C", 
  "Down_exon" = "#2D336B",    
  "Down_intron" = "#7886C7",  
  "Down_intergenic" = "#A9B5DF" 
)

p1 <- ggplot(region_count_data_complete, aes(x = Species, y = Count, fill = fill_group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = fill_colors,
                    name = "",
                    # labels = c("Up_exon", "Up_intron", "Up_intergenic",
                    #            "Down_exon", "Down_intron", "Down_intergenic")) +
                    labels = c("Up_exon","Up_intron", "Up_intergenic",
                              "Down_exon","Down_intron", "Down_intergenic")) +
  scale_y_continuous(
    breaks = seq(-1000, 1000, 200),
    labels = function(x) abs(x), 
    limits = c(-1000, 1000)  
    
  ) +
  theme_minimal() +
  labs(
    title = "Number of differentially expressed TEs subfamily by genomic feature(padj0.05log2FC1)",
    x = "Species",
    y = "No. of TE subfamily feature"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "right",
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.05, "cm")
  ) +
  geom_hline(yintercept = 0, color = "black", size = 0.5)
p1
label_data <- region_count_data_complete %>% 
  group_by(Species, Direction) %>% 
  summarise(
    Total = sum(abs(Count)),
    y_position = ifelse(Direction[1] == "Up", 
                        sum(Count) + 15, 
                        sum(Count) - 15), 
    .groups = "drop"
  )


p2<-p1 + geom_text(
  data = label_data,
  mapping = aes(x = Species, y = y_position, label = Total),
  size = 2, 
  # fontface = "bold",
  inherit.aes = FALSE
)
p2
extract_cell_type <- function(col_name) {
  if (startsWith(col_name, "B_cell")) {
    return("B_cell")
  } else if (startsWith(col_name, "T_cell")) {
    return("T_cell")
  } else {
    return(sub("_.*", "", col_name))
  }
}

celltype_data <- data.frame(
  Species = unique(region_count_data_complete$Species),
  celltype = sapply(as.character(unique(region_count_data_complete$Species)), extract_cell_type),
  stringsAsFactors = FALSE 
)
all_species <- unique(region_count_data_complete$Species)
celltype_data$position <- match(celltype_data$Species, all_species)

celltype_positions <- celltype_data %>%
  group_by(celltype) %>%
  summarise(
    start = min(position),
    end = max(position),
    center = (start + end) / 2,
    .groups = "drop"
  )

p3 <- p2 +
  geom_vline(
    data = celltype_positions,
    aes(xintercept = end + 0.5),
    linetype = "dashed",
    color = "gray60",
    size = 0.5
  ) +
  annotate(
    "text",
    x = celltype_positions$center,
    y = max(region_count_data$Count) * 1.45, 
    label = celltype_positions$celltype,
    size = 2.5,
    fontface = "bold"
  )

ggsave("Fig2_all_TE_subfamily_feature_expression_dataset_by_feature_padj0.05log2FC1_v2.pdf", p3, width = 10, height = 6)
# ggsave("Fig2_no_exon_TE_subfamily_feature_expression_dataset_by_feature_padj0.05log2FC1_v2.pdf", p3, width = 10, height = 6)
# ggsave("Fig2_no_exon_TE_subfamily_feature_expression_dataset_by_feature_padj0.05log2FC1.pdf", p3, width = 10, height = 6)

write.csv(region_count_data_complete,file = "Fig2_TE_subfamily_expression_dataset_by_feature_region_count_data_complete.csv",quote = F,row.names = F)
