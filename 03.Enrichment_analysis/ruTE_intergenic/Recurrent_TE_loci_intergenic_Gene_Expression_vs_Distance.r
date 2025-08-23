library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
setwd("/data2t_2/pathogen_TE_2025_New/03.Enrichment_analysis/")
matrix <- read.csv("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_Gene/plots/All_species_split_Gene_padj0.05log2FC1.csv", header = T, row.names = NULL)
colnames(matrix)[colnames(matrix) == "GeneID"] <- "Subfamily"
data <- tidyr::pivot_longer(matrix,
                            cols = -Subfamily,
                            names_to = "Species",
                            values_to = "Value")
############################################################
data <- data %>%
  mutate(
    Pathogen = case_when(
      str_detect(Species, "^Huh-7_") ~ str_replace(Species, "^Huh-7_", ""),
      str_detect(Species, "^T_cell_") ~ str_replace(Species, "^T_cell_", ""),
      str_detect(Species, "^B_cell_") ~ str_replace(Species, "^B_cell_", ""),
      TRUE ~ str_replace(Species, "^[^_]+_", "")
    )
  )

Up_pathogen_summary <- data %>%
  filter(Value == "Up") %>%
  group_by(Subfamily) %>%
  summarise(
    SumPathogens = n(),               
    UniqPathogens = n_distinct(Pathogen)  
  ) %>%
  filter(SumPathogens >= 1) %>%
  arrange(desc(SumPathogens), desc(UniqPathogens))
Down_pathogen_summary <- data %>%
  filter(Value == "Down") %>%
  group_by(Subfamily) %>%
  summarise(
    SumPathogens = n(),               
    UniqPathogens = n_distinct(Pathogen) 
  ) %>%
  filter(SumPathogens >= 1) %>%
  arrange(desc(SumPathogens), desc(UniqPathogens))
#############################################
all_genes <- full_join(
  Up_pathogen_summary %>% dplyr::select(Subfamily) %>% mutate(Up = TRUE),
  Down_pathogen_summary %>% dplyr::select(Subfamily) %>% mutate(Down = TRUE),
  by = "Subfamily"
) %>%
  mutate(
    Up = ifelse(is.na(Up), FALSE, Up),
    Down = ifelse(is.na(Down), FALSE, Down),
    Expression_Type = case_when(
      Up & Down ~ "Mixed",
      Up ~ "Up",
      Down ~ "Down",
      TRUE ~ "/"  
    )
  )

TE_closest_genes<-read.csv("Recurrent_TE_Loci_gt5_intergenic_n838_closest_genes.bed",sep = "\t",header = F)
colnames(TE_closest_genes) <- c(
  "chr1", "start", "end", "TE_id", "UniqPathogens","TE_location" ,"class", "family","feature",
  "chr2", "gene_start", "gene_end", "gene_width", "strand",  "gene_info","distance"
)
# TE_closest_genes<-TE_closest_genes[TE_closest_genes$feature_type=="intergenic",]

extract_gene_id <- function(gene_info) {
  gene_id <- str_extract(gene_info, "gene_id\\s+([^;]+)")
  gene_id <- str_replace(gene_id, "gene_id\\s+", "")
  return(gene_id)
}

te_gene_data <- TE_closest_genes %>%
  mutate(gene_id = sapply(gene_info, extract_gene_id))

te_gene_data <- te_gene_data %>%
  left_join(all_genes, by = c("gene_info" = "Subfamily")) %>%
  mutate(relationship = paste0("UpTE-", ifelse(is.na(Expression_Type), "/", Expression_Type), "Gene"))

custom_breaks <- c(0,1000, 2000, 3000,4000, 5000, 6000,7000,8000,9000,10000, 20000,30000, 60000, Inf)
custom_labels <- c("<=1kb", "1kb-2kb", "2kb-3kb", "3kb-4kb",
                   "4kb-5kb", "5kb-6kb","6kb-7kb","7kb-8kb","8kb-9kb", "9kb-10kb",
                   "10kb-20kb","20kb-30kb", "30kb-60kb", ">60kb")

te_gene_data <- te_gene_data %>%
  mutate(
    distance_range = cut(
      distance,
      breaks = custom_breaks,
      labels = custom_labels,
      include.lowest = TRUE,
      right = FALSE
    )
  )

relationship_counts <- te_gene_data %>%
  filter(relationship %in% c("UpTE-UpGene", "UpTE-DownGene", "UpTE-MixedGene", "UpTE-not sigGene")) %>%
  group_by(distance_range, relationship) %>%
  summarise(count = n(), .groups = "drop")

all_combinations <- expand.grid(
  distance_range = custom_labels,
  relationship = c("UpTE-UpGene", "UpTE-DownGene", "UpTE-MixedGene", "UpTE-not sigGene")
)

relationship_counts <- all_combinations %>%
  left_join(relationship_counts, by = c("distance_range", "relationship")) %>%
  mutate(count = ifelse(is.na(count), 0, count))

relationship_counts$distance_range <- factor(
  relationship_counts$distance_range,
  levels = custom_labels
)

relationship_counts$relationship <- factor(relationship_counts$relationship, 
                                           levels = c("UpTE-UpGene", "UpTE-DownGene", "UpTE-MixedGene", "UpTE-not sigGene"))
p <- ggplot(relationship_counts, aes(x = distance_range, y = count, group = relationship, color = relationship)) +
  geom_line(size = 0.5) +
  geom_point(size = 1.5) +  
  labs(
    title = "intergenic Recurrent TE_loci-Gene Expression Relationships with Distance",
    x = "absolute distance to TE instance(bp)",
    y = "Count",
    color = ""
  ) +
  theme_minimal() +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    # legend.position = "right",
    legend.position = c(0.7, 0.6),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    axis.title = element_text(size = 7,face = "bold"),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white"),
    axis.line = element_line(colour = "black"),  
    axis.ticks = element_line(colour = "black"),  
    axis.ticks.length = unit(0.1, "cm")  
  ) +
  scale_color_manual(
    values = c("UpTE-UpGene" = "#88304E", "UpTE-DownGene" = "#003161", "UpTE-MixedGene" = "#78A083", "UpTE-not sigGene" = "#DCD7C9")
  ) +
  scale_y_continuous(
    breaks = seq(0, max(relationship_counts$count, na.rm = TRUE) + 100, by = 100),  
    expand = c(0.02, 0)  
  ) +
  scale_y_continuous(
    breaks = seq(0, max(relationship_counts$count, na.rm = TRUE) + 1, by = 20),  
    expand = c(0.02, 0)  
  ) +
  scale_x_discrete(
    expand = c(0.01, 0.05)  
  ) +
  coord_cartesian(clip = "off") 

print(p)
ggsave("/data2t_2/pathogen_TE_2025_New/03.Enrichment_analysis/Fig3_Recurrent_TE_Loci_gt5_intergenic_Gene_Expression_vs_Distance.pdf", p, width = 5, height = 3)