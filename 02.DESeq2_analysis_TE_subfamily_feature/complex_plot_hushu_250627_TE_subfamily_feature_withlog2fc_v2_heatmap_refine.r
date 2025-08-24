library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid) 
library(gtable)
library(RColorBrewer)
library(stringr)
library(data.table)
library(patchwork)
setwd("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_subfamily_feature/plots/")
matrix <- read.csv("All_species_split_Up_DE-TEs_feature_padj0.05log2FC1.csv", header = T, row.names = NULL)
colnames(matrix)[colnames(matrix) == "A549_Virus_SARSCoV"] <- "A549_Virus_SARSCoV2"
colnames(matrix)[colnames(matrix) == "Macrophages_Virus_SARSCoV"] <- "Macrophages_Virus_SARSCoV2"
colnames(matrix)[colnames(matrix) == "T_cell_HIV"] <- "T_cell_Virus_HIV1"
logfc<-read.csv("All_celltypes_combined_Average_log2FoldChange_sig_DE-TEs_feature_padj0.05log2FC1.csv")
# origin_data <- read.csv("/home/pathogen_TE_2025_New/pathogen_TE_2025/02.DESeq2_analysis/plots/Origin_data.csv", header = F)
colnames(matrix)[colnames(matrix) == "GeneID"] <- "Subfamily"
colnames(logfc)[colnames(logfc) == "GeneID"] <- "Subfamily"

TE_Loci_info<-fread("/data2t_2/pathogen_TE_2025_New/01.Genomic_features_Gencode/hg38_TE_anno_custom_v20240110_0based_with_feature_summary_v2.bed")
TE_Loci_info <- TE_Loci_info[, c("gene_id", "transcript_id", "class_id", "family_id", "prioritized_feature")]
TE_Loci_info$region <- TE_Loci_info$prioritized_feature
TE_Loci_info$region[TE_Loci_info$prioritized_feature %in%
                      c("3UTR", "5UTR", "CDS", "noncoding_exon")] <- "exon"
# matrix <- matrix[!grepl("_exon$|_intron$", matrix$Subfamily), ]
# matrix <- matrix[!grepl("_exon$", matrix$Subfamily), ]
############################################################

data <- tidyr::pivot_longer(matrix,
                            cols = -Subfamily,
                            names_to = "Species",
                            values_to = "Value")

data_with_pathogens <- data %>%
  filter(Value == 1) %>%
  mutate(
    Pathogen = case_when(
      str_detect(Species, "^Huh-7_") ~ str_replace(Species, "^Huh-7_", ""),
      str_detect(Species, "^T_cell_") ~ str_replace(Species, "^T_cell_", ""),
      str_detect(Species, "^B_cell_") ~ str_replace(Species, "^B_cell_", ""),
      TRUE ~ str_replace(Species, "^[^_]+_", "")
    )
  )

pathogen_counts <- data_with_pathogens %>%
  group_by(Subfamily) %>%
  summarise(UniqPathogens = n_distinct(Pathogen))
# write.csv(pathogen_counts,file = "Recurrent_Up_TEs_feature_intergenic_UniqPathogens_gt1.csv",quote = F,row.names = F)
# write.csv(pathogen_counts,file = "Recurrent_Up_TEs_feature_intergenic_intron_UniqPathogens_gt1.csv",quote = F,row.names = F)
# write.csv(pathogen_counts,file = "Recurrent_Up_TEs_feature_intergenic_intron_exon_UniqPathogens_gt1.csv",quote = F,row.names = F)

# #top100
# sorted_subfamilies <- pathogen_counts %>%
#   arrange(desc(UniqPathogens)) %>%
#   head(100)

sorted_subfamilies <- pathogen_counts %>%
  filter(UniqPathogens >= 5) %>%
  arrange(desc(UniqPathogens))

# write.csv(sorted_subfamilies,"Fig2_All_species_Up_DE-TEs_feature_padj0.05log2FC1_v2_intergenic_top5.csv",quote = F,row.names = F)
# write.csv(sorted_subfamilies,"Fig2s_All_species_Up_DE-TEs_feature_padj0.05log2FC1_v2_intergenic_intron_top100.csv",quote = F,row.names = F)
# write.csv(sorted_subfamilies,"Fig2s_All_species_Up_DE-TEs_feature_padj0.05log2FC1_v2_intergenic_intron_exon_top100.csv",quote = F,row.names = F)
ordered_subfamilies_by_pathogens <- sorted_subfamilies$Subfamily

data_filtered <- data %>%
  filter(Subfamily %in% ordered_subfamilies_by_pathogens)

species_sum <- data_filtered %>% 
  filter(Value == 1) %>%
  group_by(Species) %>%
  summarise(Sum = sum(Value))

logfc_long <- logfc %>%
  tidyr::pivot_longer(cols = -Subfamily, 
                      names_to = "Species", 
                      values_to = "LogFC")

plot_data <- data_filtered %>%
  filter(Value == 1) %>%
  left_join(logfc_long, by = c("Subfamily", "Species")) %>%
  left_join(pathogen_counts, by = "Subfamily")

plot_data$Subfamily <- factor(plot_data$Subfamily, levels = rev(ordered_subfamilies_by_pathogens))

right_bar_df <- data.frame(
  Subfamily = ordered_subfamilies_by_pathogens
) %>%
  left_join(pathogen_counts, by = "Subfamily")

right_bar_df$UniqPathogens[is.na(right_bar_df$UniqPathogens)] <- 0


create_spaced_grouped_heatmap <- function(data) {
  data_grouped <- data %>%
    arrange(UniqPathogens, Subfamily)
  unique_groups <- sort(unique(data_grouped$UniqPathogens))
  
  processed_data <- data.frame()
  group_info <- data.frame()
  current_y_position <- 1
  
  for(i in seq_along(unique_groups)) {
    group_data <- data_grouped %>% 
      filter(UniqPathogens == unique_groups[i])
    unique_subfamilies <- unique(group_data$Subfamily)
    group_start_y <- current_y_position
    
    for(j in seq_along(unique_subfamilies)) {
      subfamily_data <- group_data %>% 
        filter(Subfamily == unique_subfamilies[j])
      
      subfamily_data$y_position <- current_y_position
      subfamily_data$group_id <- i
      subfamily_data$subfamily_for_label <- unique_subfamilies[j]
      
      processed_data <- rbind(processed_data, subfamily_data)
      current_y_position <- current_y_position + 1
    }
    
    group_end_y <- current_y_position - 1
    group_info <- rbind(group_info, data.frame(
      group_id = i,
      UniqPathogens = unique_groups[i],
      y_min = group_start_y - 0.5,
      y_max = group_end_y + 0.5,
      x_min = 0.5,
      x_max = length(unique(data_grouped$Species)) + 0.5
    ))
    if(i < length(unique_groups)) {
      current_y_position <- current_y_position + 1  
    }
  }
  
  return(list(data = processed_data, group_info = group_info))
}

processed_result <- create_spaced_grouped_heatmap(plot_data)
plot_data_spaced <- processed_result$data
group_info <- processed_result$group_info

y_breaks <- unique(plot_data_spaced$y_position)
y_labels <- plot_data_spaced %>% 
  select(y_position, subfamily_for_label) %>% 
  distinct() %>% 
  arrange(y_position) %>% 
  pull(subfamily_for_label)


all_species <- unique(plot_data_spaced$Species)
all_y_positions <- unique(plot_data_spaced$y_position)

grid_data <- expand.grid(
  Species = all_species,
  y_position = all_y_positions,
  stringsAsFactors = FALSE
)

logfc_range <- range(plot_data_spaced$LogFC, na.rm = TRUE)
size_range <- c(0.5, 2.5)  
clean_labels <- function(labels) {
  # cleaned <- gsub("_intergenic$|_intron$", "", labels)
  cleaned <- gsub("", "", labels)
  return(cleaned)
}
heatmap_grouped <- ggplot() +
  geom_tile(data = grid_data, 
            aes(x = Species, y = y_position), 
            fill = "white", color = "#bfbfbf", size = 0.1) +
  
  geom_point(data = plot_data_spaced,
             aes(x = Species, y = y_position, size = LogFC), 
             color = "#0D47A1", alpha = 0.8) +
  
  scale_size_continuous(range = size_range, 
                        name = "log2FC",
                        guide = guide_legend(override.aes = list(alpha = 0.8))) +
  
  geom_rect(data = group_info, 
            aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max),
            fill = NA, color = "black", size = 0.5, inherit.aes = FALSE) +
  
  scale_y_continuous(breaks = y_breaks, labels = clean_labels(y_labels), expand = c(0.02, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(x = "Species", y = "TE_Subfamily")

heatmap_grouped

heatmap_no_legend <- heatmap_grouped + 
  theme(legend.position = "none",
        plot.margin = margin(5, 5, 5, 5))

species_order <- levels(factor(plot_data_spaced$Species))
species_sum <- data.frame(Species = species_order)

species_sum$Color <- ifelse(grepl("A549", species_sum$Species), "#8DD1C6",
                            ifelse(grepl("B_cell", species_sum$Species), "#FEFEB3",
                                   ifelse(grepl("Dendritic", species_sum$Species), "#BBB8D9",
                                          ifelse(grepl("hSAEC", species_sum$Species), "#A8E6CF",
                                                 ifelse(grepl("Huh.7", species_sum$Species), "#F8BBD9",  
                                                        ifelse(grepl("Liver", species_sum$Species), "#FA7F73",
                                                               ifelse(grepl("Macrophages", species_sum$Species), "#7FAFD1",
                                                                      ifelse(grepl("Monocyte", species_sum$Species), "#FCB264",
                                                                             ifelse(grepl("PBMCs", species_sum$Species), "#B1DD6D",
                                                                                    ifelse(grepl("T_cell", species_sum$Species), "#DCD7C9", "gray")))))))))) 



species_sum$Species <- factor(species_sum$Species, levels = species_order)

color_bar <- ggplot(species_sum, aes(x = Species, y = 1, fill = Color)) +
  geom_tile() +
  scale_fill_identity() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  # theme(plot.margin = margin(5, 5, -2, 80))
  theme(plot.margin = margin(5, -8, -2, 110))
#margin(top, right, bottom, left)

y_breaks <- unique(plot_data_spaced$y_position)
y_range <- range(plot_data_spaced$y_position)

create_color_data <- function(loci_info, column_name, colors_map) {
  color_data <- data.frame()
  
  for(y_pos in y_breaks) {
    subfamily <- plot_data_spaced$subfamily_for_label[plot_data_spaced$y_position == y_pos][1]
    subfamily_clean <- gsub("_intron$|_intergenic$|_exon$", "", subfamily)
    type_value <- loci_info[[column_name]][loci_info$gene_id == subfamily_clean]
    if(length(type_value) == 0) type_value <- "unknown"
    
    color_data <- rbind(color_data, data.frame(
      y_position = y_pos,
      type = type_value[1],
      color = colors_map[type_value[1]]
    ))
  }
  
  color_data$color[is.na(color_data$color)] <- "white"
  
  return(color_data)
}

colors <- c("DNA" = "#f8d196", "LINE" = "#d7a9cb", "LTR" = "#6fa4af",
            "SINE" = "#8290bb", "Retroposon" = "#277899")
subfamily_families <- data.frame()
unique_subfamilies <- unique(plot_data_spaced$subfamily_for_label)

for(subfamily in unique_subfamilies) {
  subfamily_clean <- gsub("_intron$|_intergenic$|_exon$", "", subfamily)
  
  family <- TE_Loci_info$family_id[TE_Loci_info$gene_id == subfamily_clean]
  if(length(family) > 0) {
    subfamily_families <- rbind(subfamily_families, data.frame(
      subfamily = subfamily,  
      family_id = family[1]
    ))
  } else {
    subfamily_families <- rbind(subfamily_families, data.frame(
      subfamily = subfamily,  # 
      family_id = "not_found"
    ))
  }
}
table(subfamily_families$family_id)
colors_family <- c("ERV1" = "#a2d2e7", "ERVL" = "#dc4aa8", "L1" = "#ca9600",  
                   "Alu" = "#619c60", "L2" = "#3955a1", "ERVL-MaLR" = "#fff59b", 
                   "hAT-Tip100" = "#f36569", "hAT-Blackjack" = "#4d6d7a", "TcMar-Tigger" = "#9a7fbd", 
                   "hAT-Charlie" = "#ff8831",  "MIR" = "#706D54", "RTE-X"="#c35338",
                   "CR1" = "#FFDA76","hAT"="#dba9a8","Gypsy"="#e99b78","PiggyBac"="#393E46",
                   "TcMar-Tc2"="#c2bb82","ERVK"="#eaa944","TcMar-Mariner" = "#872341","SVA"="#97866A","TcMar"="#A62C2C")
 
# colors_family <- c("ERV1" = "#a2d2e7", "ERVL" = "#dc4aa8", "L1" = "#ca9600", "TcMar-Tc2" = "#c2bb82",
#                    "hAT-Charlie" = "#6fb3a8", "Alu" = "#619c60", "L2" = "#3955a1", "ERVL-MaLR" = "#fff59b",
#                    "hAT-Tip100" = "#f36569", "Gypsy" = "#e99b78", "ERVK" = "#eaa944", "hAT-Blackjack" = "#704ba3",
#                    "TcMar-Tigger" = "#9a7fbd", "hAT" = "#ff8831", "SVA" = "#dba9a8", "MIR" = "#706D54", "RTE-X"="#c35338",
#                    "TcMar-Mariner" = "#872341", "CR1" = "#FFDA76", "MULE-MuDR" = "#393E46", "unknown" = "white")

feature_colors <- c(
  "exon" = "#b3a8c8",
  "intron" = "#cbd7e9",       
  "intergenic" = "#6c9fa9"
)

te_color_data <- create_color_data(TE_Loci_info, "class_id", colors)
family_color_data <- create_color_data(TE_Loci_info, "family_id", colors_family)

right_color_bar_with_legend <- ggplot(te_color_data, aes(x = 1, y = y_position, fill = type)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_manual(values = colors, name = "TE Class") +
  scale_y_continuous(expand = c(0.02, 0), limits = c(y_range[1] - 0.5, y_range[2] + 0.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(
    plot.margin = margin(5, 0, 5, 0),
    legend.position = "right",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6)
  )
right_color_bar2_with_legend <- ggplot(family_color_data, aes(x = 1, y = y_position, fill = type)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_manual(values = colors_family, name = "TE Family") +
  scale_y_continuous(expand = c(0.02, 0), limits = c(y_range[1] - 0.5, y_range[2] + 0.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(
    plot.margin = margin(5, 5, 5, 0),
    legend.position = "right",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6)
  )

create_region_color_data <- function() {
  color_data <- data.frame()
  for(y_pos in y_breaks) {
    subfamily <- plot_data_spaced$subfamily_for_label[plot_data_spaced$y_position == y_pos][1]
    region_type <- case_when(
      str_detect(subfamily, "_intron$") ~ "intron",
      str_detect(subfamily, "_intergenic$") ~ "intergenic", 
      str_detect(subfamily, "_exon$") ~ "exon",
      TRUE ~ "unknown"  
    )
    color_data <- rbind(color_data, data.frame(
      y_position = y_pos,
      region = region_type,
      color = feature_colors[region_type]
    ))
  }
  color_data$color[is.na(color_data$color)] <- "white"
  return(color_data)
}
feature_color_data <- create_region_color_data()

right_color_bar3_with_legend <- ggplot(feature_color_data, aes(x = 1, y = y_position, fill = region)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_manual(values = feature_colors, name = "Genomic Region") +
  scale_y_continuous(expand = c(0.02, 0), limits = c(y_range[1] - 0.5, y_range[2] + 0.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(
    plot.margin = margin(5, 5, 5, -5),
    legend.position = "right",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6)
  )

heatmap_legend <- get_legend(heatmap_grouped)
te_class_legend <- get_legend(right_color_bar_with_legend)
te_family_legend <- get_legend(right_color_bar2_with_legend)
region_legend <- get_legend(right_color_bar3_with_legend)
right_color_bar_no_legend <- right_color_bar_with_legend + theme(legend.position = "none")
right_color_bar2_no_legend <- right_color_bar2_with_legend + theme(legend.position = "none")
right_color_bar3_no_legend <- right_color_bar3_with_legend + theme(legend.position = "none")

top_row <- plot_grid(color_bar, NULL, NULL, NULL, 
                     ncol = 4, rel_widths = c(1, 0.05, 0.05, 0.05),
                     align = "h")

# bottom_row <- plot_grid(heatmap_no_legend, 
#                         right_color_bar_no_legend, 
#                         right_color_bar2_no_legend,
#                         right_color_bar3_no_legend,
#                         ncol = 4, rel_widths = c(1, 0.05, 0.05, 0.05),
#                         align = "h")
bottom_row <- plot_grid(heatmap_no_legend,
                        right_color_bar_no_legend,
                        right_color_bar2_no_legend,
                        ncol = 3, rel_widths = c(1, 0.05, 0.05),
                        align = "h")

main_plot <- plot_grid(top_row, bottom_row,
                       ncol = 1, rel_heights = c(0.02, 1),
                       align = "v")

all_legends <- plot_grid(heatmap_legend, 
                         te_class_legend, 
                         te_family_legend, 
                         ncol = 1, 
                         align = "v")

final_plot_with_all_legends <- plot_grid(main_plot, all_legends, 
                                         ncol = 2, rel_widths = c(1, 0.2))

print(final_plot_with_all_legends)
# ggsave('Fig2_All_species_Up_DE-TEs_feature_padj0.05log2FC1_v2_intergenic_top5.pdf', final_plot_with_legend, width = 10, height = 13, units = "cm")
# ggsave('Fig2s_All_species_Up_DE-TEs_feature_padj0.05log2FC1_v2_intergenic_intron_top100.pdf', final_plot_with_legend, width = 12, height = 26, units = "cm")
ggsave('Fig2s_All_species_Up_DE-TEs_feature_padj0.05log2FC1_v2_intergenic_intron_exon_top5.pdf', final_plot_with_all_legends, width = 13, height = 30, units = "cm")
# write.table(plot_data_spaced,"Fig2s_All_species_Up_DE-TEs_feature_padj0.05log2FC1_v2_intergenic_intron_e xon_top100_plot_data_spaced.csv",quote = F,row.names = F,sep = "\t")
# write.table(plot_data_spaced,"Fig2s_All_species_Up_DE-TEs_feature_padj0.05log2FC1_v2_intergenic_intron_top100_plot_data_spaced.csv",quote = F,row.names = F,sep = "\t")
write.csv(plot_data_spaced,"Fig2_All_species_Up_DE-TEs_feature_padj0.05log2FC1_v2_top5_plot_data_spaced.csv",quote = F,row.names = F)
