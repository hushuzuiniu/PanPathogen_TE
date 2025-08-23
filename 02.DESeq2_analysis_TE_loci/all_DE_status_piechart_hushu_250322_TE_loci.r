library(ggplot2)
library(dplyr)
library(readr)
library(scales)
library(gridExtra)
library(cowplot)
library(stringr)
library(ggrepel)
library(data.table)
setwd("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_loci/plots/")
p1 <- fread("All_DE_status_sig_DE-TEs_Loci_padj0.05log2FC1.csv")

# p1_status_counts <- p1 %>%
#   count(DE_status) %>%
#   mutate(percentage = n / sum(n) * 100,
#          label = paste0(DE_status, "\n", round(percentage, 2), "%"),
#          category = ifelse(DE_status == "/", "non-DE_TE_Loci", "DE_TE_Loci"))
# 
# p1_counts <- p1_status_counts %>%
#   group_by(category) %>%
#   summarise(n = sum(n)) %>%
#   mutate(percentage = n / sum(n) * 100,
#          label = paste0(category, "\n", round(percentage, 2), "%"))
# 
# p1_fig <- ggplot(p1_counts, aes(x = "", y = n, fill = category)) +
#   geom_bar(stat = "identity", width = 1) +
#   coord_polar(theta = "y") +
#   scale_fill_manual(values = c("non-DE_TE_Loci" = "#D2D8E7", "DE_TE_Loci" = "#B187AB")) +
#   labs(title = "TE_DE_Loci Distribution(padj0.05log2FC1)",
#        # subtitle = paste0("Total records: ", sum(p0.5_counts$n)),
#        x = NULL, y = NULL, fill = "Category") +
#   geom_text(aes(label = paste0(category, "\n", n, " (", round(percentage, 2), "%)")), 
#             position = position_stack(vjust = 0.5),
#             color = "white", fontface = "bold", size = 5) +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
#         legend.position = "none",
#         legend.title = element_text(face = "bold"),
#         panel.grid = element_blank(),
#         axis.text = element_blank())
# p1_fig

p1_status_counts <- p1 %>%
  count(DE_status) %>%
  mutate(percentage = n / sum(n) * 100,
         display_label = ifelse(DE_status == "/", "Non-DE", DE_status))


p1_fig <- ggplot(p1_status_counts, aes(x = 2, y = n, fill = display_label)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +  
  scale_fill_manual(values = c(
    "Non-DE" = "#D2D8E7",     
    "Up" = "#811d3f",         
    "Down" = "#2c6aa0",       
    "Mixed" = "#F39C12"       
  )) +
  labs(title = "TE Differential Expression Status Distribution",
       subtitle = paste0("Total: ", sum(p1_status_counts$n), " loci (padj<0.05, log2FC>1)"),
       x = NULL, y = NULL, fill = "DE Status") +
  geom_text(aes(label = paste0(display_label, "\n", n, " (", round(percentage, 2), "%)")), 
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 1)+
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 7),
    plot.subtitle = element_text(hjust = 0.5, size = 4),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )
p1_fig
ggsave("Fig2_All_DE_status_sig_DE_TEs_Loci_padj0.05log2FC1_pies.pdf", p1_fig, width = 2, height = 2)

#######################################################################################
#####  v2 2025-07-27 ##########################################################
####################################################################################### 
setwd("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_loci/plots/")
p1 <- fread("Recurrent_Up_TE_Loci_intergenic_intron_exon_UniqPathogens_gt1.csv")
p1<-p1[p1$UniqPathogens>=9,]
te_bed <-fread("/data2t_2/pathogen_TE_2025_New/01.Genomic_features_Gencode/hg38_TE_anno_custom_v20240110_0based_with_feature_summary_v2.bed")
setDT(te_bed)
setDT(p1)

p1[, c("te_type", "chrom", "position") := {
  parts <- strsplit(Subfamily, "_chr")
  te_type <- sapply(parts, function(x) x[1])  
  chrom_pos <- sapply(parts, function(x) paste0("chr", x[2]))  
  chrom_pos_split <- strsplit(chrom_pos, ":")
  list(
    te_type = te_type,
    chrom = sapply(chrom_pos_split, function(x) x[1]),
    position = as.integer(sapply(chrom_pos_split, function(x) x[2]))-1
  )
}]
result <- te_bed[p1, on = .(gene_id = te_type, chrom, start = position), nomatch = 0]

final_result <- result[, .(Subfamily, UniqPathogens, transcript_id, chrom, start, end, 
                           class_id, family_id, prioritized_feature)]

up_feature_counts <- final_result %>%
  count(prioritized_feature, name = "n")

total_te_bed <- nrow(te_bed)
total_up <- nrow(final_result)
non_up_count <- total_te_bed - total_up

# p1_status_counts <- tibble(
#   DE_status = c(paste0("Up_", up_feature_counts$prioritized_feature), "Non-Up"),
#   n = c(up_feature_counts$n, non_up_count)
# ) %>%
#   mutate(
#     percentage = n / sum(n) * 100,
#     display_label = DE_status
#   )

p1_status_counts <- tibble(
  DE_status = paste0("Up_", up_feature_counts$prioritized_feature),
  n = up_feature_counts$n
) %>%
  mutate(
    percentage = n / sum(n) * 100,
    display_label = DE_status
  )

color_palette <- c(
  "Up_exon" = "#B82132",        
  "Up_intron" = "#D2665A",      
  "Up_intergenic" = "#F2B28C",  
  "Non-Up" = "#D2D8E7"          
)

p1_fig <- ggplot(p1_status_counts, aes(x = 2, y = n, fill = display_label)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +  
  scale_fill_manual(values = color_palette) +
  # labs(title = "Recurrent Up TE Loci(intergenic) Unique Pathogens>=1",
       labs(title = "Recurrent Up TE Loci(intergenic) Unique Pathogens>=9",
       subtitle = paste0("Total: ", sum(p1_status_counts$n), " loci (padj<0.05, log2FC>1)"),
       x = NULL, y = NULL, fill = "DE Status") +
  geom_text(aes(label = paste0(display_label, "\n", n, " (", round(percentage, 2), "%)")), 
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 1)+
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 7),
    plot.subtitle = element_text(hjust = 0.5, size = 4),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )
p1_fig
# ggsave("Fig2_recurrent_up_TE_loci_all_Up_status_padj0.05log2FC1_pies_gt5.pdf", p1_fig, width = 2, height = 2)
ggsave("Fig2_recurrent_up_TE_loci_all_Up_status_padj0.05log2FC1_pies_gt9.pdf", p1_fig, width = 2, height = 2)
# ggsave("Fig2_recurrent_up_TE_loci_all_Up_status_padj0.05log2FC1_pies_gt1.pdf", p1_fig, width = 2, height = 2)
# ggsave("Fig2_recurrent_up_TE_loci_intergenic_status_padj0.05log2FC1_pies_gt5.pdf", p1_fig, width = 2, height = 2)
# ggsave("Fig2_recurrent_up_TE_loci_intergenic_intron_status_padj0.05log2FC1_pies.pdf", p1_fig, width = 2, height = 2)
# ggsave("Fig2_recurrent_up_TE_loci_intergenic_intron_status_padj0.05log2FC1_pies_gt5.pdf", p1_fig, width = 2, height = 2)

final_result<-final_result[final_result$prioritized_feature=="intergenic",]
class_counts <- table(final_result$class_id)

class_percentages <- round(100 * class_counts / sum(class_counts), 2)
class_summary <- data.frame(class_counts) %>%
  rename(class_id = Var1, count = Freq) %>%
  mutate(percentage = round(100 * count / sum(count), 2),
         label = paste0(class_id, ": ", percentage, "%"))

class_summary <- class_summary %>%
  arrange(desc(class_id)) %>%
  mutate(ypos = cumsum(count) - 0.5 * count)

colors <- c("DNA" = "#f8d196", "LINE" = "#d7a9cb", "LTR" = "#6fa4af",
            "SINE" = "#8290bb", "Retroposon" = "#277899")

p2 <- ggplot(class_summary, aes(x = "", y = count, fill = class_id)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
    legend.position = "none"  
  ) +
  labs(
    # title = "Recurrent Up TE Loci(intergenic) Unique Pathogens>=1 Class Distribution") +
    title = "Recurrent Up TE Loci(intergenic) Unique Pathogens>=9 Class Distribution") +
  geom_text(
    aes(y = ypos, label = paste0(class_id, "\n", count, " (", percentage, "%)")),
    color = "black",
    size = 1
  ) +
  scale_fill_manual(values = colors)

p2
# ggsave("Fig2_recurrent_up_TE_loci_intergenic_Class_Distribution_padj0.05log2FC1_gt5.pdf", p2, width = 2, height = 2)
ggsave("Fig2_recurrent_up_TE_loci_intergenic_Class_Distribution_padj0.05log2FC1_gt9.pdf", p2, width = 2, height = 2)
# ggsave("Fig2_recurrent_up_TE_loci_intergenic_Class_Distribution_padj0.05log2FC1_gt1.pdf", p2, width = 2, height = 2)
# ggsave("Fig2_recurrent_up_TE_loci_intergenic_intron_Class_Distribution_padj0.05log2FC1.pdf", p2, width = 2, height = 2)
# ggsave("Fig2_recurrent_up_TE_loci_intergenic_intron_Class_Distribution_padj0.05log2FC1_gt5.pdf", p2, width = 2, height = 2)
