library(ggplot2)
library(dplyr)
library(readr)
library(scales)
library(gridExtra)
library(cowplot)
library(stringr)
library(ggrepel)
library(data.table)
setwd("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_subfamily_feature/DE_results/")
raw_count <- read.csv("/data2t_2/hushu/02.DESeq2_analysis_TE_subfamily_feature/raw_data/new_all_sample_readscounts_matrix_combined_n4231_TE_subfamily_feature.txt",sep = "\t",stringsAsFactors = FALSE)
all_files <- list.files(pattern = "sig_DE-TEs_feature_padj0.05log2FC1\\.csv$", recursive = TRUE)
files <- all_files[!grepl("^01\\.", basename(all_files)) & !grepl("^All", basename(all_files))]

res_list <- lapply(files, function(f) {
  df <- read.csv(f, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
})

all_genes <- union(raw_count$GeneID, unique(unlist(lapply(res_list, rownames))))

sig_mat <- sapply(res_list, function(df) {
  out <- rep(NA, length(all_genes))
  names(out) <- all_genes
  common <- intersect(all_genes, rownames(df))
  out[common] <- df[common, "sig"]
  return(out)
})

gene_status_map <- apply(sig_mat, 1, function(statuses) {
  statuses <- statuses[!is.na(statuses)]
  if (length(statuses) == 0) {
    return("/")
  } else if (all(statuses == "Up")) {
    return("Up")
  } else if (all(statuses == "Down")) {
    return("Down")
  } else {
    return("Mixed")
  }
})


final_df <- data.frame(
  GeneID    = raw_count$GeneID,
  DE_status = gene_status_map[ raw_count$GeneID ],
  stringsAsFactors = FALSE
)

write.csv(final_df, "/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_subfamily_feature/plots/All_DE_status_sig_DE-TEs_feature_padj0.05log2FC1.csv", row.names = FALSE)

setwd("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_subfamily_feature/plots/")

p1 <- read_csv("All_DE_status_sig_DE-TEs_feature_padj0.05log2FC1.csv")
p1 <- p1 %>%filter(!grepl("_intron$|_exon$", GeneID))
# p1_status_counts <- p1 %>%
#   dplyr::count(DE_status) %>%
#   dplyr::mutate(percentage = n / sum(n) * 100,
#          label = paste0(DE_status, "\n", round(percentage, 2), "%"),
#          category = ifelse(DE_status == "/", "non-DE_TE_Loci", "DE_TE_Loci"))


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
  labs(title = "TE subfamily Differential Expression Status Distribution",
       subtitle = paste0("Total: ", sum(p1_status_counts$n), " loci (padj<0.05, log2FC>1)"),
       x = NULL, y = NULL, fill = "DE Status") +
  geom_text(aes(label = paste0(display_label, "\n", n, " (", round(percentage, 2), "%)")), 
            position = position_stack(vjust = 0.5),
            color = "black", size = 1)+
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 7),
    plot.subtitle = element_text(hjust = 0.5, size = 4),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )
p1_fig
ggsave("Fig2_All_DE_status_sig_DE_TEs_feature_padj0.05log2FC1_pie_intergenic.pdf", p1_fig, width = 2, height = 2)

#######################################################################################
#####  v2 2025-07-28 ##########################################################
####################################################################################### 
setwd("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_subfamily_feature/plots/")
p1 <- fread("Recurrent_Up_TEs_feature_intergenic_intron_exon_UniqPathogens_gt1.csv")
p1<-p1[p1$UniqPathogens>=5,]
raw_count <- read.csv("/data2t_2/hushu/02.DESeq2_analysis_TE_subfamily_feature/raw_data/new_all_sample_readscounts_matrix_combined_n4231_TE_subfamily_feature.txt",sep = "\t",stringsAsFactors = FALSE)
p1_parsed <- p1 %>%
  mutate(
    region_type = str_extract(Subfamily, "_(exon|intron|intergenic)$"),
    region_type = str_remove(region_type, "^_"),
    te_type = str_remove(Subfamily, "_(exon|intron|intergenic)$")
  )

p1_counts <- p1_parsed %>%
  count(region_type, name = "up_count")

raw_counts <- raw_count %>%
  mutate(region_type = str_extract(GeneID, "_(exon|intron|intergenic)$")) %>%
  mutate(region_type = str_remove(region_type, "^_")) %>%
  count(region_type, name = "total_count")


combined_counts <- raw_counts %>%
  left_join(p1_counts, by = "region_type") %>%
  mutate(
    up_count = ifelse(is.na(up_count), 0, up_count),
    non_up_count = total_count - up_count,
    up_label = paste0("Up_", region_type),
    non_up_label = "Non-Up"
  )

final_stats <- tibble(
  category = c(combined_counts$up_label),
  count = c(combined_counts$up_count)
) %>%
  mutate(
    percentage = count / sum(count) * 100
  )


color_palette <- c(
  "Up_exon" = "#B82132",        
  "Up_intron" = "#D2665A",      
  "Up_intergenic" = "#F2B28C",  
  "Non-Up" = "#D2D8E7"          
)
final_stats <- final_stats %>%
  mutate(display_label = category)  

p1_fig <- ggplot(final_stats, aes(x = 2, y = count, fill = category)) +  
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +  
  scale_fill_manual(values = color_palette) +
  # labs(title = "Recurrent Up TE Loci(intergenic) Unique Pathogens>=1",
  labs(title = "Recurrent Up TE Loci(intergenic) Unique Pathogens>=5",
       subtitle = paste0("Total: ", sum(final_stats$count), " loci (padj<0.05, log2FC>1)"),
       x = NULL, y = NULL, fill = "DE Status") +
  geom_text(aes(label = paste0(category, "\n", scales::comma(count), " (", round(percentage, 1), "%)")), 
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 1)+
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 7),
    plot.subtitle = element_text(hjust = 0.5, size = 4),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

ggsave("Fig2_recurrent_up_TE_loci_all_Up_status_padj0.05log2FC1_pies_gt5.pdf", p1_fig, width = 2, height = 2)
# ggsave("Fig2_recurrent_up_TE_loci_all_Up_status_padj0.05log2FC1_pies_gt1.pdf", p1_fig, width = 2, height = 2)

p1_parsed<-p1_parsed[p1_parsed$region_type=="intergenic",]
te_bed<-fread("/data2t_2/pathogen_TE_2025_New/01.Genomic_features_Gencode/hg38_TE_anno_custom_v20240110_0based_with_feature_summary_v2.bed")
gene_class_mapping <- te_bed %>%
  select(gene_id, class_id) %>%
  distinct()

head(gene_class_mapping)

p1_with_class <- p1_parsed %>%
  left_join(gene_class_mapping, by = c("te_type" = "gene_id"))


unmatched <- p1_with_class %>% filter(is.na(class_id))

class_percentages <- round(100 * class_counts / sum(class_counts), 2)
class_summary <- data.frame(class_counts) %>%
  rename(class_id = Var1, count = Freq) %>%
  mutate(
    percentage = round(100 * count / sum(count), 2),
    label = paste0(class_id, ": ", percentage, "%")
  ) %>%
  arrange(desc(count)) %>%  
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
    title = "Recurrent Up TE Loci(intergenic) Unique Pathogens>=5 Class Distribution") +
  geom_text(
    aes(y = ypos, label = paste0(class_id, "\n", count, " (", percentage, "%)")),
    color = "black",
    size = 1
  ) +
  scale_fill_manual(values = colors)

p2
# ggsave("Fig2_recurrent_up_TEs_feature_intergenic_Class_Distribution_padj0.05log2FC1_gt1.pdf", p2, width = 2, height = 2)
ggsave("Fig2_recurrent_up_TEs_feature_intergenic_Class_Distribution_padj0.05log2FC1_gt5.pdf", p2, width = 2, height = 2)

