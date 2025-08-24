library(readxl)
library(ggplot2)
library(dplyr)
library(DESeq2)

old_count <- read.csv("/data2t_2/hushu/02.DESeq2_analysis_TE_subfamily_feature_v2/raw_data/new_all_sample_readscounts_matrix_combined_n4231_TE_subfamily_feature_v2.txt", sep = "\t")
new_count<-read.csv("/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data_TE_subfamily_feature_v2/new_add_all_sample_readscounts_matrix_combined_n102_TE_subfamily_feature_v2.txt",sep = "\t")

count<-cbind(old_count,new_count)
metadata <- read_excel("/data2t_2/pathogen_TE_2025_New/01.new_raw_data/all_sample.xlsx")

gene_id_col <- colnames(count)[1]
count_data <- count[, !(colnames(count) %in% gene_id_col)]
gene_ids <- count[[gene_id_col]]

sample_data <- data.frame(
  Sample = metadata$Sample_Name,
  Condition = metadata$Infection_State,
  stringsAsFactors = FALSE)

rownames(sample_data) <- sample_data$Sample

common_samples <- intersect(colnames(count_data), sample_data$Sample)

count_data <- count_data[, common_samples]
sample_data_filtered <- sample_data[sample_data$Sample %in% common_samples, ]

sample_data_filtered <- sample_data_filtered[match(colnames(count_data), sample_data_filtered$Sample), ]
rownames(sample_data_filtered) <- sample_data_filtered$Sample

count_matrix <- as.matrix(count_data)
rownames(count_matrix) <- gene_ids

sample_data_filtered$Condition <- factor(sample_data_filtered$Condition, levels = c("Pre", "Post"))

dataset_dds1 <- DESeqDataSetFromMatrix(
  countData = count_matrix, 
  colData = sample_data_filtered, 
  design = ~ Condition
)

dataset_dds2 <- estimateSizeFactors(dataset_dds1)

dataset_count_norm <- counts(dataset_dds2, normalized = TRUE)

top10_subfamilies <- c("AluYk11_intergenic","MamGypLTR2b_intergenic","MamRep38_intergenic","L1MA10_intergenic","MamGypLTR3_intergenic",
                       "LTR54_intergenic","THE1C-int_intergenic","HERV9NC-int_intergenic","LTR46-int_intergenic","MER110A_intergenic","MER5B_intergenic")

combined_plot_data <- data.frame()

for(subfamily in top10_subfamilies) {
  subfamily_idx <- which(rownames(dataset_count_norm) == subfamily)
  
  subfamily_norm_counts <- dataset_count_norm[subfamily_idx, ]
  
  plot_data <- data.frame(
    Sample = colnames(dataset_count_norm),
    Normalized_Count = as.numeric(subfamily_norm_counts),
    Condition = sample_data_filtered$Condition,
    Subfamily = subfamily
  )
  
  plot_data <- plot_data[plot_data$Normalized_Count > 0, ]
  
  plot_data$Log2_Count <- log2(plot_data$Normalized_Count + 1)
  
  combined_plot_data <- rbind(combined_plot_data, plot_data)
}

combined_plot_data$Subfamily_short <- gsub("_intergenic", "", combined_plot_data$Subfamily)

combined_plot_data$Subfamily_short <- factor(combined_plot_data$Subfamily_short, 
                                             levels = gsub("_intergenic", "", top10_subfamilies))

p_values_data <- combined_plot_data %>%
  group_by(Subfamily_short) %>%
  summarise(
    p_value = ifelse(length(unique(Condition)) == 2, 
                     t.test(Log2_Count ~ Condition)$p.value, 
                     NA),
    y_max = max(Log2_Count),
    .groups = 'drop'
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value >= 0.05 ~ "ns",
      TRUE ~ ""
    )
  )

library(gghalves)

main_plot <- ggplot(combined_plot_data, aes(x = Subfamily_short, y = Log2_Count, fill = Condition)) +
  geom_half_violin(data = combined_plot_data[combined_plot_data$Condition == "Pre",], 
                   aes(fill = Condition), side = "l", alpha = 0.7,width = 1, color = NA) +
  geom_half_violin(data = combined_plot_data[combined_plot_data$Condition == "Post",], 
                   aes(fill = Condition), side = "r", alpha = 0.7,width = 1,color = NA) +
  stat_summary(data = combined_plot_data[combined_plot_data$Condition == "Pre",],
               fun = median, geom = "point", 
               position = position_nudge(x = -0.07), 
               color = "black", size = 0.3) +
  stat_summary(data = combined_plot_data[combined_plot_data$Condition == "Post",],
               fun = median, geom = "point", 
               position = position_nudge(x = 0.07), 
               color = "black", size = 0.3) +
  stat_summary(data = combined_plot_data[combined_plot_data$Condition == "Pre",],
               fun.min = function(x) quantile(x, 0.25),
               fun.max = function(x) quantile(x, 0.75),
               geom = "linerange", 
               position = position_nudge(x = -0.07),
               color = "black", size = 0.2) +
  stat_summary(data = combined_plot_data[combined_plot_data$Condition == "Post",],
               fun.min = function(x) quantile(x, 0.25),
               fun.max = function(x) quantile(x, 0.75),
               geom = "linerange", 
               position = position_nudge(x = 0.07),
               color = "black", size = 0.2) +
  scale_fill_manual(values = c("Pre" = "#4c9bcf", "Post" = "#d00732")) +
  labs(
    y = "log2(DESeq2 Normalized Counts + 1)",
    x = "TE Subfamily",
    title = "Top 10 Intergenic TE Subfamily Expression Comparison (Pre vs Post)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 7, angle = 90,hjust = 1,vjust = 0.5),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(size = 0.3),  
    axis.line.y = element_line(size = 0.3),  
    legend.position = "top",
    legend.title = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 6)
  )

main_plot
for(i in 1:nrow(p_values_data)) {
  subfamily <- p_values_data$Subfamily_short[i]
  significance <- p_values_data$significance[i]
  y_pos <- p_values_data$y_max[i] * 1.05
  
  main_plot <- main_plot +
    annotate("text", x = i, y = y_pos, 
             label = significance, 
             size = 3)
}

ggsave("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_subfamily_feature/plots/Fig2_Top10_TE_subfamilies_feature_Pre_Post_Normalized_Counts_violin_v2.pdf",main_plot, width = 2.5, height = 3)
