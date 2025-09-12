library(readxl)
library(ggplot2)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(gridExtra)
library(cowplot)

count1 <- read.csv("/data2t_2/hushu/02.DESeq2_analysis_TE_subfamily_feature/raw_data/new_all_sample_readscounts_matrix_combined_n4231_TE_subfamily_feature.txt", sep = "\t")
count2 <- read.csv("/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data_TE_subfamily_feature_v2/new_add_all_sample_readscounts_matrix_combined_n102_TE_subfamily_feature_v2.txt",sep = "\t")
count<-cbind(count1,count2)
metadata<-read_excel("/data2t_2/pathogen_TE_2025_New/01.new_raw_data/all_sample.xlsx")
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
top10_subfamilies <- c("LTR101_Mam_intergenic", "THE1C-int_intergenic", "BLACKJACK_intergenic", "MER5B_intergenic", "L1MA10_intergenic")
all_plot_data <- data.frame()
##########################
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
  all_plot_data <- rbind(all_plot_data, plot_data)
}

all_plot_data$Subfamily <- factor(all_plot_data$Subfamily, levels = top10_subfamilies)

p <- ggplot(all_plot_data, aes(x = Subfamily, y = Log2_Count, fill = Condition, color = Condition)) +
  geom_boxplot(
    outlier.size = 0.8, 
    width = 0.6, 
    position = position_dodge(0.7),
    alpha = 0,  
    outlier.alpha = 1  
  ) +
  scale_fill_manual(values = c("Pre" = "#4C9BCF", "Post" = "#D00732")) +
  scale_color_manual(values = c("Pre" = "#4C9BCF", "Post" = "#D00732")) + 
  scale_x_discrete(
    expand = expansion(mult = c(0.15, 0.15))  
  ) +
  scale_y_continuous(
    limits = c(0, 13),
    breaks = seq(0, 12, by = 2)
  ) +
  labs(
    title = "",
    x = "",
    y = "log2(Normalized Counts)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    axis.line = element_line(size = 0.5),
    axis.ticks = element_line(size = 0.5),
    plot.margin = margin(10, 10, 10, 10, "pt"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks.length.x = unit(0.1, "cm")
  )

for(i in 1:(length(top10_subfamilies)-1)) {
  p <- p + geom_vline(
    xintercept = i + 0.5, 
    color = "gray70", 
    size = 0.5, 
    linetype = "longdash"
  )
}

p_values <- c()
y_positions <- c()
line_colors <- c()

for(i in 1:length(top10_subfamilies)) {
  subfamily <- top10_subfamilies[i]
  subfamily_data <- all_plot_data[all_plot_data$Subfamily == subfamily, ]
  
  if(nrow(subfamily_data) > 0) {
    count_pre <- sum(subfamily_data$Condition == "Pre")
    count_post <- sum(subfamily_data$Condition == "Post")
    
    if(count_pre > 1 && count_post > 1) {
      stat_test <- t.test(Log2_Count ~ Condition, data = subfamily_data)
      p_value <- stat_test$p.value
      p_values[i] <- p_value
      max_y <- max(subfamily_data$Log2_Count, na.rm = TRUE)
      y_positions[i] <- min(max_y + 1, 11)
      line_colors[i] <- "black"
    } else {
      p_values[i] <- NA
      y_positions[i] <- 11
      line_colors[i] <- "lightgray"
    }
  }
}

for(i in 1:length(top10_subfamilies)) {
  if(!is.na(p_values[i])) {
    p <- p + annotate("segment", 
                      x = i - 0.35, xend = i + 0.35, 
                      y = y_positions[i] + 1, yend = y_positions[i] + 1,
                      size = 0.5, color = line_colors[i])
    p <- p + annotate("segment", 
                      x = i - 0.35, xend = i - 0.35, 
                      y = y_positions[i] + 0.8, yend = y_positions[i] + 1,
                      size = 0.5, color = line_colors[i])
    p <- p + annotate("segment", 
                      x = i + 0.35, xend = i + 0.35, 
                      y = y_positions[i] + 0.8, yend = y_positions[i] + 1,
                      size = 0.5, color = line_colors[i])
    p <- p + annotate("text", 
                      x = i, y = y_positions[i] + 1.6, 
                      label = format(p_values[i], scientific = TRUE, digits = 2),
                      size = 2.5, color = line_colors[i])
  }
}

ggsave(p,file="/data2t_2/pathogen_TE_2025_New/14.Whole_proteome/MHC/result-all-hla-191mb/infection_only_peptide_all_samples_TE_subfamilies_feature_Pre_Post_Normalized_Counts_boxplot_combine.pdf",width = 4, height = 3.5)
############################################################
library(dplyr)
library(readr)
library(stringr)

setwd("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_subfamily_feature/DE_results")

find_command <- 'find . -name "*all_res_TEs_feature.csv"'
file_paths <- system(find_command, intern = TRUE)

parse_file_info <- function(file_path) {
  celltype <- str_extract(file_path, "(?<=\\./)[^/]+")
  filename <- basename(file_path)
  parts <- str_split(str_remove(filename, "_all_res_TEs_feature\\.csv$"), "_")[[1]]

  pathogen_type <- parts[1]
  species <- parts[2]
  dataset_id <- parts[3]
  
  post_time <- str_extract(filename, "Post_([0-9]+m)") %>% str_remove("Post_")
  pre_time <- str_extract(filename, "Pre_([0-9]+m)") %>% str_remove("Pre_")
  
  return(data.frame(
    file_path = file_path,
    celltype = celltype,
    comparison = paste0(pathogen_type, "_", species, "_", dataset_id, "_", 
                        "Post_", post_time, "_vs_Pre_", pre_time),
    stringsAsFactors = FALSE
  ))
}

file_info <- do.call(rbind, lapply(file_paths, parse_file_info))

read_and_annotate <- function(file_path, file_info_row) {
  data <- read_csv(file_path, show_col_types = FALSE)
  data$celltype <- file_info_row$celltype
  data$comparison <- file_info_row$comparison
  return(data)
}

all_data_list <- list()
for(i in 1:nrow(file_info)) {
  tryCatch({
    all_data_list[[i]] <- read_and_annotate(file_info$file_path[i], file_info[i,])
  }, error = function(e) {
    cat(file_info$file_path[i], e$message, "\n")
  })
}

all_data_list <- all_data_list[!sapply(all_data_list, is.null)]
combined_data <- do.call(rbind, all_data_list)
colnames(combined_data)[1] <- "TEs_feature"

#LTR101_Mam_intergenic
LTR101_Mam_intergenic_df<-combined_data[combined_data$TEs_feature=="LTR101_Mam_intergenic",]
Up_LTR101_Mam_intergenic_df <- LTR101_Mam_intergenic_df[LTR101_Mam_intergenic_df$log2FoldChange > 1 & LTR101_Mam_intergenic_df$pvalue <= 0.05, ]
#THE1C-int_intergenic
THE1C_int_intergenic_df<-combined_data[combined_data$TEs_feature=="THE1C-int_intergenic",]
Up_THE1C_int_intergenic_df <- THE1C_int_intergenic_df[THE1C_int_intergenic_df$log2FoldChange > 1 & THE1C_int_intergenic_df$pvalue <= 0.05, ]
#BLACKJACK_intergenic
BLACKJACK_intergenic_df<-combined_data[combined_data$TEs_feature=="MLT-int_intergenic",]
Up_BLACKJACK_intergenic_df <- BLACKJACK_intergenic_df[BLACKJACK_intergenic_df$log2FoldChange > 1 & BLACKJACK_intergenic_df$pvalue <= 0.05, ]
#MER5B_intergenic
MER5B_intergenic_df<-combined_data[combined_data$TEs_feature=="MER5B_intergenic",]
Up_MER5B_intergenic_df <- MER5B_intergenic_df[MER5B_intergenic_df$log2FoldChange > 1 & MER5B_intergenic_df$pvalue <= 0.05, ]
#L1MA10_intergenic
L1MA10_intergenic_df<-combined_data[combined_data$TEs_feature=="L1MA10_intergenic",]
Up_L1MA10_intergenic_df <- L1MA10_intergenic_df[L1MA10_intergenic_df$log2FoldChange > 1 & L1MA10_intergenic_df$pvalue <= 0.05, ]

Virus_IAV_11_Post_720m_vs_Pre_0m<-combined_data[combined_data$comparison=="Virus_IAV_11_Post_720m_vs_Pre_0m",]
Bacteria_NTHi_1_Post_360m_vs_Pre_360m<-combined_data[combined_data$comparison=="Bacteria_NTHi_1_Post_360m_vs_Pre_360m",]
Bacteria_SalTy_4_Post_120m_vs_Pre_120m<-combined_data[combined_data$comparison=="Bacteria_SalTy_4_Post_120m_vs_Pre_120m",]
Virus_SARSCoV2_1_Post_2240m_vs_Pre_2240m<-combined_data[combined_data$comparison=="Virus_SARSCoV2_1_Post_2240m_vs_Pre_2240m",]

# res_final_ord<-Virus_IAV_11_Post_720m_vs_Pre_0m
# res_final_ord<-Bacteria_NTHi_1_Post_360m_vs_Pre_360m
# res_final_ord<-Bacteria_SalTy_4_Post_120m_vs_Pre_120m
res_final_ord<-Virus_SARSCoV2_1_Post_2240m_vs_Pre_2240m
res_intergenic <- res_final_ord[grep("_intergenic$", res_final_ord$TEs_feature), ]
res_intergenic$TEs_feature <- gsub("_intergenic$", "", res_intergenic$TEs_feature)

te_classification <- read.csv("/data2t_2/hushu/02.DESeq2_analysis/plots/Origin_data.csv", header = FALSE, stringsAsFactors = FALSE)
colnames(te_classification) <- c("Classification")
te_classification$TE_Subfamily <- str_extract(te_classification$Classification, "^[^:]+")
te_classification$TE_Family <- str_extract(te_classification$Classification, "(?<=:)[^:]+(?=:[^:]*$)")
te_classification$TE_Class <- str_extract(te_classification$Classification, "[^:]+$")

res_intergenic_class <- merge(res_intergenic, te_classification,  by.x="TEs_feature",by.y="TE_Subfamily", all.x = TRUE)
colors <- c("DNA" = "#f8d196","LINE" = "#d7a9cb", "LTR" = "#6fa4af","SINE" = "#8290bb","Retroposon" = "#277899")

fixed_tes_to_label <- c("L1MA10", "LTR101_Mam", "MER5B", "BLACKJACK", "THE1C-int")

res_intergenic_class$significant <- ifelse(
  res_intergenic_class$padj < 0.05 & abs(res_intergenic_class$log2FoldChange) > 1,
  "Significant",
  "Not Significant"
)

top_to_label <- res_intergenic_class %>%
  filter(significant == "Significant") %>%
  arrange(padj) %>%
  head(15)  
res_intergenic_class$label <- ifelse(
  res_intergenic_class$TEs_feature %in% fixed_tes_to_label,
  res_intergenic_class$TEs_feature,
  ""
)

res_intergenic_class$padj_adjusted <- ifelse(res_intergenic_class$padj == 0 | is.na(res_intergenic_class$padj), 1e-300, res_intergenic_class$padj)

p <- ggplot(res_intergenic_class, aes(x = log2FoldChange, y = -log10(padj_adjusted))) +
  geom_point(data = subset(res_intergenic_class, significant == "Not Significant"),
             color = "grey80", size = 0.8, alpha = 0.5) +
  geom_point(data = subset(res_intergenic_class, significant == "Significant"),
             aes(color = TE_Class), size = 1.2, alpha = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", alpha = 0.5) +
  geom_text_repel(aes(label = label), 
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = 'grey50',
                  max.overlaps = Inf,
                  size = 3) +
  scale_color_manual(values = colors, name = "TE Class") +
  labs(x = "Log2 Infection/Mock",
       y = "-Log10 padj",
      title = "Virus_SARSCoV2_1_Post_2240m_vs_Pre_2240m intergenic TE") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        title = element_text(size = 10),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))

print(p)