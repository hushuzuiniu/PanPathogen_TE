library(ggplot2)
library(reshape2)
library(patchwork)
library(scales)
library(dplyr)
All_TEs_fc<-fread("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_subfamily_feature/plots/All_celltypes_combined_log2FoldChange_sig_DE-TEs_feature.csv")
All_Gene_fc<-fread("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_Gene/plots/All_celltypes_combined_log2FoldChange_sig_DE-Gene.csv")
gtf<-fread("/data2t_2/hushu/01.Genomic_features_Gencode/hg38.p13.gene.anno.gtf")
extract_gene_info <- function(string) {
  gene_id_pattern <- "gene_id \\\"(.*?)\\\";"
  gene_id <- sub(gene_id_pattern, "\\1", regmatches(string, regexpr(gene_id_pattern, string)))
  gene_name_pattern <- "gene_name \\\"(.*?)\\\";"
  gene_name_match <- regexpr(gene_name_pattern, string)
  if (gene_name_match != -1) {
    gene_name <- sub(gene_name_pattern, "\\1", regmatches(string, gene_name_match))
  } else {
    gene_name <- NA
  }
  return(data.frame(gene_id = gene_id, gene_name = gene_name, stringsAsFactors = FALSE))
}

gene_info_rows <- !duplicated(sub("gene_id \\\"(.*?)\\\";.*", "\\1", gtf$V9))
gene_info <- do.call(rbind, lapply(gtf$V9[gene_info_rows], extract_gene_info))
gene_id_to_name <- setNames(gene_info$gene_name, gene_info$gene_id)
All_Gene_fc$GeneName <- gene_id_to_name[All_Gene_fc$GeneID]

ZNF287_Gene_fc<-All_Gene_fc[All_Gene_fc$GeneName=="ZNF287",]
L1MA10_intergenic_fc<-All_TEs_fc[All_TEs_fc$GeneID=="L1MA10_intergenic",]

gene_values <- as.numeric(ZNF287_Gene_fc[1, -c(1, ncol(ZNF287_Gene_fc)), with = FALSE])
te_values <- as.numeric(L1MA10_intergenic_fc[1, -1])

dataset_names <- colnames(ZNF287_Gene_fc)[-c(1, ncol(ZNF287_Gene_fc))]  
plot_data <- data.frame(
  Dataset = rep(dataset_names, 2),
  log2FC = c(gene_values, te_values),
  Type = rep(c("ZNF287_Gene", "L1MA10_TE"), each = length(dataset_names))
)

plot_data_clean <- plot_data[!is.na(plot_data$log2FC), ]

gene_data <- data.frame(
  Dataset = dataset_names,
  Gene_log2FC = gene_values,
  TE_log2FC = te_values
)

paired_data <- gene_data[!is.na(gene_data$Gene_log2FC) & !is.na(gene_data$TE_log2FC), ]
paired_data <- paired_data[order(paired_data$TE_log2FC), ]
paired_data$Index <- 1:nrow(paired_data)

plot_data_long <- data.frame(
  Index = rep(paired_data$Index, 2),
  Dataset = rep(paired_data$Dataset, 2),
  log2FC = c(paired_data$TE_log2FC, paired_data$Gene_log2FC),
  Type = rep(c("L1MA10_TE", "ZNF287_Gene"), each = nrow(paired_data))
)

linear_trend_plot <- ggplot(plot_data_long, aes(x = Index, y = log2FC, color = Type)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +  
  scale_color_manual(values = c("L1MA10_TE" = "#c98d22", "ZNF287_Gene" = "#756fae")) +
  labs(
    title = "ZNF287 Gene vs L1MA10 TE Linear Trends",
    subtitle = "Ordered by L1MA10 TE log2FoldChange",
    x = "Sample Index (ordered by TE expression)",
    y = "log2FoldChange",
    color = "Gene/TE Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "top"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

print(linear_trend_plot)

te_lm <- lm(TE_log2FC ~ Index, data = paired_data)
gene_lm <- lm(Gene_log2FC ~ Index, data = paired_data)

gene_data <- plot_data[plot_data$Type == "ZNF287_Gene", ]
te_data <- plot_data[plot_data$Type == "L1MA10_TE", ]

paired_data <- merge(gene_data[, c("Dataset", "log2FC")], 
                     te_data[, c("Dataset", "log2FC")], 
                     by = "Dataset", 
                     suffixes = c("_Gene", "_TE"))

paired_data_clean <- paired_data[!is.na(paired_data$log2FC_Gene) & 
                                   !is.na(paired_data$log2FC_TE), ]

paired_t_test <- t.test(paired_data_clean$log2FC_Gene, 
                        paired_data_clean$log2FC_TE, 
                        paired = TRUE)
cor_test_pearson <- cor.test(paired_data_clean$log2FC_Gene, 
                             paired_data_clean$log2FC_TE, 
                             method = "pearson")

cor_test_spearman <- cor.test(paired_data_clean$log2FC_Gene, 
                              paired_data_clean$log2FC_TE, 
                              method = "spearman")


paired_data_clean$Difference <- paired_data_clean$log2FC_Gene - paired_data_clean$log2FC_TE
paired_data_clean$Index <- 1:nrow(paired_data_clean)


paired_plot_data <- data.frame(
  Index = rep(1:nrow(paired_data_clean), 2),
  Dataset = rep(paired_data_clean$Dataset, 2),
  log2FC = c(paired_data_clean$log2FC_Gene, paired_data_clean$log2FC_TE),
  Type = rep(c("ZNF287_Gene", "L1MA10_TE"), each = nrow(paired_data_clean))
)

paired_line_plot <- ggplot(paired_plot_data, aes(x = Type, y = log2FC, group = Index)) +
  geom_line(alpha = 0.8, color = "gray60") +
  geom_point(aes(color = Type), size = 0.5, alpha = 0.8) +
  scale_color_manual(values = setNames(c("#756fae", "#c98d22"), 
                                       c(paste0("ZNF287_Gene"), paste0("L1MA10_TE")))) +
  coord_cartesian(ylim = c(-max(abs(range(paired_plot_data$log2FC, na.rm = TRUE))), 
                           max(abs(range(paired_plot_data$log2FC, na.rm = TRUE))))) +
  labs(
    title = paste("Paired Expression Comparison ZNF287 vs L1MA10"),
    subtitle = paste("Mean difference =", round(paired_t_test$estimate, 3), 
                     ", p =", format.pval(paired_t_test$p.value, digits = 4)),
    x = "",
    y = "log2FoldChange",
    caption = "Lines connect paired observations from same dataset"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 7, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 6),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    legend.position = "none",
    axis.line = element_line(color = "black", size = 0.5),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5),
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

print(paired_line_plot)


# ggsave("/data2t_2/pathogen_TE_2025_New/08.TF_prediction/ZNF287_vs_L1MA10_paired_logfc_plot.pdf", paired_line_plot, width = 2, height = 2)


analyze_gene_te_pair <- function(gene_name, te_name, 
                                 gene_data = All_Gene_fc, 
                                 te_data = All_TEs_fc,
                                 output_dir = NULL,
                                 sort_by_correlation = TRUE) { 
  
  gene_fc <- gene_data[gene_data$GeneName == gene_name, ]
  te_fc <- te_data[te_data$GeneID == te_name, ]

  gene_values <- as.numeric(gene_fc[1, -c(1, ncol(gene_fc)), with = FALSE])
  te_values <- as.numeric(te_fc[1, -1])
  
  dataset_names <- colnames(gene_fc)[-c(1, ncol(gene_fc))]
  
  gene_data_df <- data.frame(
    Dataset = dataset_names,
    Gene_log2FC = gene_values,
    TE_log2FC = te_values
  )
  
  paired_data_clean <- gene_data_df[!is.na(gene_data_df$Gene_log2FC) & 
                                      !is.na(gene_data_df$TE_log2FC), ]
  
  

  paired_t_test <- t.test(paired_data_clean$Gene_log2FC, 
                          paired_data_clean$TE_log2FC, 
                          paired = TRUE)
  
  cor_test_pearson <- cor.test(paired_data_clean$Gene_log2FC, 
                               paired_data_clean$TE_log2FC, 
                               method = "pearson")
  
  cor_test_spearman <- cor.test(paired_data_clean$Gene_log2FC, 
                                paired_data_clean$TE_log2FC, 
                                method = "spearman")
  
  te_clean <- gsub("_intergenic$", "", te_name)
  te_clean <- gsub("-", "", te_clean)
  
  if(sort_by_correlation) {
    lm_model <- lm(Gene_log2FC ~ TE_log2FC, data = paired_data_clean)
    predicted_gene <- predict(lm_model, newdata = paired_data_clean)
    residuals <- paired_data_clean$Gene_log2FC - predicted_gene
    
    paired_data_ordered <- paired_data_clean[order(residuals), ]
    subtitle_text <- paste("rho=", round(cor_test_spearman$estimate, 3), ", p=", format.pval(cor_test_spearman$p.value, digits = 3))
  } else {
    paired_data_ordered <- paired_data_clean[order(paired_data_clean$TE_log2FC), ]
    subtitle_text <- paste("Ordered by", te_clean, "TE log2FoldChange")
  }
  
  paired_data_ordered$Index <- 1:nrow(paired_data_ordered)
  
  plot_data_long <- data.frame(
    Index = rep(paired_data_ordered$Index, 2),
    Dataset = rep(paired_data_ordered$Dataset, 2),
    log2FC = c(paired_data_ordered$TE_log2FC, paired_data_ordered$Gene_log2FC),
    Type = rep(c(paste0(te_clean, "_TE"), paste0(gene_name, "_Gene")), 
               each = nrow(paired_data_ordered))
  )
  
  linear_trend_plot <- ggplot(plot_data_long, aes(x = Index, y = log2FC, color = Type)) +
    geom_point(size = 0.5, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_color_manual(values = setNames(c("#c98d22", "#756fae"), 
                                         c(paste0(te_clean, "_TE"), paste0(gene_name, "_Gene")))) +
    coord_cartesian(ylim = c(-max(abs(range(plot_data_long$log2FC, na.rm = TRUE))), 
                             max(abs(range(plot_data_long$log2FC, na.rm = TRUE))))) +
    labs(
      title = paste(gene_name, "Gene vs", te_clean, "TE Linear Trends"),
      subtitle = subtitle_text,
      x = "Sample Index (ordered by correlation strength)",
      y = "log2FoldChange",
      color = "Gene/TE Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 7, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 6),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6),
      legend.position = "none",
      axis.line = element_line(color = "black", size = 0.5),
      axis.line.x = element_line(color = "black", size = 0.5),
      axis.line.y = element_line(color = "black", size = 0.5),
      panel.grid.major = element_line(color = "gray90", size = 0.2),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 5)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

  paired_data_clean$Difference <- paired_data_clean$Gene_log2FC - paired_data_clean$TE_log2FC
  paired_data_clean$Index <- 1:nrow(paired_data_clean)
  
  paired_plot_data <- data.frame(
    Index = rep(1:nrow(paired_data_clean), 2),
    Dataset = rep(paired_data_clean$Dataset, 2),
    log2FC = c(paired_data_clean$Gene_log2FC, paired_data_clean$TE_log2FC),
    Type = rep(c(paste0(gene_name, "_Gene"), paste0(te_clean, "_TE")), 
               each = nrow(paired_data_clean))
  )
  print(paired_data_clean)
  paired_line_plot <- ggplot(paired_plot_data, aes(x = Type, y = log2FC, group = Index)) +
    geom_line(alpha = 0.8, color = "gray60") +
    geom_point(aes(color = Type), size = 0.5, alpha = 0.8) +
    scale_color_manual(values = setNames(c("#756fae", "#c98d22"), 
                                         c(paste0(gene_name, "_Gene"), paste0(te_clean, "_TE")))) +
    scale_x_discrete(limits = c(paste0(te_clean, "_TE"), paste0(gene_name, "_Gene")),
                     labels = c(te_clean, gene_name)) +
    coord_cartesian(ylim = c(-max(abs(range(paired_plot_data$log2FC, na.rm = TRUE))), 
                             max(abs(range(paired_plot_data$log2FC, na.rm = TRUE))))) +
    labs(
      title = paste("Paired Expression Comparison:", gene_name, "vs", te_clean),
      subtitle = paste("Mean difference =", round(paired_t_test$estimate, 3), 
                       ", p =", format.pval(paired_t_test$p.value, digits = 4)),
      x = "",
      y = "log2FoldChange",
      caption = "Lines connect paired observations from same dataset"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 7, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 6),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6),
      legend.position = "none",
      axis.line = element_line(color = "black", size = 0.5),
      axis.line.x = element_line(color = "black", size = 0.5),
      axis.line.y = element_line(color = "black", size = 0.5),
      panel.grid.major = element_line(color = "gray90", size = 0.2),
      panel.grid.minor = element_blank()
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)
  
  if(!is.null(output_dir)) {
    if(!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    file_prefix <- paste0(gene_name, "_vs_", te_clean)
    
    ggsave(filename = file.path(output_dir, paste0(file_prefix, "_trend_logfc_plot.pdf")), 
           plot = linear_trend_plot, width = 2, height = 2, dpi = 300)
    
    ggsave(filename = file.path(output_dir, paste0(file_prefix, "_paired_logfc_plot.pdf")), 
           plot = paired_line_plot, width = 1.5, height = 2, dpi = 300)
    
  }
  
  results <- list(
    gene_name = gene_name,
    te_name = te_name,
    sample_size = nrow(paired_data_clean),
    paired_t_test = paired_t_test,
    paired_data = paired_data_clean,
    plots = list(
      trend_plot = linear_trend_plot,
      paired_plot = paired_line_plot
    )
  )
  
  return(invisible(results))
}

analyze_gene_te_pair("ZNF182", "AluYk11_intergenic", 
                     output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/",sort_by_correlation = TRUE)
analyze_gene_te_pair("ZFP28", "L1PBa_intergenic", 
                     output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/",sort_by_correlation = TRUE)
analyze_gene_te_pair("ZNF287", "L1PBa_intergenic", 
                     output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/",sort_by_correlation = TRUE)
analyze_gene_te_pair("ZNF418", "LTR54_intergenic", 
                     output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/",sort_by_correlation = TRUE)

analyze_gene_te_pair("ZNF287", "L1MA10_intergenic", 
                     output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/",sort_by_correlation = TRUE)
analyze_gene_te_pair("ZNF287", "L1P4d_intergenic", 
                     output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/",sort_by_correlation = TRUE)
analyze_gene_te_pair("ZNF287", "L1PBb_intergenic", 
                     output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/",sort_by_correlation = TRUE)
analyze_gene_te_pair("ZNF432", "L1PBb_intergenic", 
                     output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/",sort_by_correlation = TRUE)
analyze_gene_te_pair("ZNF432", "L1P4d_intergenic", 
                     output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/",sort_by_correlation = TRUE)
analyze_gene_te_pair("ZNF432", "L1MA10_intergenic", 
                     output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/",sort_by_correlation = TRUE)
analyze_gene_te_pair("ZNF432", "L1PBa_intergenic", 
                     output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/",sort_by_correlation = TRUE)
analyze_gene_te_pair("ZFP28", "L1MA10_intergenic", 
                     output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/",sort_by_correlation = TRUE)
analyze_gene_te_pair("ZFP28", "L1PBb_intergenic", 
                     output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/",sort_by_correlation = TRUE)
analyze_gene_te_pair("ZFP28", "L1P4d_intergenic", 
                     output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/",sort_by_correlation = TRUE)
analyze_gene_te_pair("ZNF287", "L1PBa_intergenic", 
                     output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/",sort_by_correlation = TRUE)

analyze_gene_te_pair("IRF1", "UCON33_intergenic", 
                     output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_subfamily_Correlation/",sort_by_correlation = TRUE)
analyze_gene_te_pair("STAT2", "L1MA10_intergenic", 
                     output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_subfamily_Correlation/",sort_by_correlation = TRUE)
analyze_gene_te_pair("NFAT5", "MER61F_intergenic", 
                     output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_subfamily_Correlation/",sort_by_correlation = TRUE)

analyze_gene_te_pair("ZNF287", "L1MA10_intergenic", 
                                output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction")

analyze_gene_te_pair("ZNF287", "L1PBa_intergenic", 
                                output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction")




