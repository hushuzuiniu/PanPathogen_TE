Recurrent_Gene<-read_csv("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_Gene/plots/All_species_split_DE-Gene_padj0.05log2FC1.csv")
TEs_TF<-read_csv("/data2t_2/pathogen_TE_2025_New/08.TF_prediction/Fig3_TE_subfamily_all_fimo_results_coverage_percent_gt20_heatmap_unsort_no_pvalue.csv")
TE_loci_TF<-read_csv("/data2t_2/pathogen_TE_2025_New/08.TF_prediction/Fig3_TE_loci_fimo_results_heatmap_unsort_no_pvalue.csv")

non_ZNF_TEs_TF <- TEs_TF[!TEs_TF$TF_type %in% c("KRAB", "C2H2 ZF"), ]
non_ZNF_TEs_TF<-unique(non_ZNF_TEs_TF$motif_alt_id)
non_ZNF_TE_loci_TF <- TE_loci_TF[!TE_loci_TF$TF_type %in% c("KRAB", "C2H2 ZF"), ]
non_ZNF_TE_loci_TF<-unique(non_ZNF_TE_loci_TF$motif_alt_id)

non_ZNF <- unique(c(non_ZNF_TE_loci_TF, non_ZNF_TEs_TF))
Recurrent_non_ZNF <- Recurrent_Gene[
  tolower(Recurrent_Gene$GeneName) %in% tolower(non_ZNF),]
write.csv(Recurrent_non_ZNF,file = "/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_Gene/plots/non_ZNF_TF_All_split_DE-Gene_padj0.05log2FC1.csv",row.names = F) 
################################################################################################################################################
library(ggplot2)
library(reshape2)
library(patchwork)
library(scales)
library(dplyr)
setwd("/data2t_2/pathogen_TE_2025_New/08.TF_prediction/")
load("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_loci/plots/All_celltypes_combined_log2FoldChange_sig_DE-TEs.rdata")
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

########################################################################
NFKB1_Gene_fc<-All_Gene_fc[All_Gene_fc$GeneName=="NFKB1",]
THE1B_c4799565_fc<-All_TEs_fc[All_TEs_fc$GeneID=="THE1B_c4799565",]

gene_values <- as.numeric(NFKB1_Gene_fc[1, -c(1, ncol(NFKB1_Gene_fc)), with = FALSE])
te_values <- as.numeric(THE1B_c4799565_fc[1, -1])

dataset_names <- colnames(NFKB1_Gene_fc)[-c(1, ncol(NFKB1_Gene_fc))]  
plot_data <- data.frame(
  Dataset = rep(dataset_names, 2),
  log2FC = c(gene_values, te_values),
  Type = rep(c("NFKB1_Gene", "THE1B_c4799565_TE"), each = length(dataset_names))
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
  Type = rep(c("THE1B_c4799565_TE", "NFKB1_Gene"), each = nrow(paired_data))
)

te_lm <- lm(TE_log2FC ~ Index, data = paired_data)
gene_lm <- lm(Gene_log2FC ~ Index, data = paired_data)

gene_data <- plot_data[plot_data$Type == "NFKB1_Gene", ]
te_data <- plot_data[plot_data$Type == "THE1B_c4799565_TE", ]

paired_data <- merge(gene_data[, c("Dataset", "log2FC")], 
                     te_data[, c("Dataset", "log2FC")], 
                     by = "Dataset", 
                     suffixes = c("_Gene", "_TE"))

paired_data_clean <- paired_data[!is.na(paired_data$log2FC_Gene) & 
                                   !is.na(paired_data$log2FC_TE), ]

paired_t_test <- t.test(paired_data_clean$log2FC_Gene, 
                        paired_data_clean$log2FC_TE, 
                        paired = TRUE)


paired_data_clean$Difference <- paired_data_clean$log2FC_Gene - paired_data_clean$log2FC_TE
paired_data_clean$Index <- 1:nrow(paired_data_clean)

paired_plot_data <- data.frame(
  Index = rep(1:nrow(paired_data_clean), 2),
  Dataset = rep(paired_data_clean$Dataset, 2),
  log2FC = c(paired_data_clean$log2FC_Gene, paired_data_clean$log2FC_TE),
  Type = rep(c("NFKB1_Gene", "THE1B_c4799565_TE"), each = nrow(paired_data_clean))
)

paired_line_plot <- ggplot(paired_plot_data, aes(x = Type, y = log2FC, group = Index)) +
  geom_line(alpha = 0.8, color = "gray60") +
  geom_point(aes(color = Type), size = 0.5, alpha = 0.8) +
  scale_color_manual(values = setNames(c("#756fae", "#c98d22"), 
                                       c(paste0("NFKB1_Gene"), paste0("THE1B_c4799565_TE")))) +
  scale_x_discrete(limits = c("THE1B_c4799565_TE", "NFKB1_Gene"),
                   labels = c("THE1B_c4799565", "NFKB1")) + 
  coord_cartesian(ylim = c(-max(abs(range(paired_plot_data$log2FC, na.rm = TRUE))), 
                           max(abs(range(paired_plot_data$log2FC, na.rm = TRUE))))) +
  labs(
    title = paste("Paired Expression Comparison NFKB1 vs THE1B_c4799565"),
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
ggsave("/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/NFKB1_vs_THE1B_c4799565_paired_logfc_plot.pdf", paired_line_plot, width = 2, height = 2)

cor_test_spearman <- cor.test(paired_data_clean$log2FC_Gene, 
                              paired_data_clean$log2FC_TE, 
                              method = "spearman")
cor_test_spearman

scatter_plot <- ggplot(paired_data_clean, aes(x = log2FC_TE, y = log2FC_Gene)) +
  geom_point(size = 0.5, alpha = 0.8, color = "#2c3e50") +
  geom_smooth(method = "lm", se = TRUE, color = "#e74c3c") +
  labs(
    title = paste("NFKB1 Gene vs THE1B_c4799565 log2FoldChange Correlation"),
    subtitle = paste("rho =", round(cor_test_spearman$estimate, 3), 
                     ", p =", format(cor_test_spearman$p.value, scientific = TRUE)),
    x = "THE1B_c4799565 TE log2FoldChange",
    y = "NFKB1 Gene log2FoldChange"
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
  ) 
scatter_plot
ggsave("/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/NFKB1_vs_THE1B_c4799565_scatter_plot.pdf", scatter_plot, width = 2, height = 2)

paired_data_clean$mean_expression <- (paired_data_clean$log2FC_Gene + paired_data_clean$log2FC_TE) / 2
correlation_ordered <- paired_data_clean[order(paired_data_clean$mean_expression), ]
correlation_ordered$Index <- 1:nrow(correlation_ordered)

plot_data_corr <- data.frame(
  Index = rep(correlation_ordered$Index, 2),
  log2FC = c(correlation_ordered$log2FC_Gene, correlation_ordered$log2FC_TE),
  Type = rep(c("NFKB1_Gene", "THE1B_c4799565_TE"), each = nrow(correlation_ordered)),
  Dataset = rep(correlation_ordered$Dataset, 2)
)

correlation_plot <- ggplot(plot_data_corr, aes(x = Index, y = log2FC, color = Type)) +
  geom_point(size = 0.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = c("THE1B_c4799565_TE" = "#756fae", "NFKB1_Gene" = "#c98d22")) +
  labs(
    title = "NFKB1 Gene vs THE1B_c4799565 TE - Ordered by Combined Expression",
    subtitle = paste("rho =", round(cor_test_spearman$estimate, 3), 
                     ", p =", format(cor_test_spearman$p.value, scientific = TRUE)),
    x = "Dataset",
    y = "log2FoldChange"
  ) +
  theme_minimal()+
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
  ) 
correlation_plot
ggsave("/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/NFKB1_vs_THE1B_c4799565_correlation_plot.pdf", correlation_plot, width = 2, height = 2)
########################################################################################################################
analyze_gene_te_correlation <- function(gene_name, te_name, gene_data = All_Gene_fc, te_data = All_TEs_fc,output_dir = NULL,save_plots = TRUE) {
  gene_fc <- gene_data[gene_data$GeneName == gene_name, ]
  te_fc <- te_data[te_data$GeneID == te_name, ]
  if(nrow(gene_fc) == 0) {stop(gene_name)}; if(nrow(te_fc) == 0) {stop(te_name)}

  gene_values <- as.numeric(gene_fc[1, -c(1, ncol(gene_fc)), with = FALSE])
  te_values <- as.numeric(te_fc[1, -1])
  
  dataset_names <- colnames(gene_fc)[-c(1, ncol(gene_fc))]
  plot_data <- data.frame(
    Dataset = rep(dataset_names, 2),
    log2FC = c(gene_values, te_values),
    Type = rep(c(paste0(gene_name, "_Gene"), paste0(te_name, "_TE")), each = length(dataset_names)))
    plot_data_clean <- plot_data[!is.na(plot_data$log2FC), ]
  
  gene_data_df <- data.frame(
    Dataset = dataset_names,
    Gene_log2FC = gene_values,
    TE_log2FC = te_values
  )

  gene_plot_data <- plot_data[plot_data$Type == paste0(gene_name, "_Gene"), ]
  te_plot_data <- plot_data[plot_data$Type == paste0(te_name, "_TE"), ]
  
  paired_data <- merge(gene_plot_data[, c("Dataset", "log2FC")], 
                       te_plot_data[, c("Dataset", "log2FC")], 
                       by = "Dataset", 
                       suffixes = c("_Gene", "_TE"))
  
  paired_data_clean <- paired_data[!is.na(paired_data$log2FC_Gene) & 
                                     !is.na(paired_data$log2FC_TE), ]
  print(paired_data_clean)
  
  if(nrow(paired_data_clean) < 3) {
    warning("(n < 3)")}

  paired_t_test <- t.test(paired_data_clean$log2FC_Gene, 
                          paired_data_clean$log2FC_TE, 
                          paired = TRUE)
  
  cor_test_spearman <- cor.test(paired_data_clean$log2FC_Gene, 
                                paired_data_clean$log2FC_TE, 
                                method = "spearman")
  

  te_clean <- gsub("_intergenic$", "", te_name)
  te_clean <- gsub("-", "", te_clean)
  
 
  paired_data_clean$Difference <- paired_data_clean$log2FC_Gene - paired_data_clean$log2FC_TE
  paired_data_clean$Index <- 1:nrow(paired_data_clean)
  
  paired_plot_data <- data.frame(
    Index = rep(1:nrow(paired_data_clean), 2),
    Dataset = rep(paired_data_clean$Dataset, 2),
    log2FC = c(paired_data_clean$log2FC_Gene, paired_data_clean$log2FC_TE),
    Type = rep(c(paste0(gene_name, "_Gene"), paste0(te_clean, "_TE")), 
               each = nrow(paired_data_clean))
  )
  
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
      title = paste("Paired Expression Comparison", gene_name, "vs", te_clean),
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
  
  scatter_plot <- ggplot(paired_data_clean, aes(x = log2FC_TE, y = log2FC_Gene)) +
    geom_point(size = 0.5, alpha = 0.8, color = "#2c3e50") +
    geom_smooth(method = "lm", se = TRUE, color = "#e74c3c") +
    labs(
      title = paste(gene_name, "Gene vs", te_clean, "log2FoldChange Correlation"),
      subtitle = paste("rho =", round(cor_test_spearman$estimate, 3), 
                       ", p =", format(cor_test_spearman$p.value, scientific = TRUE)),
      x = paste(te_clean, "TE log2FoldChange"),
      y = paste(gene_name, "Gene log2FoldChange")
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
    )
  

  paired_data_clean$mean_expression <- (paired_data_clean$log2FC_Gene + paired_data_clean$log2FC_TE) / 2
  correlation_ordered <- paired_data_clean[order(paired_data_clean$mean_expression), ]
  correlation_ordered$Index <- 1:nrow(correlation_ordered)
  
  plot_data_corr <- data.frame(
    Index = rep(correlation_ordered$Index, 2),
    log2FC = c(correlation_ordered$log2FC_Gene, correlation_ordered$log2FC_TE),
    Type = rep(c(paste0(gene_name, "_Gene"), paste0(te_clean, "_TE")), 
               each = nrow(correlation_ordered)),
    Dataset = rep(correlation_ordered$Dataset, 2)
  )
  
  correlation_plot <- ggplot(plot_data_corr, aes(x = Index, y = log2FC, color = Type)) +
    geom_point(size = 0.5, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_color_manual(values = setNames(c("#c98d22", "#756fae"), 
                                         c(paste0(te_clean, "_TE"), paste0(gene_name, "_Gene")))) +
    labs(
      title = paste(gene_name, "Gene vs", te_clean, "TE - Ordered by Combined Expression"),
      subtitle = paste("rho =", round(cor_test_spearman$estimate, 3), 
                       ", p =", format(cor_test_spearman$p.value, scientific = TRUE)),
      x = "Dataset",
      y = "log2FoldChange"
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
    )
  
  if(save_plots && !is.null(output_dir)) {
    if(!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    file_prefix <- paste0(gene_name, "_vs_", te_clean)
    
    ggsave(filename = file.path(output_dir, paste0(file_prefix, "_paired_logfc_plot.pdf")), 
           plot = paired_line_plot, width = 1.5, height = 2, dpi = 300)
    
    ggsave(filename = file.path(output_dir, paste0(file_prefix, "_scatter_plot.pdf")), 
           plot = scatter_plot, width = 2, height = 2, dpi = 300)
    
    ggsave(filename = file.path(output_dir, paste0(file_prefix, "_trend_plot.pdf")), 
           plot = correlation_plot, width = 2, height = 2, dpi = 300)
    
  }
  
  results <- list(
    gene_name = gene_name,
    te_name = te_name,
    sample_size = nrow(paired_data_clean),
    paired_t_test = paired_t_test,
    correlation_test = cor_test_spearman,
    paired_data = paired_data_clean,
    plots = list(
      paired_plot = paired_line_plot,
      scatter_plot = scatter_plot,
      correlation_plot = correlation_plot
    )
  )
  
  return(invisible(results))
}
analyze_gene_te_correlation("NFKB2", "AluSx1_c5104405",
                            output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/")
analyze_gene_te_correlation("HES1", "AluYg6_c1738062",
                            output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/")
analyze_gene_te_correlation("NFKB2", "AluYg6_c1738062",
                            output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/")
analyze_gene_te_correlation("NFKB1", "AluYg6_c1738062",
                            output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/")


analyze_gene_te_correlation("NFKB1", "THE1B_c4799565",
                            output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/")
analyze_gene_te_correlation("NFKB1", "AluYg6_c2536804",
                            output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/")
analyze_gene_te_correlation("NFKB1", "AluYg6_c3729459",
                            output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/")
analyze_gene_te_correlation("NFKB2", "AluYg6_c2536804",
                            output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/")
analyze_gene_te_correlation("NFKB2", "AluYg6_c3729459",
                            output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/")
analyze_gene_te_correlation("NFKB2", "AluSx1_c5104399",
                            output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/")
analyze_gene_te_correlation("NFKB2", "THE1B_c4799565",
                            output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/")

analyze_gene_te_correlation("HEY1", "AluYg6_c1738062",
                            output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/")
analyze_gene_te_correlation("HEY1", "AluYg6_c2536804",
                            output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/")
analyze_gene_te_correlation("HEY1", "AluYg6_c3729459",
                            output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/")
analyze_gene_te_correlation("HEY1", "L1PA7_c2329676",
                            output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/")

analyze_gene_te_correlation("HES1", "AluYg6_c2536804",
                            output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/")
analyze_gene_te_correlation("HES1", "AluYg6_c3729459",
                            output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/")
analyze_gene_te_correlation("HES1", "L1PA7_c2329676",
                            output_dir = "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/")






