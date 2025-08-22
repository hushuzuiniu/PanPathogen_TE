library(data.table)
library(DESeq2)
library(ggplot2)
library(dplyr)


raw_count_Gene <- fread("/data2t_2/hushu/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt")
# raw_count_Gene <- fread("/data2t_2/hushu/02.DESeq2_analysis_Gene/new_add_raw_data/new_add_sample_readscounts_matrix_combined_n90_Gene.txt")
raw_count_Gene$GeneName <- gene_id_to_name[raw_count_Gene$GeneID]
raw_count_TEs <- fread("/data2t_2/hushu/02.DESeq2_analysis_TE_loci/raw_data/Bacteria_SalTy_4_readscounts_matrix_TE_Loci.txt")
# raw_count_TEs <- fread("/data2t_2/hushu/02.DESeq2_analysis_TE_subfamily_feature_v2/raw_data/new_add_all_sample_readscounts_matrix_combined_n102_TE_subfamily_feature_v2.txt")
all_sample<-read_excel("/data2t_2/pathogen_TE_2025_New/01.new_raw_data/all_sample.xlsx")
Bacteria_SalTy_4<-all_sample[all_sample$Dataset=="Bacteria_SalTy_4",]

gene_columns <- colnames(raw_count_Gene)
te_columns <- colnames(raw_count_TEs)
bacteria_gene_cols <- gene_columns[gene_columns %in% Bacteria_SalTy_4$Sample_Name]
bacteria_te_cols <- te_columns[te_columns %in% Bacteria_SalTy_4$Sample_Name]

# bacteria_gene_cols <- grep("Bacteria_SalTy_4", gene_columns, value = TRUE)
# bacteria_te_cols <- grep("Bacteria_SalTy_4", te_columns, value = TRUE)
common_bacteria_cols <- intersect(bacteria_gene_cols, bacteria_te_cols)

pre_samples <- grep("_Pre_", common_bacteria_cols, value = TRUE)
post_samples <- grep("_Post_", common_bacteria_cols, value = TRUE)

gene_bacteria_data <- raw_count_Gene[, c("GeneID", "GeneName", common_bacteria_cols), with = FALSE]
te_bacteria_data <- raw_count_TEs[, c("Geneid", common_bacteria_cols), with = FALSE]

gene_count_matrix <- as.matrix(gene_bacteria_data[, -c("GeneID", "GeneName"), with = FALSE])
rownames(gene_count_matrix) <- gene_bacteria_data$GeneID

gene_sample_info <- data.frame(
  sample_id = colnames(gene_count_matrix),
  condition = ifelse(grepl("_Pre_", colnames(gene_count_matrix)), "Pre", "Post"),
  row.names = colnames(gene_count_matrix)
)

gene_dds <- DESeqDataSetFromMatrix(
  countData = gene_count_matrix,
  colData = gene_sample_info,
  design = ~ condition
)

keep_genes <- rowSums(counts(gene_dds) >= 10) >= ncol(gene_dds) * 0.25
gene_dds_filtered <- gene_dds[keep_genes, ]

gene_dds_filtered <- estimateSizeFactors(gene_dds_filtered)
gene_normalized_counts <- counts(gene_dds_filtered, normalized = TRUE)

te_count_matrix <- as.matrix(te_bacteria_data[, -"Geneid", with = FALSE])
rownames(te_count_matrix) <- te_bacteria_data$Geneid

te_sample_info <- data.frame(
  sample_id = colnames(te_count_matrix),
  condition = ifelse(grepl("_Pre_", colnames(te_count_matrix)), "Pre", "Post"),
  row.names = colnames(te_count_matrix)
)

te_dds <- DESeqDataSetFromMatrix(
  countData = te_count_matrix,
  colData = te_sample_info,
  design = ~ condition
)

# keep_tes <- rowSums(counts(te_dds) >= 10) >= ncol(te_dds) * 0.25
keep_tes <- rep(TRUE, nrow(counts(te_dds)))
te_dds_filtered <- te_dds[keep_tes, ]

te_dds_filtered <- estimateSizeFactors(te_dds_filtered)
te_normalized_counts <- counts(te_dds_filtered, normalized = TRUE)

NFKB1_id <- "ENSG00000109320.14"
NFKB1_present <- NFKB1_id %in% rownames(gene_normalized_counts)

AluYg6_c1738062_id <- "AluYg6_c1738062"
AluYg6_c1738062_present <- AluYg6_c1738062_id %in% rownames(te_normalized_counts)

NFKB1_normalized <- gene_normalized_counts[NFKB1_id, ]
AluYg6_c1738062_normalized <- te_normalized_counts[AluYg6_c1738062_id, ]

# log2
NFKB1_log2 <- log2(NFKB1_normalized + 1)
AluYg6_c1738062_log2 <- log2(AluYg6_c1738062_normalized + 1)

NFKB1_scaled <- scale(NFKB1_log2)[, 1]
AluYg6_c1738062_scaled <- scale(AluYg6_c1738062_log2)[, 1]

sample_names <- names(NFKB1_scaled)
condition_labels <- ifelse(grepl("_Pre_", sample_names), "Pre", "Post")

correlation_data <- data.frame(
  Sample = sample_names,
  Condition = condition_labels,
  NFKB1_raw = NFKB1_normalized[sample_names],
  AluYg6_c1738062_raw = AluYg6_c1738062_normalized[sample_names],
  NFKB1_log2 = NFKB1_log2[sample_names],
  AluYg6_c1738062_log2 = AluYg6_c1738062_log2[sample_names],
  NFKB1_scaled = NFKB1_scaled[sample_names],
  AluYg6_c1738062_scaled = AluYg6_c1738062_scaled[sample_names]
)

correlation_data_clean <- correlation_data[complete.cases(correlation_data), ]

gene_ranks <- rank(correlation_data_clean$NFKB1_scaled)
te_ranks <- rank(correlation_data_clean$AluYg6_c1738062_scaled)
rank_diff <- gene_ranks - te_ranks
n <- length(gene_ranks)
rho_overall <- 1 - (6 * sum(rank_diff^2)) / (n * (n^2 - 1))

spearman_overall <- cor.test(correlation_data_clean$NFKB1_scaled, 
                             correlation_data_clean$AluYg6_c1738062_scaled, 
                             method = "spearman", exact = FALSE)

pearson_overall <- cor.test(correlation_data_clean$NFKB1_scaled, 
                            correlation_data_clean$AluYg6_c1738062_scaled, 
                            method = "pearson")

pre_data <- correlation_data_clean[correlation_data_clean$Condition == "Pre", ]
post_data <- correlation_data_clean[correlation_data_clean$Condition == "Post", ]

if(nrow(pre_data) >= 3) {
  pre_spearman <- cor.test(pre_data$NFKB1_scaled, pre_data$AluYg6_c1738062_scaled, method = "spearman")
}

if(nrow(post_data) >= 3) {
  post_spearman <- cor.test(post_data$NFKB1_scaled, post_data$AluYg6_c1738062_scaled, method = "spearman")
}

colors <- c("Pre" = "#4c9bcf", "Post" = "#d00732")  
scatter_plot <- ggplot(correlation_data_clean, aes(x = NFKB1_scaled, y = AluYg6_c1738062_scaled, color = Condition)) +
  geom_point(size = 0.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid") +
  scale_color_manual(values = colors) +
  labs(
    x = "NFKB1 (scaled log2 expression)",
    y = "AluYg6_c1738062 (scaled log2 expression)",
    title = "NFKB1 vs AluYg6_c1738062 Expression Correlation\n(Bacteria_SalTy_4 samples, DESeq2 normalized)",
    color = "Condition"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 7, face = "bold"),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    legend.position = "none"
  ) +
  annotate("text", 
           x = max(correlation_data_clean$NFKB1_scaled) * 0.6, 
           y = max(correlation_data_clean$AluYg6_c1738062_scaled) * 0.9,
           label = paste("rho =", round(rho_overall, 3), ", p =", format(spearman_overall$p.value, digits = 3)),
           size = 2, fontface = "bold")

print(scatter_plot)
ggsave("/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_TE_Loci_Correlation/Bacteria_SalTy_4_NFKB1_vs_AluYg6_c1738062_expression_correlation.pdf",scatter_plot,width = 2,height = 2)
