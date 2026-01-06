library(dorothea)
library(viper)
library(dplyr)
library(data.table)
library(readxl)

setwd("/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_prediction_add_bg_v3/TF_NES_vs_ruTE")

# ==================== Step 1: 获取NFKB1和RELA的regulon ====================
cat("Step 1: 获取TF regulon...\n")
regulons <- dorothea_hs %>% filter(confidence %in% c("A", "B", "C"))
cat("A+B+C级别regulon总数:", nrow(regulons), "\n")
cat("NFKB1 targets:", sum(regulons$tf == "NFKB1"), "\n")
cat("RELA targets:", sum(regulons$tf == "RELA"), "\n")

# 转换为viper格式的regulon
regulon_viper <- dorothea::df2regulon(regulons)

# ==================== Step 2: 读取并合并基因表达矩阵 ====================
cat("\nStep 2: 读取基因表达矩阵...\n")

raw_count_Gene_1 <- fread("/data2t_2/hushu/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt")
raw_count_Gene_2 <- fread("/data2t_2/hushu/02.DESeq2_analysis_Gene/new_add_raw_data/new_add_sample_readscounts_matrix_combined_n90_Gene.txt")

cat("矩阵1:", ncol(raw_count_Gene_1) - 1, "样本\n")
cat("矩阵2:", ncol(raw_count_Gene_2) - 1, "样本\n")

# 合并两个矩阵
raw_count_Gene <- merge(raw_count_Gene_1, raw_count_Gene_2, by = "GeneID", all = TRUE)
for (col in names(raw_count_Gene)[-1]) {
  raw_count_Gene[[col]][is.na(raw_count_Gene[[col]])] <- 0
}

cat("合并后:", nrow(raw_count_Gene), "基因 x", ncol(raw_count_Gene) - 1, "样本\n")

# ==================== Step 3: ENSG ID转换为Gene Symbol ====================
cat("\nStep 3: 从GTF文件提取ID映射...\n")

gtf_file <- "/data2t_2/pathogen_TE_2025_New/00.ref/hg38.p13.gene.anno.gtf"
gtf_lines <- fread(cmd = paste0("grep -v '^#' ", gtf_file, " | awk -F'\\t' '$3==\"gene\"'"),
                   header = FALSE, sep = "\t")

extract_attr <- function(attr_string, attr_name) {
  pattern <- paste0(attr_name, ' "([^"]+)"')
  m <- regmatches(attr_string, regexec(pattern, attr_string))
  sapply(m, function(x) if(length(x) > 1) x[2] else NA)
}

gtf_gene_id <- extract_attr(gtf_lines$V9, "gene_id")
gtf_gene_name <- extract_attr(gtf_lines$V9, "gene_name")

id_mapping <- data.frame(
  gene_id = gtf_gene_id,
  gene_name = gtf_gene_name,
  stringsAsFactors = FALSE
)
id_mapping <- id_mapping[!duplicated(id_mapping$gene_id), ]

cat("GTF中基因数:", nrow(id_mapping), "\n")

ensg_to_symbol <- setNames(id_mapping$gene_name, id_mapping$gene_id)
gene_symbols <- ensg_to_symbol[raw_count_Gene$GeneID]

cat("成功映射的基因数:", sum(!is.na(gene_symbols)), "/", length(gene_symbols), "\n")
raw_count_Gene$Symbol <- gene_symbols

# ==================== Step 4: 读取样本信息 ====================
cat("\nStep 4: 读取样本信息...\n")

all_sample <- read_excel("/data2t_2/pathogen_TE_2025_New/01.new_raw_data/all_sample.xlsx")
cat("all_sample中样本数:", nrow(all_sample), "\n")

# 找到在表达矩阵中的样本
gene_columns <- colnames(raw_count_Gene)[!colnames(raw_count_Gene) %in% c("GeneID", "Symbol")]
common_samples <- gene_columns[gene_columns %in% all_sample$Sample_Name]
cat("表达矩阵与all_sample共同样本数:", length(common_samples), "\n")

# ==================== Step 5: 准备表达矩阵 ====================
cat("\nStep 5: 准备表达矩阵...\n")
# 过滤有Symbol的行
gene_data <- raw_count_Gene[!is.na(raw_count_Gene$Symbol), c("Symbol", common_samples), with = FALSE]
# 去掉重复的Symbol（取平均值）
gene_data <- gene_data[, lapply(.SD, mean), by = Symbol, .SDcols = common_samples]
cat("去重后基因数:", nrow(gene_data), "\n")

# 转换为矩阵
count_matrix <- as.matrix(gene_data[, -"Symbol", with = FALSE])
rownames(count_matrix) <- gene_data$Symbol

cat("表达矩阵:", nrow(count_matrix), "基因 x", ncol(count_matrix), "样本\n")

# ==================== Step 6: log2转换 ====================
cat("\nStep 6: log2(counts+1)转换...\n")
# 直接log2转换
expr_matrix <- log2(count_matrix + 1)
cat("表达矩阵:", nrow(expr_matrix), "基因 x", ncol(expr_matrix), "样本\n")

# ==================== Step 7: VIPER计算TF活性 ====================
cat("\nStep 7: VIPER计算TF活性...\n")
tf_activity <- viper(expr_matrix, regulon_viper, method = "none", minsize = 5, verbose = FALSE)
cat("TF活性矩阵:", nrow(tf_activity), "TF x", ncol(tf_activity), "样本\n")

# 提取NFKB1和RELA
tfs_of_interest <- c("NFKB1", "RELA")
available_tfs <- tfs_of_interest[tfs_of_interest %in% rownames(tf_activity)]
cat("成功计算TF活性:", paste(available_tfs, collapse = ", "), "\n")
tf_activity_subset <- tf_activity[available_tfs, , drop = FALSE]
cat("最终矩阵:", nrow(tf_activity_subset), "TF x", ncol(tf_activity_subset), "样本\n")

# ==================== Step 8: 添加样本信息并保存 ====================
cat("\nStep 8: 保存结果...\n")

# 保存2×N矩阵格式（TF为行，样本为列）
matrix_file <- "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_prediction_add_bg_v3/TF_NES_vs_ruTE/TF_activity_matrix_2xN.csv"
write.csv(tf_activity_subset, matrix_file, row.names = TRUE)
cat("矩阵格式保存至:", matrix_file, "\n")

# 转换为数据框（带样本信息的长格式）
result_df <- as.data.frame(t(tf_activity_subset))
result_df$Sample <- rownames(result_df)
colnames(result_df) <- c("NFKB1_NES", "RELA_NES", "Sample")

# 添加样本信息
sample_info_df <- all_sample[all_sample$Sample_Name %in% result_df$Sample,
                              c("Sample_Name", "Dataset", "CellType_1", "Condition")]
colnames(sample_info_df)[1] <- "Sample"
result_df <- merge(result_df, sample_info_df, by = "Sample", all.x = TRUE)
# 重新排列列顺序
result_df <- result_df[, c("Sample", "Dataset", "CellType_1", "Condition", "NFKB1_NES", "RELA_NES")]
# 保存
output_file <- "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_prediction_add_bg_v3/TF_NES_vs_ruTE/TF_activity_all_samples.csv"
write.csv(result_df, output_file, row.names = FALSE)
