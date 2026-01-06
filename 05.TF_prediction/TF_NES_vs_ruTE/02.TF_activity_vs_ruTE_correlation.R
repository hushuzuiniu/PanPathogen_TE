library(data.table)
library(dplyr)
library(ggplot2)

setwd("/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_prediction_add_bg_v3/TF_NES_vs_ruTE")

# ==================== Step 1: 定义TF-TE对 ====================
cat("Step 1: 定义TF-TE对...\n")

tf_te_pairs <- list(
  c("RELA", "HERV3-int_intergenic"),
  c("NFKB1", "HERV9NC-int_intergenic"),
  c("RELA", "HERV9NC-int_intergenic"),
  c("RELA", "HERVIP10FH-int_intergenic"),
  c("RELA", "L1P4d_intergenic"),
  c("NFKB1", "L1PBa_intergenic"),
  c("RELA", "L1PBa_intergenic"),
  c("NFKB1", "L1PBb_intergenic"),
  c("RELA", "L1PBb_intergenic"),
  c("RELA", "LTR27E_intergenic"),
  c("RELA", "LTR54_intergenic"),
  c("RELA", "MER4-int_intergenic"),
  c("RELA", "THE1C-int_intergenic")
)

cat("共", length(tf_te_pairs), "个TF-TE对\n")

# ==================== Step 2: 读取TE表达矩阵 ====================
cat("\nStep 2: 读取TE表达矩阵...\n")

raw_count_TEs_1 <- fread("/data2t_2/hushu/02.DESeq2_analysis_TE_subfamily_feature_v2/raw_data/new_add_all_sample_readscounts_matrix_combined_n102_TE_subfamily_feature_v2.txt")
raw_count_TEs_2 <- fread("/data2t_2/hushu/02.DESeq2_analysis_TE_subfamily_feature_v2/raw_data/new_all_sample_readscounts_matrix_combined_n4231_TE_subfamily_feature_v2.txt")

cat("TE矩阵1: ", ncol(raw_count_TEs_1) - 1, "样本\n")
cat("TE矩阵2: ", ncol(raw_count_TEs_2) - 1, "样本\n")

raw_count_TEs <- merge(raw_count_TEs_1, raw_count_TEs_2, by = "GeneID", all = TRUE)
for (col in names(raw_count_TEs)[-1]) {
  raw_count_TEs[[col]][is.na(raw_count_TEs[[col]])] <- 0
}

cat("合并后: ", nrow(raw_count_TEs), "TE x ", ncol(raw_count_TEs) - 1, "样本\n")

# ==================== Step 3: 读取TF活性数据 ====================
cat("\nStep 3: 读取TF活性数据...\n")

tf_activity_file <- "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_prediction_add_bg_v3/TF_NES_vs_ruTE/TF_activity_all_samples.csv"
tf_all_samples <- fread(tf_activity_file)

cat("TF活性数据: ", nrow(tf_all_samples), "样本\n")

# ==================== Step 4: 读取TE上调的Dataset信息 ====================
cat("\nStep 4: 读取TE上调的Dataset信息...\n")

te_up_file <- "/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_subfamily_feature/plots/All_species_split_DE-TEs_feature_padj0.05log2FC1.csv"
te_up_df <- fread(te_up_file)

cat("TE上调信息: ", nrow(te_up_df), "TEs x", ncol(te_up_df) - 1, "datasets\n")

# 获取该TE上调的Dataset列表
# 定义已知的CellType前缀（按长度降序排列，优先匹配较长的）
known_celltypes <- c("Macrophages", "Dendritic", "Monocyte", "B_cell", "T_cell",
                      "A549", "PBMCs", "Liver")

# CellType 映射: te_up_df 列名前缀 -> tf_all_samples 中的 CellType_1 可能值
celltype_mapping <- list(
  "Macrophages" = c("primary macrophages", "MDM", "Fetal placental macrophages",
                    "stem cell-derived macrophage", "alveolar macrophages"),
  "Dendritic" = c("Dendritic_cell"),
  "Monocyte" = c("Monocyte"),
  "B_cell" = c("B_cell"),
  "T_cell" = c("T_cell", "CD4+_T_cell"),
  "A549" = c("A549"),
  "PBMCs" = c("PBMCs"),
  "Liver" = c("Hepatocyte", "Huh-7", "Huh-7.5.1")
)

# 从列名解析CellType和Dataset
parse_column_name <- function(col_name) {
  for (ct in known_celltypes) {
    if (startsWith(col_name, paste0(ct, "_"))) {
      dataset <- sub(paste0("^", ct, "_"), "", col_name)
      return(list(CellType = ct, Dataset = dataset))
    }
  }
  # 如果没有匹配到已知的CellType，用第一个_分割
  parts <- strsplit(col_name, "_")[[1]]
  return(list(CellType = parts[1], Dataset = paste(parts[-1], collapse = "_")))
}

get_te_up_datasets <- function(te_name, te_up_df) {
  if (!(te_name %in% te_up_df$GeneID)) {
    return(data.frame())
  }

  te_row <- te_up_df[te_up_df$GeneID == te_name, ]
  dataset_cols <- colnames(te_up_df)[-1]  # 除GeneID外的所有列

  up_info <- data.frame(ColName = character(), CellType = character(), Dataset = character(),
                        stringsAsFactors = FALSE)

  for (col in dataset_cols) {
    if (te_row[[col]] == "Up") {
      parsed <- parse_column_name(col)
      up_info <- rbind(up_info, data.frame(ColName = col,
                                            CellType = parsed$CellType,
                                            Dataset = parsed$Dataset,
                                            stringsAsFactors = FALSE))
    }
  }

  return(up_info)
}

# ==================== Step 5: 创建输出目录 ====================
output_dir <- "/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_prediction_add_bg_v3/TF_NES_vs_ruTE/TF_activity_vs_ruTE_correlation_UpOnly"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==================== Step 6: 对每个TF-TE对、每个Dataset分别计算相关性 ====================
cat("\nStep 6: 计算相关性 (按TF-TE对和Dataset分开)...\n")

all_correlation_results <- list()
per_dataset_results <- list()

for (pair in tf_te_pairs) {
  tf_name <- pair[1]
  te_name <- pair[2]
  te_clean <- gsub("-", "_", gsub("_intergenic$", "", te_name))

  cat("\n========================================\n")
  cat("处理:", tf_name, "vs", te_name, "\n")

  tf_col <- paste0(tf_name, "_NES")

  # 检查TF列是否存在
  if (!(tf_col %in% colnames(tf_all_samples))) {
    cat("警告: 未找到", tf_col, "列\n")
    next
  }

  # 检查TE是否在表达矩阵中
  if (!(te_name %in% raw_count_TEs$GeneID)) {
    cat("警告: 未找到TE:", te_name, "\n")
    next
  }

  # 获取该TE上调的Dataset列表 (返回data.frame: ColName, CellType, Dataset)
  up_info <- get_te_up_datasets(te_name, te_up_df)

  cat("TE上调的Dataset数:", nrow(up_info), "\n")

  if (nrow(up_info) == 0) {
    cat("警告: 该TE没有上调的Dataset\n")
    next
  }

  # 获取TE表达数据
  te_row <- raw_count_TEs[raw_count_TEs$GeneID == te_name, ]
  te_columns <- colnames(raw_count_TEs)[-1]

  # ==================== 对每个CellType-Dataset组合分别计算 ====================
  for (i in 1:nrow(up_info)) {
    col_name <- up_info$ColName[i]
    celltype <- up_info$CellType[i]
    dataset <- up_info$Dataset[i]

    # 获取对应的 tf_all_samples 中的 CellType_1 值
    mapped_celltypes <- celltype_mapping[[celltype]]
    if (is.null(mapped_celltypes)) {
      mapped_celltypes <- celltype  # 如果没有映射，使用原值
    }

    # 筛选该Dataset和CellType的样本
    tf_dataset <- tf_all_samples[tf_all_samples$Dataset == dataset &
                                  tf_all_samples$CellType_1 %in% mapped_celltypes, ]

    if (nrow(tf_dataset) == 0) {
      next
    }

    # 找到共同样本
    common_samples <- tf_dataset$Sample[tf_dataset$Sample %in% te_columns]

    if (length(common_samples) == 0) {
      next
    }

    # 构建合并数据表
    merged_data <- data.frame(Sample = common_samples, stringsAsFactors = FALSE)

    # 添加TF活性
    tf_subset <- tf_dataset[tf_dataset$Sample %in% common_samples, ]
    merged_data$TF_NES <- tf_subset[[tf_col]][match(merged_data$Sample, tf_subset$Sample)]
    merged_data$Condition <- tf_subset$Condition[match(merged_data$Sample, tf_subset$Sample)]
    merged_data$CellType <- tf_subset$CellType_1[match(merged_data$Sample, tf_subset$Sample)]

    # 添加TE表达量
    te_expr <- as.numeric(te_row[, merged_data$Sample, with = FALSE])
    merged_data$TE_raw <- te_expr
    merged_data$TE_log2 <- log2(merged_data$TE_raw + 1)

    # 移除NA
    merged_data <- merged_data[complete.cases(merged_data), ]


    n_pre <- sum(grepl("^[Pp]re", merged_data$Condition))
    n_post <- sum(grepl("^[Pp]ost", merged_data$Condition))

    # 计算Spearman相关性
    cor_result <- cor.test(merged_data$TF_NES, merged_data$TE_log2,
                            method = "spearman", exact = FALSE)

    rho <- cor_result$estimate
    p_value <- cor_result$p.value

    cat("  ", col_name, ": n=", nrow(merged_data), ", rho=", round(rho, 3), ", p=", format(p_value, digits = 2), "\n")

    # 保存结果 (使用col_name区分不同CellType-Dataset组合)
    result_key <- paste0(tf_name, "_", te_name, "_", col_name)
    per_dataset_results[[result_key]] <- data.frame(
      TF = tf_name,
      CellType = celltype,
      TE = te_name,
      Dataset = dataset,
      N_samples = nrow(merged_data),
      N_Pre = n_pre,
      N_Post = n_post,
      Spearman_rho = rho,
      P_value = p_value,
      Significant = ifelse(p_value < 0.05, "Yes", "No")
    )

    # 创建子目录
    pair_dir <- file.path(output_dir, paste0(tf_name, "_vs_", te_clean))
    dir.create(pair_dir, showWarnings = FALSE, recursive = TRUE)

    # 绘制散点图
    merged_data$ConditionGroup <- ifelse(grepl("^[Pp]re", merged_data$Condition), "Pre", "Post")
    colors <- c("Pre" = "#6c9dc0", "Post" = "#c1543d")

    p <- ggplot(merged_data, aes(x = TF_NES, y = TE_log2, color = ConditionGroup)) +
      geom_point(size = 1, alpha = 0.8) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid", linewidth = 0.5) +
      scale_color_manual(values = colors) +
      labs(
        x = paste0(tf_name, " Activity (NES)"),
        y = paste0(te_clean, " (log2(count+1))"),
        title = paste0(tf_name, " vs ", te_clean, "\n", col_name),
        color = "Condition"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 6, face = "bold"),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 5),
        legend.position = "none"
      ) +
      annotate("text",
               x = min(merged_data$TF_NES) + (max(merged_data$TF_NES) - min(merged_data$TF_NES)) * 0.5,
               y = max(merged_data$TE_log2) * 0.95,
               label = paste0("n=", nrow(merged_data), ", rho=", round(rho, 2), ", p=", format(p_value, digits = 2)),
               size = 1.8, fontface = "bold")

    ggsave(filename = file.path(pair_dir, paste0(col_name, "_scatter.pdf")),
           plot = p, width = 2, height = 2)
  }

  # ==================== 计算该TF-TE对的总体相关性 ====================
  tf_up_only <- tf_all_samples[tf_all_samples$Dataset %in% up_info$Dataset, ]
  common_samples <- tf_up_only$Sample[tf_up_only$Sample %in% te_columns]

  if (length(common_samples) >= 10) {
    merged_data <- data.frame(Sample = common_samples, stringsAsFactors = FALSE)
    tf_subset <- tf_up_only[tf_up_only$Sample %in% common_samples, ]
    merged_data$TF_NES <- tf_subset[[tf_col]][match(merged_data$Sample, tf_subset$Sample)]
    merged_data$Condition <- tf_subset$Condition[match(merged_data$Sample, tf_subset$Sample)]
    merged_data$Dataset <- tf_subset$Dataset[match(merged_data$Sample, tf_subset$Sample)]
    merged_data$CellType <- tf_subset$CellType_1[match(merged_data$Sample, tf_subset$Sample)]

    te_expr <- as.numeric(te_row[, merged_data$Sample, with = FALSE])
    merged_data$TE_raw <- te_expr
    merged_data$TE_log2 <- log2(merged_data$TE_raw + 1)
    merged_data <- merged_data[complete.cases(merged_data), ]

    n_pre <- sum(grepl("^[Pp]re", merged_data$Condition))
    n_post <- sum(grepl("^[Pp]ost", merged_data$Condition))

    cor_result <- cor.test(merged_data$TF_NES, merged_data$TE_log2,
                            method = "spearman", exact = FALSE)
    rho <- cor_result$estimate
    p_value <- cor_result$p.value

    cat("  总体: n=", nrow(merged_data), ", rho=", round(rho, 3), ", p=", format(p_value, digits = 2), "\n")

    all_correlation_results[[paste0(tf_name, "_", te_name)]] <- data.frame(
      TF = tf_name,
      TE = te_name,
      N_Up_Datasets = nrow(up_info),
      N_samples = nrow(merged_data),
      N_Pre = n_pre,
      N_Post = n_post,
      Spearman_rho = rho,
      P_value = p_value,
      Significant = ifelse(p_value < 0.05, "Yes", "No")
    )

    # 保存总的merged_data
    merged_data$TF <- tf_name
    merged_data$TE <- te_name
    merged_data <- merged_data[, c("Sample", "Dataset", "CellType", "Condition", "TF", "TF_NES", "TE", "TE_raw", "TE_log2")]
    pair_dir <- file.path(output_dir, paste0(tf_name, "_vs_", te_clean))
    write.csv(merged_data, file.path(pair_dir, paste0(tf_name, "_vs_", te_clean, ".csv")), row.names = FALSE)

    # 总体散点图
    merged_data$ConditionGroup <- ifelse(grepl("^[Pp]re", merged_data$Condition), "Pre", "Post")
    colors <- c("Pre" = "#6c9dc0", "Post" = "#c1543d")

    p <- ggplot(merged_data, aes(x = TF_NES, y = TE_log2, color = ConditionGroup)) +
      geom_point(size = 0.5, alpha = 0.8) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid", linewidth = 0.5) +
      scale_color_manual(values = colors) +
      labs(
        x = paste0(tf_name, " Activity (NES)"),
        y = paste0(te_clean, " (log2(count+1))"),
        title = paste0(tf_name, " vs ", te_clean, " (All)"),
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
               x = min(merged_data$TF_NES) + (max(merged_data$TF_NES) - min(merged_data$TF_NES)) * 0.5,
               y = max(merged_data$TE_log2) * 0.95,
               label = paste0("n=", nrow(merged_data), ", rho=", round(rho, 2), ", p=", format(p_value, digits = 2)),
               size = 2, fontface = "bold")

    pair_dir <- file.path(output_dir, paste0(tf_name, "_vs_", te_clean))
    ggsave(filename = file.path(pair_dir, "All_datasets_scatter.pdf"),
           plot = p, width = 2, height = 2)
  }
}

# ==================== Step 7: 汇总所有相关性结果 ====================
cat("\n========================================\n")
cat("Step 7: 汇总结果...\n")

# 保存总体结果
if (length(all_correlation_results) > 0) {
  summary_df <- do.call(rbind, all_correlation_results)
  summary_df <- summary_df[order(-abs(summary_df$Spearman_rho)), ]

  write.csv(summary_df, file.path(output_dir, "all_TF_TE_correlation_summary.csv"), row.names = FALSE)

  cat("\n==================== 总体结果 ====================\n")
  print(summary_df)

  cat("\n显著相关 (p<0.05):", sum(summary_df$P_value < 0.05), "/", nrow(summary_df), "\n")
  cat("正相关:", sum(summary_df$Spearman_rho > 0 & summary_df$P_value < 0.05), "\n")
  cat("负相关:", sum(summary_df$Spearman_rho < 0 & summary_df$P_value < 0.05), "\n")
}

# 保存每个Dataset的结果
if (length(per_dataset_results) > 0) {
  per_dataset_df <- do.call(rbind, per_dataset_results)
  per_dataset_df <- per_dataset_df[order(per_dataset_df$TF, per_dataset_df$TE, -abs(per_dataset_df$Spearman_rho)), ]

  write.csv(per_dataset_df, file.path(output_dir, "per_dataset_TF_TE_correlation.csv"), row.names = FALSE)

  cat("\n==================== 每个Dataset结果 ====================\n")
  cat("共", nrow(per_dataset_df), "个TF-TE-Dataset组合\n")
  cat("显著相关 (p<0.05):", sum(per_dataset_df$P_value < 0.05), "/", nrow(per_dataset_df), "\n")
  cat("正相关:", sum(per_dataset_df$Spearman_rho > 0 & per_dataset_df$P_value < 0.05), "\n")
  cat("负相关:", sum(per_dataset_df$Spearman_rho < 0 & per_dataset_df$P_value < 0.05), "\n")
}

cat("\n结果保存在:", output_dir, "\n")
cat("========================================\n")
