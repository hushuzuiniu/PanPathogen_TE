# HOMER Per-Subfamily TF Enrichment 点图热图可视化
# 使用P-value < 1e-5 和 Fold >= 3 筛选
# TF类型从注释文件获取
# Date: 2024-12-09

library(tidyverse)
library(data.table)
library(RColorBrewer)
library(patchwork)
library(readxl)

setwd("/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_prediction_add_bg_v3/")

# =============================================================================
# 1. 读取数据
# =============================================================================

cat("=== Step 1: 读取数据 ===\n")

# 读取合并的HOMER结果
combined_df <- fread("combined_HOMER_results_all_subfamilies.txt", header = TRUE)

cat("原始数据维度:", nrow(combined_df), "行\n")
cat("Subfamilies:", length(unique(combined_df$Subfamily)), "\n")
cat("TF motifs:", length(unique(combined_df$TF_Name)), "\n")

# =============================================================================
# 2. 筛选显著结果 (P < 1e-5, Fold >= 3)
# =============================================================================

cat("\n=== Step 2: 筛选数据 (P < 1e-5, Fold >= 3) ===\n")

sig_df <- combined_df %>%
  filter(P_value < 1e-5) %>%
  filter(Fold_Enrichment >= 3) %>%
  filter(Target_with_Motif_pct >= 10)

cat("筛选后组合数:", nrow(sig_df), "\n")

# =============================================================================
# 3. 简化TF_Name并匹配TF类型
# =============================================================================

cat("\n=== Step 3: 简化TF名称并匹配类型 ===\n")

sig_df$TF_Name <- as.character(sig_df$TF_Name)

# 简化TF_Name用于显示：
# OCT:OCT(POU,Homeobox,IR1)/NPC-Brn2-ChIP-Seq(GSE35496)/Homer -> OCT:OCT(POU,Homeobox,IR1)
# Pax7(Paired,Homeobox),longest/Myoblast-Pax7-ChIP-Seq(GSE25064)/Homer -> Pax7(Paired,Homeobox),longest
# 保留TF名称和domain信息，去掉实验细节
simplify_tf_name_display <- function(tf_name) {
  # 取第一个 / 之前的部分（TF名称+domain信息）
  simple <- sub("/.*", "", tf_name)
  return(simple)
}

# 简化TF_Name用于匹配TF类型数据库：
# OCT:OCT(POU,Homeobox,IR1) -> OCT:OCT
# Pax7(Paired,Homeobox),longest -> Pax7
simplify_tf_name_match <- function(tf_name) {
  # 先用display简化
  simple <- sub("/.*", "", tf_name)
  # 再去掉括号内容和逗号后缀用于匹配
  simple <- sub("\\(.*", "", simple)
  simple <- sub(",.*", "", simple)
  return(simple)
}

# 从Motif_Name括号中提取TF类型（作为fallback）
# Brn1(POU,Homeobox)/... -> POU
# FXR(NR),IR1/... -> NR
# GATA(Zf),IR3/... -> Zf
extract_tf_type_from_motif <- function(motif_name) {
  # 提取第一个括号内的内容
  match <- regmatches(motif_name, regexpr("\\([^)]+\\)", motif_name))
  if (length(match) == 0 || match == "") return(NA)
  # 去掉括号
  type_str <- gsub("[()]", "", match)
  # 取第一个逗号之前的部分（主要类型）
  type_str <- sub(",.*", "", type_str)
  return(type_str)
}

# HOMER类型到标准类型的映射
homer_type_mapping <- c(
  "bZIP" = "bZIP",
  "bHLH" = "bHLH",
  "NR" = "NR",
  "Zf" = "C2H2 ZF",
  "ETS" = "Ets",
  "IRF" = "IRF",
  "Forkhead" = "Forkhead",
  "Homeobox" = "Homeodomain",
  "POU" = "Homeodomain; POU",
  "GATA" = "GATA",
  "HMG" = "HMG/Sox",
  "Stat" = "STAT",
  "T-box" = "T-box",
  "RHD" = "Rel",
  "MAD" = "SMAD",
  "HSF" = "HSF",
  "EBF" = "EBF",
  "Paired" = "Paired box",
  "?" = NA,
  "AP2" = "AP-2",
  "p53" = "p53",
  "E2F" = "E2F",
  "Myb" = "Myb/SANT",
  "RFX" = "RFX",
  "TEA" = "TEA",
  "Runt" = "Runt",
  "CP2" = "CP2"
)

# 用于显示的简化名称（保留domain信息，区分变体）
# 注意：从 Motif_Name 提取，而不是 TF_Name，因为变体信息只在 Motif_Name 中
sig_df$TF_Name_display <- sapply(sig_df$Motif_Name, simplify_tf_name_display)
# 用于匹配TF类型的简化名称
sig_df$TF_Name_match <- sapply(sig_df$Motif_Name, simplify_tf_name_match)

cat("原始唯一TF数:", length(unique(sig_df$TF_Name)), "\n")
cat("显示用唯一TF数:", length(unique(sig_df$TF_Name_display)), "\n")
cat("匹配用唯一TF数:", length(unique(sig_df$TF_Name_match)), "\n")

# 读取注释文件
cell_anno <- read_excel("/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_annotation/DatabaseExtract_v_1.01.xlsx")
cell_anno <- cell_anno[, c("HGNC symbol", "DBD")]

KRAB_list <- read_csv("/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_annotation/KRFP_list_combined_uniq_n397.csv", show_col_types = FALSE)

# 添加TF_type列
sig_df$TF_type <- NA

# 步骤1: 用匹配名称匹配 KRAB_list（忽略大小写）
krab_tfs_upper <- toupper(KRAB_list$TF)

for(i in 1:nrow(sig_df)) {
  if(toupper(sig_df$TF_Name_match[i]) %in% krab_tfs_upper) {
    sig_df$TF_type[i] <- "KRAB"
  }
}

# 步骤2: 对于没有匹配到 KRAB 的，用匹配名称匹配 cell_anno 的 DBD（忽略大小写）
cell_anno_lookup <- cell_anno
cell_anno_lookup$`HGNC symbol` <- toupper(cell_anno_lookup$`HGNC symbol`)

for(i in 1:nrow(sig_df)) {
  if(is.na(sig_df$TF_type[i])) {
    motif_upper <- toupper(sig_df$TF_Name_match[i])
    match_idx <- which(cell_anno_lookup$`HGNC symbol` == motif_upper)
    if(length(match_idx) > 0) {
      sig_df$TF_type[i] <- cell_anno$DBD[match_idx[1]]
    }
  }
}

# 简化一些类型名称
sig_df$TF_type[sig_df$TF_type == "Nuclear receptor"] <- "NR"

# Fallback: 从Motif_Name括号中提取TF类型
cat("数据库匹配后NA数量:", sum(is.na(sig_df$TF_type)), "\n")

for(i in 1:nrow(sig_df)) {
  if(is.na(sig_df$TF_type[i])) {
    # 从括号中提取类型
    homer_type <- extract_tf_type_from_motif(sig_df$Motif_Name[i])
    if(!is.na(homer_type) && homer_type %in% names(homer_type_mapping)) {
      sig_df$TF_type[i] <- homer_type_mapping[homer_type]
    }
  }
}

cat("Fallback后NA数量:", sum(is.na(sig_df$TF_type)), "\n")

# 查看匹配情况
cat("TF类型匹配情况:\n")
print(table(sig_df$TF_type, useNA = "ifany"))

# =============================================================================
# 4. 选择top TFs (按Fold Enrichment排序)
# =============================================================================

cat("\n=== Step 4: 选择显示的TFs ===\n")

# 统计每个TF的最大Fold Enrichment
tf_counts <- sig_df %>%
  group_by(TF_Name) %>%
  summarise(
    n_subfamilies = n_distinct(Subfamily),
    max_fold = max(Fold_Enrichment, na.rm = TRUE),
    mean_fold = mean(Fold_Enrichment, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(max_fold), desc(mean_fold))

# 选择top 100 TF (按最大Fold Enrichment排序)
top_tfs <- tf_counts %>%
  head(100) %>%
  pull(TF_Name)

cat("选择的TFs数量:", length(top_tfs), "\n")

# 筛选数据
heatmap_data <- sig_df %>%
  filter(TF_Name %in% top_tfs)

cat("筛选后 TE subfamilies:", length(unique(heatmap_data$Subfamily)), "\n")
cat("筛选后 TFs (原始):", length(unique(heatmap_data$TF_Name)), "\n")
cat("筛选后 TFs (显示):", length(unique(heatmap_data$TF_Name_display)), "\n")
cat("去重前 TE-TF组合数:", nrow(heatmap_data), "\n")

# 去重：同一个 Subfamily + TF_Name_display，保留最显著的结果（P值最小）
# 这处理的是同一TF来自不同实验的情况（如 HRE(HSF) 来自不同ChIP-seq实验）
heatmap_data <- heatmap_data %>%
  group_by(Subfamily, TF_Name_display) %>%
  slice_min(P_value, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("去重后 TE-TF组合数:", nrow(heatmap_data), "\n")

# =============================================================================
# 5. 排序
# =============================================================================

cat("\n=== Step 5: 排序 ===\n")

# 使用 TF_Name_display 用于显示（保留domain信息，区分变体）
# 按TF的最大Fold Enrichment排序
tf_order <- heatmap_data %>%
  group_by(TF_Name_display) %>%
  summarise(max_fold = max(Fold_Enrichment, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(max_fold)) %>%
  pull(TF_Name_display)

# 按subfamily有多少个显著TF排序，相同时按平均Fold Enrichment降序
te_order <- heatmap_data %>%
  group_by(Subfamily) %>%
  summarise(n = n(),
            mean_fold = mean(Fold_Enrichment, na.rm = TRUE),
            .groups = 'drop') %>%
  arrange(desc(n), desc(mean_fold)) %>%
  pull(Subfamily)

heatmap_data$Subfamily <- factor(heatmap_data$Subfamily, levels = rev(te_order))
heatmap_data$TF_Name_display <- factor(heatmap_data$TF_Name_display, levels = tf_order)

# 计算-log10(P)用于颜色
heatmap_data$neg_log10_p <- -log10(heatmap_data$P_value + 1e-300)
# Cap at 50 for visualization
heatmap_data$neg_log10_p_capped <- pmin(heatmap_data$neg_log10_p, 50)

# =============================================================================
# 6. TF类型颜色
# =============================================================================

tf_type_colors <- c(
  "KRAB" = "#9E9DCB",
  "C2H2 ZF" = "#C1747B",
  "NR" = "#E1C270",
  "AP-2" = "#929F74",
  "bHLH" = "#DEA368",
  "Ets" = "#57B3C3",
  "EBF" = "#CF8DB4",
  "IRF" = "#636363",
  "Rel" = "#FD79A8",
  "SMAD" = "#A0E7E5",
  "STAT" = "#2A4492",
  "Forkhead" = "#87CEEB",
  "HMG/Sox" = "#F0E68C",
  "GATA" = "#FFA07A",
  "Homeodomain" = "#DDA0DD",
  "Homeodomain; POU" = "#C9A0DC",
  "T-box" = "#98FB98",
  "p53" = "#FF6B6B",
  "E2F" = "#4ECDC4",
  "CP2" = "#95E1D3",
  "Viral" = "#FF4757",
  "bZIP" = "#F8B500",
  "Paired box" = "#A8E6CF",
  "RFX" = "#DDA15E",
  "Myb/SANT" = "#BC6C25"
)

# =============================================================================
# 7. 绘制点图热图 (基本版) - 点大小=Fold Enrichment, 颜色=-log10(P)
# =============================================================================

cat("\n=== Step 6: 绘制点图热图 ===\n")

p_heatmap <- ggplot(heatmap_data, aes(x = TF_Name_display, y = Subfamily)) +
  geom_point(aes(size = Fold_Enrichment,
                 fill = neg_log10_p_capped),
             shape = 21,
             stroke = 0.1) +
  scale_size_continuous(name = "Fold Enrichment",
                        range = c(1, 4),
                        breaks = c(3, 5, 10, 50, 100),
                        limits = c(2, NA)) +
  scale_fill_gradient(name = "-log10(P)",
                      low = "#fcf3ef",
                      high = "#6A0B10",
                      limits = c(1, NA)) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = unit(0.1, "cm"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 9, color = "black", face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    panel.grid.major = element_line(color = "#ebe7e6", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "HOMER TF Motif Enrichment per TE Subfamily (vs Background)",
    x = "Transcription Factor",
    y = "TE Subfamily"
  )

ggsave("HOMER_per_subfamily_dotplot_heatmap_basic.pdf", p_heatmap, width = 14, height = 5)
cat("已保存: HOMER_per_subfamily_dotplot_heatmap_basic.pdf\n")

# =============================================================================
# 8. 绘制带TF type注释的热图
# =============================================================================

cat("\n=== Step 7: 绘制带注释的热图 ===\n")

# 创建TF注释数据
tf_annotation <- heatmap_data %>%
  dplyr::select(TF_Name_display, TF_type) %>%
  dplyr::distinct() %>%
  dplyr::arrange(match(TF_Name_display, tf_order))

tf_annotation$TF_Name_display <- factor(tf_annotation$TF_Name_display, levels = tf_order)

# 主热图
p_main <- ggplot(heatmap_data, aes(x = TF_Name_display, y = Subfamily)) +
  geom_point(aes(size = Fold_Enrichment,
                 fill = neg_log10_p_capped),
             shape = 21,
             stroke = 0.1) +
  scale_size_continuous(name = "Fold Enrichment",
                        range = c(1, 4),
                        breaks = c(3, 5, 10, 50, 100),
                        limits = c(2, NA)) +
  scale_fill_gradient(name = "-log10(P)",
                      low = "#fcf3ef",
                      high = "#6A0B10",
                      limits = c(1, NA)) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = unit(0.1, "cm"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 8, color = "black", face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    panel.grid.major = element_line(color = "#ebe7e6", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 2, r = 5, b = 5, l = 5)
  ) +
  labs(
    x = "Transcription Factor",
    y = "TE Subfamily"
  )

# 顶部注释条
p_anno <- ggplot(tf_annotation, aes(x = TF_Name_display, y = 1)) +
  geom_tile(aes(fill = TF_type), height = 1) +
  scale_fill_manual(values = tf_type_colors,
                    name = "TF Type",
                    na.value = "grey80",
                    guide = guide_legend(nrow = 1)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.x = unit(0.2, "cm"),
    plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
  )

# 组合
p_combined <- p_anno / p_main +
  plot_layout(heights = c(0.05, 1))

ggsave("HOMER_per_subfamily_dotplot_heatmap_with_anno.pdf", p_combined, width = 14, height = 5)
cat("已保存: HOMER_per_subfamily_dotplot_heatmap_with_anno.pdf\n")

# =============================================================================
# 9. 按TF type分面
# =============================================================================

cat("\n=== Step 8: 绘制分面热图 ===\n")

p_facet <- ggplot(heatmap_data, aes(x = TF_Name_display, y = Subfamily)) +
  geom_point(aes(size = Fold_Enrichment,
                 fill = neg_log10_p_capped),
             shape = 21,
             stroke = 0.1) +
  scale_size_continuous(name = "Fold Enrichment",
                        range = c(1, 4),
                        breaks = c(3, 5, 10, 50, 100),
                        limits = c(2, NA)) +
  scale_fill_gradient(name = "-log10(P)",
                      low = "#fcf3ef",
                      high = "#6A0B10",
                      limits = c(1, NA)) +
  facet_grid(. ~ TF_type, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = unit(0.1, "cm"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 8, color = "black", face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    panel.grid.major = element_line(color = "#ebe7e6", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 6, face = "bold",
                              margin = margin(t = 1.5, b = 1.5, unit = "pt")),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    strip.placement = "outside",
    panel.spacing = unit(0, "lines"),
    panel.spacing.x = unit(0, "lines"),
    plot.margin = margin(t = 2, r = 5, b = 5, l = 5)
  ) +
  labs(
    title = "HOMER TF Motif Enrichment by TF Type",
    x = "Transcription Factor",
    y = "TE Subfamily"
  )

ggsave("HOMER_per_subfamily_dotplot_heatmap_facet.pdf", p_facet, width = 14, height = 5)
cat("已保存: HOMER_per_subfamily_dotplot_heatmap_facet.pdf\n")

# =============================================================================
# 10. 保存数据
# =============================================================================

cat("\n=== Step 9: 保存数据 ===\n")

write_csv(heatmap_data, "HOMER_per_subfamily_heatmap_data.csv")
cat("已保存: HOMER_per_subfamily_heatmap_data.csv\n")

# 汇总
cat("\n")
cat("=", rep("=", 60), "\n", sep = "")
cat("VISUALIZATION SUMMARY\n")
cat("=", rep("=", 60), "\n", sep = "")
cat("\nTE subfamilies:", length(unique(heatmap_data$Subfamily)))
cat("\nTFs:", length(unique(heatmap_data$TF_Name)))
cat("\nTE-TF combinations:", nrow(heatmap_data))
cat("\n\nTF Type分布:\n")
print(table(heatmap_data$TF_type, useNA = "ifany"))
cat("\n完成!\n")
