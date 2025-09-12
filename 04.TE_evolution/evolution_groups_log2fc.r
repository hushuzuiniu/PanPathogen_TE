library(data.table)
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)

setwd("/data2t_2/pathogen_TE_2025_New/09.TE_evolution/")
#TE Loci DE results
auto_file_paths <- list.files(path = "/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_loci/DE_results/", 
                              pattern = "*all_res_TEs\\.csv$", 
                              recursive = TRUE, 
                              full.names = TRUE)

if (length(auto_file_paths) > 0) {
  TE_loci_all_data <- do.call(rbind, lapply(auto_file_paths, function(file) {
    df <- read.csv(file, stringsAsFactors = FALSE)
    df$file_source <- file
    df$file_name <- basename(file)
    return(df)
  }))
}

# filer pvalue>=0.05
TE_loci_all_data_filter<-TE_loci_all_data[TE_loci_all_data$pvalue<=0.05,]

#TE subfamily DE results
TEs_auto_file_paths <- list.files(path = "/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_subfamily_feature/DE_results/", 
                                  pattern = "*all_res_TEs_feature\\.csv$", 
                                  recursive = TRUE, 
                                  full.names = TRUE)

if (length(TEs_auto_file_paths) > 0) {
  TEs_all_data <- do.call(rbind, lapply(TEs_auto_file_paths, function(file) {
    df <- read.csv(file, stringsAsFactors = FALSE)
    df$file_source <- file
    df$file_name <- basename(file)
    return(df)
  }))
}
# filer pvalue>=0.05
TEs_all_data_filter<-TEs_all_data[TEs_all_data$pvalue<=0.05,]

##TE loci level
TE_loci<-read.csv("Fig2_All_species_Up_DE-TEs_Loci_padj0.05log2FC1_v2_intergenic_top9.csv")
TE_loci<-TE_loci$loci_name
TE_loci_all_data_filter<-TE_loci_all_data_filter[TE_loci_all_data_filter$X %in% TE_loci, ]

##TE subfamily features level
TE_subfamily<-read.csv("Fig2_All_species_Up_DE-TEs_feature_padj0.05log2FC1_v2_intergenic_top5.csv")
TE_subfamily<-TE_subfamily$Subfamily
TEs_all_data_filter<-TEs_all_data_filter[TEs_all_data_filter$X %in% TE_subfamily,]


TE_loci_clean <- TE_loci_all_data_filter %>%
  mutate(TE_clean = str_replace(X, "_c\\d+$", ""))  #

TEs_clean <- TEs_all_data_filter %>%
  mutate(
    TE_clean = str_replace_all(X, "_(intergenic|intron)$", ""))  

TE_evolution_group<-read.csv("/data2t_2/pathogen_TE_2025_New/09.TE_evolution/TE_n55_age_group_v1.csv")
# TE_evolution_group<-read.csv("/data2t_2/pathogen_TE_2025_New/09.TE_evolution/TE_n55_age_group_v2.csv")
# TE_evolution_group<-read.csv("/data2t_2/pathogen_TE_2025_New/09.TE_evolution/TE_n55_age_group_v3.csv")
all_TE_data <- bind_rows(
  TE_loci_clean %>% select(TE_clean, log2FoldChange, file_name) %>% mutate(data_type = "TE_loci"),
  TEs_clean %>% select(TE_clean, log2FoldChange, file_name) %>% mutate(data_type = "TEs")
)

TE_with_group <- all_TE_data %>%
  left_join(TE_evolution_group, by = c("TE_clean" = "TE")) 

TE_with_group_TE_loci<-TE_with_group[TE_with_group$data_type=="TE_loci",]
TE_with_group_TEs<-TE_with_group[TE_with_group$data_type=="TEs",]
####################################TE subfamily####################################
p1 <- ggplot(TE_with_group_TEs, aes(x = log2FoldChange, y = factor(Group, levels = rev(levels(factor(Group)))), fill = factor(Group))) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.3,outlier.alpha = 0.6,size = 0.2,width = 0.7) +
  geom_jitter(alpha = 0.3, size = 0.3, height = 0.1,color = "gray40") +
  geom_vline(xintercept = 0,linetype = "dashed", color = "black", size = 0.3) +
  # coord_cartesian(xlim = c(-5, 10)) +
  labs(
    title = "TE Subfamily Evolution Groups Expression Distribution",
    y = "Evolution Group",
    x = "Log2 Fold Change") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 8, face = "bold", margin = margin(b = 15)),
    axis.title.x = element_text(size = 7, face = "bold", margin = margin(t = 8)),
    axis.title.y = element_text(size = 7, face = "bold", margin = margin(r = 8)),
    axis.text.x = element_text(size = 6, color = "black"),
    axis.text.y = element_text(size = 6, color = "black", face = "bold"),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.3),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(t = 10, r =10, b = 10, l = 10)  
  ) +
  scale_fill_manual(values = c("#C4E1E6", "#A4CCD9", "#8DBCC7", "#648DB3", "#7F8CAA")) +
  scale_x_continuous(breaks = seq(-4, 10, 2), expand = c(0.02, 0))
p1
kruskal_result <- kruskal.test(log2FoldChange ~ Group, data = TE_with_group_TEs)
pairwise_result <- pairwise.wilcox.test(TE_with_group_TEs$log2FoldChange, 
                                        TE_with_group_TEs$Group, 
                                        p.adjust.method = "BH")
print(pairwise_result)
get_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

create_bracket_annotations <- function(pairwise_result, y_levels) {
  p_matrix <- pairwise_result$p.value
  significant_comparisons <- data.frame()
  
  for(i in 1:nrow(p_matrix)) {
    for(j in 1:ncol(p_matrix)) {
      p_val <- p_matrix[i, j]
      stars <- get_stars(p_val)
      if(!is.na(p_val) && stars != "ns" && stars != "") {
        group1 <- rownames(p_matrix)[i]
        group2 <- colnames(p_matrix)[j]
        
        y1 <- which(y_levels == group1)
        y2 <- which(y_levels == group2)
        
        significant_comparisons <- rbind(significant_comparisons, data.frame(
          y1 = y1,
          y2 = y2,
          group1 = group1,
          group2 = group2,
          p_value = p_val,
          stars = stars,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  significant_comparisons <- significant_comparisons[order(significant_comparisons$p_value), ]
  
  return(significant_comparisons)
}

x_range <- range(TE_with_group_TEs$log2FoldChange, na.rm = TRUE)
x_max <- x_range[2]
x_min <- x_range[1]

y_levels <- rev(levels(factor(TE_with_group_TEs$Group)))
bracket_data <- create_bracket_annotations(pairwise_result, y_levels)

print(bracket_data[, c("group1", "group2", "p_value", "stars")])

p1_with_brackets <- p1
if(nrow(bracket_data) > 0) {
  bracket_spacing <- 0.5 
  start_x <- x_max + 1     
  
  for(i in 1:nrow(bracket_data)) {
    y_pos <- c(bracket_data$y1[i], bracket_data$y2[i])
    x_pos <- start_x + (i-1) * bracket_spacing
    
    y_min <- min(y_pos)
    y_max <- max(y_pos)
    
    p1_with_brackets <- p1_with_brackets +
      annotate("segment", 
               x = x_pos, xend = x_pos,
               y = y_min, yend = y_max,
               size = 0.2, color = "#4d4d4d") +
      annotate("segment", 
               x = x_pos - 0.15, xend = x_pos,
               y = y_min, yend = y_min,
               size = 0.2, color = "#4d4d4d") +
      annotate("segment", 
               x = x_pos - 0.15, xend = x_pos,
               y = y_max, yend = y_max,
               size = 0.2, color = "#4d4d4d") +
      annotate("text", 
               x = x_pos + 0.2, 
               y = (y_min + y_max) / 2,
               label = bracket_data$stars[i],
               size = 2, color = "#4d4d4d", fontface = "bold")
  }
  
  # p1_with_brackets <- p1_with_brackets +
  #   annotate("text", 
  #            x = start_x, 
  #            y = length(y_levels) + 0.3,
  #            label = paste("Kruskal-Wallis\np =", format(kruskal_result$p.value, digits = 3)),
  #            size = 2.5, 
  #            fontface = "bold",
  #            hjust = 0)
  
  max_bracket_x <- start_x + (nrow(bracket_data)-1) * bracket_spacing + 1
  p1_with_brackets <- p1_with_brackets +
    coord_cartesian(xlim = c(x_min, max_bracket_x))
}
print(p1_with_brackets)

ggsave("./03.expression_foldchange/TE_subfamily_evolution_groups_boxplot_v1.pdf",p1_with_brackets,width = 3.2,height = 2.5)
# ggsave("./03.expression_foldchange/TE_subfamily_evolution_groups_boxplot_v2.pdf",p1,width = 3,height = 2.5)
# ggsave("./03.expression_foldchange/TE_subfamily_evolution_groups_boxplot_v3.pdf",p1,width = 3,height = 2.5)
####################################TE Loci####################################
p2 <- ggplot(TE_with_group_TE_loci, aes(x = log2FoldChange, y = factor(Group, levels = rev(levels(factor(Group)))), fill = factor(Group))) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.3,outlier.alpha = 0.6,size = 0.2,width = 0.7) +
  geom_jitter(alpha = 0.3, size = 0.3, height = 0.1,color = "gray40") +
  geom_vline(xintercept = 0,linetype = "dashed", color = "black", size = 0.3) +
  coord_cartesian(xlim = c(-5, 12)) +
  labs(
    title = "TE Loci Evolution Groups Expression Distribution",
    y = "Evolution Group",
    x = "Log2 Fold Change") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 8, face = "bold", margin = margin(b = 15)),
    axis.title.x = element_text(size = 7, face = "bold", margin = margin(t = 8)),
    axis.title.y = element_text(size = 7, face = "bold", margin = margin(r = 8)),
    axis.text.x = element_text(size = 6, color = "black"),
    axis.text.y = element_text(size = 6, color = "black", face = "bold"),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.3),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(t = 10, r =10, b = 10, l = 10)  
  ) +
  scale_fill_manual(values = c("#C4E1E6", "#A4CCD9", "#8DBCC7", "#648DB3", "#7F8CAA")) +
  scale_x_continuous(breaks = seq(-4, 10, 2), expand = c(0.02, 0))
p2
kruskal_result <- kruskal.test(log2FoldChange ~ Group, data = TE_with_group_TE_loci)
cat("Kruskal-Wallis检验:\n")
cat("H =", round(kruskal_result$statistic, 3), ", p-value =", format(kruskal_result$p.value, scientific = TRUE), "\n\n")
cat("成对比较 (Wilcoxon rank-sum test with Benjamini-Hochberg correction):\n")
pairwise_result <- pairwise.wilcox.test(TE_with_group_TE_loci$log2FoldChange, 
                                        TE_with_group_TE_loci$Group, 
                                        p.adjust.method = "BH")
print(pairwise_result)
x_range <- range(TE_with_group_TE_loci$log2FoldChange, na.rm = TRUE)
x_max <- x_range[2]
x_min <- x_range[1]

y_levels <- rev(levels(factor(TE_with_group_TE_loci$Group)))
bracket_data <- create_bracket_annotations(pairwise_result, y_levels)

print(bracket_data[, c("group1", "group2", "p_value", "stars")])

p2_with_brackets <- p2
if(nrow(bracket_data) > 0) {
  bracket_spacing <- 0.5 
  start_x <- x_max + 1     
  start_x<-12
  for(i in 1:nrow(bracket_data)) {
    y_pos <- c(bracket_data$y1[i], bracket_data$y2[i])
    x_pos <- start_x + (i-1) * bracket_spacing
    y_min <- min(y_pos)
    y_max <- max(y_pos)
    
    p2_with_brackets <- p2_with_brackets +
      annotate("segment", 
               x = x_pos, xend = x_pos,
               y = y_min, yend = y_max,
               size = 0.2, color = "#4d4d4d") +
      annotate("segment", 
               x = x_pos - 0.15, xend = x_pos,
               y = y_min, yend = y_min,
               size = 0.2, color = "#4d4d4d") +
      annotate("segment", 
               x = x_pos - 0.15, xend = x_pos,
               y = y_max, yend = y_max,
               size = 0.2, color = "#4d4d4d") +
      annotate("text", 
               x = x_pos + 0.2, 
               y = (y_min + y_max) / 2,
               label = bracket_data$stars[i],
               size = 2, color = "#4d4d4d", fontface = "bold")
  }
  
  # p2_with_brackets <- p2_with_brackets +
  #   annotate("text", 
  #            x = start_x, 
  #            y = length(y_levels) + 0.3,
  #            label = paste("Kruskal-Wallis\np =", format(kruskal_result$p.value, digits = 3)),
  #            size = 2.5, 
  #            fontface = "bold",
  #            hjust = 0)
  
  max_bracket_x <- start_x + (nrow(bracket_data)-1) * bracket_spacing + 1
  p2_with_brackets <- p2_with_brackets +
    coord_cartesian(xlim = c(x_min, max_bracket_x))
}
print(p2_with_brackets)


ggsave("./03.expression_foldchange/TE_Loci_evolution_groups_boxplot_v1.pdf",p2_with_brackets,width = 3.2,height = 2.5)
# ggsave("./03.expression_foldchange/TE_Loci_evolution_groups_boxplot_v2.pdf",p2,width = 3,height = 2.5)
# ggsave("./03.expression_foldchange/TE_Loci_evolution_groups_boxplot_v3.pdf",p2,width = 3,height = 2.5)

################################################################################################################
Bacteria_Lmo_1_Post_120m_vs_Pre_120m<-TE_with_group[TE_with_group$file_name=="Bacteria_Lmo_1_Post_120m_vs_Pre_120m_all_res_TEs.csv",]
p3 <- ggplot(Bacteria_Lmo_1_Post_120m_vs_Pre_120m, aes(x = factor(Group), y = log2FoldChange, fill = factor(Group))) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  # geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  # coord_cartesian(ylim = c(-10, 10)) +  
  labs(
    title = "Bacteria_Lmo_1_Post_120m_vs_Pre_120m",
    x = "Evolution Group",
    y = "Log2 Fold Change",
    fill = "Group"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0),
    legend.position = "none"
  ) +
  scale_fill_brewer(type = "qual", palette = "Set3")

print(p3)



