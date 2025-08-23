library(data.table)
library(readxl)
library(openxlsx)
library(ggplot2)
setwd("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_loci/DE_results")
df_A549 <- read.csv("A549/A549_split/01.A549_split_All_sig_DE-TEs_padj0.05log2FC1.csv", row.names = 1, stringsAsFactors = FALSE)
df_B_cell <- read.csv("B_cell/B_cell_split/01.B_cell_split_All_sig_DE-TEs_padj0.05log2FC1.csv", row.names = 1, stringsAsFactors = FALSE)
df_Dendritic <- read.csv("Dendritic/Dendritic_split/01.Dendritic_split_All_sig_DE-TEs_padj0.05log2FC1.csv", row.names = 1, stringsAsFactors = FALSE)
df_Liver <- read.csv("Liver/Liver_split/01.Liver_split_All_sig_DE-TEs_padj0.05log2FC1.csv", row.names = 1, stringsAsFactors = FALSE)
df_Macrophages <- read.csv("Macrophages/Macrophages_split/01.Macrophages_split_All_sig_DE-TEs_padj0.05log2FC1.csv", row.names = 1, stringsAsFactors = FALSE)
df_Monocyte <- read.csv("Monocyte/Monocyte_split/01.Monocyte_split_All_sig_DE-TEs_padj0.05log2FC1.csv", row.names = 1, stringsAsFactors = FALSE)
df_PBMCs <- read.csv("PBMCs/PBMCs_split/01.PBMCs_split_All_sig_DE-TEs_padj0.05log2FC1.csv", row.names = 1, stringsAsFactors = FALSE)
df_T_cell <- read.csv("T_cell/T_cell_split/01.T_cell_split_All_sig_DE-TEs_padj0.05log2FC1.csv", row.names = 1, stringsAsFactors = FALSE)
df_hSAEC <- read.csv("hSAEC/hSAEC_split/01.hSAEC_split_All_sig_DE-TEs_padj0.05log2FC1.csv", row.names = 1, stringsAsFactors = FALSE)
df_Huh_7 <- read.csv("Huh-7/Huh-7_split/01.Huh-7_split_All_sig_DE-TEs_padj0.05log2FC1.csv", row.names = 1, stringsAsFactors = FALSE)


colnames(df_A549) <- paste0("A549_", colnames(df_A549))
colnames(df_B_cell) <- paste0("B_cell_", colnames(df_B_cell))
colnames(df_Dendritic) <- paste0("Dendritic_", colnames(df_Dendritic))
colnames(df_Liver) <- paste0("Liver_", colnames(df_Liver))
colnames(df_Macrophages) <- paste0("Macrophages_", colnames(df_Macrophages))
colnames(df_Monocyte) <- paste0("Monocyte_", colnames(df_Monocyte))
colnames(df_PBMCs) <- paste0("PBMCs_", colnames(df_PBMCs))
colnames(df_T_cell) <- paste0("T_cell_", colnames(df_T_cell))
colnames(df_hSAEC) <- paste0("hSAEC_", colnames(df_hSAEC))
colnames(df_Huh_7) <- paste0("Huh-7_", colnames(df_Huh_7))


dfs <- list(df_A549, df_B_cell, df_Dendritic, df_Liver, 
            df_Macrophages, df_Monocyte, df_PBMCs, df_T_cell,df_hSAEC,df_Huh_7)

dfs <- lapply(dfs, function(df) {
  df$GeneID <- rownames(df)
  rownames(df) <- NULL
  return(df)
})

merged_df <- Reduce(function(x, y) merge(x, y, by = "GeneID", all = TRUE), dfs)
write.csv(merged_df, "/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_loci/plots/All_species_split_DE-TEs_Loci_padj0.05log2FC1.csv", row.names = F)
merged_df <- read.csv("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_loci/plots/All_species_split_DE-TEs_Loci_padj0.05log2FC1.csv")
merge_columns <- function(df) {
  result_df <- data.frame(GeneID = df$GeneID, stringsAsFactors = FALSE)
  data_cols <- colnames(df)[-1]
  
  extract_base_pattern <- function(col) {
    if (grepl("SARSCoV2", col)) {
      parts <- strsplit(col, "_")[[1]]
      return(paste(parts[1], parts[2], "SARSCoV2", sep="_"))
    }
    if (grepl("HIV1", col)) {
      parts <- strsplit(col, "_")[[1]]
      return(paste(parts[1], parts[2], "HIV1", sep="_"))
    }
    pattern <- gsub("_\\d+$", "", col)
    return(pattern)
  }
  
  patterns <- sapply(data_cols, extract_base_pattern)
  unique_patterns <- unique(patterns)
  library(data.table)
  dt <- as.data.table(df)
    for (pattern in unique_patterns) {
    matching_cols <- data_cols[patterns == pattern]
    result_df[[pattern]] <- apply(dt[, matching_cols, with=FALSE], 1, function(x) {
      if(all(x == "Up")) {
        return("Up")
      } else if(all(x == "Down")) {
        return("Down")
      } else {
        return("/")
      }
    })
  }  
  return(result_df)
}
result_df <- merge_columns(merged_df)
write.csv(result_df, "/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_loci/plots/All_species_DE-TEs_Loci_padj0.05log2FC1.csv", row.names = F)

##############################################################################
setwd("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_loci/plots/")
P1 <- read.csv("All_species_DE-TEs_Loci_padj0.05log2FC1.csv", header = TRUE)

TE_Loci_info <- fread("/data2t_2/pathogen_TE_2025_New/01.Genomic_features_Gencode/hg38_TE_anno_custom_v20240110_0based_with_feature_summary_v2.bed")
TE_Loci_info <- TE_Loci_info[, c("gene_id", "transcript_id", "class_id", "family_id", "prioritized_feature")]
TE_Loci_info$region <- TE_Loci_info$prioritized_feature
TE_Loci_info$region[TE_Loci_info$prioritized_feature %in%
                      c("3UTR", "5UTR", "CDS", "noncoding_exon")] <- "exon"

P1_long <- reshape2::melt(P1, id.vars = "GeneID", 
                          variable.name = "Species", 
                          value.name = "Direction")

P1_filtered <- P1_long %>%
  filter(Direction %in% c("Up", "Down"))

merged_data <- merge(P1_filtered, TE_Loci_info, by.x = "GeneID", by.y = "transcript_id")

merged_data$region <- factor(merged_data$region, levels = c("exon", "intron", "intergenic"))

region_count_data <- merged_data %>%
  group_by(Species, Direction, region) %>%
  summarise(Count = n(), .groups = "drop")

region_count_data$Count <- ifelse(region_count_data$Direction == "Down", 
                                  -region_count_data$Count, 
                                  region_count_data$Count)

region_count_data$fill_group <- paste(region_count_data$Direction, region_count_data$region, sep="_")

region_count_data$fill_group <- factor(region_count_data$fill_group, 
                                       levels = c("Up_exon", "Up_intron", "Up_intergenic",
                                                  "Down_exon", "Down_intron", "Down_intergenic"))


fill_colors <- c(
  "Up_exon" = "#B82132",      
  "Up_intron" = "#D2665A",    
  "Up_intergenic" = "#F2B28C", 
  "Down_exon" = "#2D336B",    
  "Down_intron" = "#7886C7",  
  "Down_intergenic" = "#A9B5DF" 
)

p1 <- ggplot(region_count_data, aes(x = Species, y = Count, fill = fill_group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = fill_colors,
                    name = "",
                    labels = c("Up_exon", "Up_intron", "Up_intergenic",
                               "Down_exon", "Down_intron", "Down_intergenic")) +
  scale_y_continuous(
    breaks = seq(-210000, 210000, 30000),  
    labels = function(x) abs(x)           
  ) +
  theme_minimal() +
  labs(
    title = "Number of differentially expressed TEs Loci by genomic feature(padj0.05&log2FC1)严格",
    x = "Species",
    y = "No. of TE Loci"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "right",
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.line = element_line(color = "black", size = 0.5)
  ) +
  geom_hline(yintercept = 0, color = "black", size = 0.5)

p1
label_data <- region_count_data %>% 
  group_by(Species, Direction) %>% 
  summarise(
    Total = sum(abs(Count)),
    y_position = ifelse(Direction[1] == "Up", 
                        sum(Count) + 5000, 
                        sum(Count) - 5000), 
    .groups = "drop"
  )


p2<-p1 + geom_text(
  data = label_data,
  mapping = aes(x = Species, y = y_position, label = Total),
  size = 2, 
  # fontface = "bold",
  inherit.aes = FALSE
)
p2

extract_cell_type <- function(col_name) {
  if (startsWith(col_name, "B_cell")) {
    return("B_cell")
  } else if (startsWith(col_name, "T_cell")) {
    return("T_cell")
  } else {
    return(sub("_.*", "", col_name))
  }
}

celltype_data <- data.frame(
  Species = unique(region_count_data$Species),
  celltype = sapply(as.character(unique(region_count_data$Species)), extract_cell_type)
)
celltype_positions <- celltype_data %>%
  group_by(celltype) %>%
  summarise(
    start = min(which(unique(region_count_data$Species) == first(Species))),
    end = max(which(unique(region_count_data$Species) == last(Species))),
    center = (start + end) / 2
  )

p3 <- p2 +
  geom_vline(
    data = celltype_positions,
    aes(xintercept = end + 0.5),
    linetype = "dashed",
    color = "gray60",
    size = 0.5
  ) +
  annotate(
    "text",
    x = celltype_positions$center,
    y = max(region_count_data$Count) * 1.55, 
    label = celltype_positions$celltype,
    size = 2.5,
    fontface = "bold"
  )

# ggsave("TE_Loci_expression_celltype&species_by_feature_padj0.05log2FC1.pdf", p3, width = 10, height = 6)