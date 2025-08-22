library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(ggplot2)
library(viridis)
library(data.table)
library(readxl)
library(patchwork)

setwd("/data2t_2/pathogen_TE_2025_New/08.TF_prediction/")

bed_file <- "/data2t_2/hushu/00.ref/hg38_TE_anno_custom_v20240110_0base.bed"
bed_data <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE,
                       col.names = c("chr", "start", "end", "name", "score", "strand"))

bed_data$length <- bed_data$end - bed_data$start
bed_data$TE_subfamily <- gsub("_c\\d+$", "", bed_data$name)

result_base_dir <- "motif_analysis_results/"
# result_base_dir <- "motif_analysis_results_bak/"
result_dirs <- list.dirs(result_base_dir, recursive = FALSE)

# all_fimo_results <- data.frame()
# 
# for (dir in result_dirs) {
#   te_name <- basename(dir)
#   te_name <- gsub("_loci$", "", te_name)  
#   fimo_file <- file.path(dir, "fimo.tsv")
#   if (file.exists(fimo_file)) {
#     fimo_data <- read.delim(fimo_file, stringsAsFactors = FALSE)
# 
#     if (nrow(fimo_data) > 0) {
#       fimo_data$TE_subfamily <- te_name
#       all_fimo_results <- rbind(all_fimo_results, fimo_data)
#     }
#   } else {
#   }
# }
# 
# save(all_fimo_results, file = "all_fimo_results.rdata")
load("all_fimo_results.rdata")

all_fimo_results <- all_fimo_results %>%
  mutate(motif_alt_id = ifelse(
    motif_alt_id == "" | is.na(motif_alt_id),  
    str_extract(motif_id, "^[^.]+"),           
    motif_alt_id                               
  ))
te_subfamily<-read.csv("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_subfamily_feature/plots/Fig2_All_species_Up_DE-TEs_feature_padj0.05log2FC1_v2_intergenic_top5.csv",header = T)
gft_info<-fread("/data2t_2/pathogen_TE_2025_New/01.Genomic_features_Gencode/hg38_TE_anno_custom_v20240110_0based_with_feature_summary_v2.bed")
result <- data.table()

for(i in 1:nrow(te_subfamily)) {
  parts <- strsplit(te_subfamily$Subfamily[i], "_")[[1]]
  gene_name <- paste(parts[1:(length(parts)-1)], collapse="_")  
  feature_type <- parts[length(parts)]  
  
  matched_rows <- gft_info[gene_id == gene_name & prioritized_feature == feature_type]
  
  if(nrow(matched_rows) > 0) {
    result <- rbind(result, matched_rows)
  }
}

result<-result[,c("prioritized_feature","gene_id","transcript_id")]

filtered_fimo_loci <- all_fimo_results %>%
  mutate(loci_name = gsub("\\([+-]\\)$", "", sequence_name)) %>%
  filter(loci_name %in% result$transcript_id)

all_fimo_results<-filtered_fimo_loci
all_fimo_results$loci_name <- gsub("\\([+-]\\)$", "", all_fimo_results$sequence_name)
all_fimo_results_add_features <- all_fimo_results %>%
  left_join(result, by = c("loci_name" = "transcript_id")) %>%
  mutate(TE_subfamily = paste(gene_id, prioritized_feature, sep = "_"))

all_fimo_results<-all_fimo_results_add_features
bed_lookup <- bed_data %>%
  select(name, TE_subfamily, length, chr, start, end) %>%
  rename(loci_name = name, loci_length = length)

fimo_with_length <- all_fimo_results %>%
  left_join(bed_lookup, by = "loci_name", suffix = c("", "_bed"))

fimo_clean <- fimo_with_length %>%
  filter(!is.na(loci_length)) %>%
  select(-TE_subfamily_bed) %>%
  rename(TE_subfamily_orig = TE_subfamily)

merged_data <- merge(bed_data, gft_info, 
                     by.x = c("chr", "start", "end"), 
                     by.y = c("chrom", "start", "end"),
                     all.x = TRUE)
setDT(merged_data)
merged_data[, TE_subfamily := paste(gene_id, prioritized_feature, sep = "_")]
bed_data <- merged_data[, c("chr", "start", "end", "name", "score", "strand.x", "length", "TE_subfamily")]
colnames(bed_data)<-c("chr", "start", "end", "name", "score", "strand", "length", "TE_subfamily")


te_stats <- bed_data %>%
  group_by(TE_subfamily) %>%
  summarise(
    total_loci = n(),                  
    total_length = sum(length),         
    .groups = 'drop'
  )


motif_stats <- fimo_clean %>%
  group_by(TE_subfamily_orig, motif_alt_id) %>%
  summarise(
    motif_instances = n(),                   
    unique_loci_with_motif = n_distinct(loci_name),  
    mean_score = mean(score),                 
    mean_pvalue = mean(p.value),             
    total_loci_length_with_motif = sum(loci_length), 
    .groups = 'drop'
  ) %>%
  rename(TE_subfamily = TE_subfamily_orig)


motif_enrichment <- motif_stats %>%
  left_join(te_stats, by = "TE_subfamily") %>%
  filter(!is.na(total_loci))  

# 计算各种标准化指标
motif_enrichment <- motif_enrichment %>%
  mutate(
    motif_density_per_kb = (motif_instances / total_length) * 1000,

    enrichment_intensity = motif_density_per_kb,
    normalized_count = motif_instances / sqrt(total_length),
    relative_enrichment = motif_density_per_kb / mean(motif_density_per_kb, na.rm = TRUE)
  )


top_motifs_per_te <- motif_enrichment %>%
  group_by(TE_subfamily) %>%
  arrange(desc(enrichment_intensity)) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 10) %>%  
  ungroup()

global_top_motifs <- motif_enrichment %>%
  arrange(desc(enrichment_intensity)) %>%
  head(20)

motif_frequency <- motif_enrichment %>%
  filter(enrichment_intensity > quantile(enrichment_intensity, 0.9, na.rm = TRUE)) %>%
  count(motif_alt_id, sort = TRUE)

te_motif_diversity <- motif_enrichment %>%
  group_by(TE_subfamily) %>%
  summarise(
    total_loci = first(total_loci),
    motif_types = n(),
    max_enrichment = max(enrichment_intensity),
    .groups = 'drop'
  )


top_te_families <- te_motif_diversity %>%
  arrange(desc(max_enrichment)) %>%
  # head(30) %>%  
  pull(TE_subfamily)


important_motifs <- motif_enrichment %>%
  filter(TE_subfamily %in% top_te_families) %>%
  filter(enrichment_intensity > quantile(enrichment_intensity, 0.8, na.rm = TRUE)) %>%
  count(motif_alt_id, sort = TRUE) %>%
  filter(n >= 2) %>%  
  head(100) %>%  
  # head(50) %>%  
  pull(motif_alt_id)

heatmap_data <- motif_enrichment %>%
  filter(TE_subfamily %in% top_te_families,
         motif_alt_id %in% important_motifs) %>%
  dplyr::select(TE_subfamily, motif_alt_id, enrichment_intensity, 
                motif_instances, unique_loci_with_motif, total_loci) %>%
  mutate(
    coverage_percent = (unique_loci_with_motif / total_loci) * 100,
    log_enrichment = log10(enrichment_intensity + 1)
  )
heatmap_data<-heatmap_data[heatmap_data$coverage_percent>=20,]

te_order <- heatmap_data %>%
  group_by(TE_subfamily) %>%
  summarise(mean_enrichment = sum(enrichment_intensity, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(mean_enrichment)) %>%
  pull(TE_subfamily)

motif_order <- heatmap_data %>%
  group_by(motif_alt_id) %>%
  summarise(mean_enrichment = sum(enrichment_intensity, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(mean_enrichment)) %>%
  pull(motif_alt_id)



cell_anno<-read_excel("TF_annotation/DatabaseExtract_v_1.01.xlsx")
cell_anno<-cell_anno[,c("HGNC symbol","DBD")]
KRAB_list<-read_csv("TF_annotation/KRFP_list_combined_uniq_n397.csv")
heatmap_data$motif_alt_id <- as.character(heatmap_data$motif_alt_id)
heatmap_data$TF_type <- NA

krab_tfs_upper <- toupper(KRAB_list$TF)

for(i in 1:nrow(heatmap_data)) {
  if(toupper(heatmap_data$motif_alt_id[i]) %in% krab_tfs_upper) {
    heatmap_data$TF_type[i] <- "KRAB"
  }
}

cell_anno_lookup <- cell_anno
cell_anno_lookup$`HGNC symbol` <- toupper(cell_anno_lookup$`HGNC symbol`)

for(i in 1:nrow(heatmap_data)) {
  if(is.na(heatmap_data$TF_type[i])) { 
    motif_upper <- toupper(heatmap_data$motif_alt_id[i])
    match_idx <- which(cell_anno_lookup$`HGNC symbol` == motif_upper)
    if(length(match_idx) > 0) {
      heatmap_data$TF_type[i] <- cell_anno$DBD[match_idx[1]]
    }
  }
}

table(heatmap_data$TF_type, useNA = "ifany")

heatmap_data$TF_type[heatmap_data$TF_type == "Nuclear receptor"] <- "NR"
heatmap_data$TF_type[heatmap_data$motif_alt_id == "STAT1::STAT2"] <- "STAT"
heatmap_data$TF_type[heatmap_data$motif_alt_id == "Pparg::Rxra"] <- "NR"


heatmap_data$motif_alt_id <- factor(heatmap_data$motif_alt_id, levels = motif_order)
heatmap_data$TE_subfamily <- factor(heatmap_data$TE_subfamily, levels = rev(te_order))
# write_csv(heatmap_data,file="/data2t_2/pathogen_TE_2025_New/08.TF_prediction/Fig3_TE_subfamily_all_fimo_results_coverage_percent_gt20_heatmap_unsort_no_pvalue.csv")

p_heatmap<- ggplot(heatmap_data, aes(x = motif_alt_id, y = TE_subfamily)) +
  geom_point(aes(size = coverage_percent, 
                 fill = log_enrichment),
             shape = 21,
             stroke = 0.1) +      
  scale_size_continuous(name = "Proportion (%)", 
                        range = c(0.5, 3),
                        breaks = c(20, 40, 60, 80,100),
                        limits = c(19, 100)) +
  # scale_fill_gradient(name = "Enrichment Intensity\n(log10)",
  scale_fill_gradient(name = "Enrichment Intensity",
                      low = "#fcf3ef",
                      high = "#6A0B10") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(0.1, "cm"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 7, color = "black",face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    panel.grid.major = element_line(color = "#ebe7e6", size = 0.3),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 7, hjust = 0.5, color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
  ) +
  labs(
    x = "Transcription factor", 
    y = "TE Subfamily"
  )
p_heatmap
# ggsave(p_heatmap,file="Fig3_TE_subfamily_all_fimo_results_coverage_percent_gt20_heatmap_unsort_no_pvalue.pdf",width = 12,height = 4)
################################################################
p_heatmap <- ggplot(heatmap_data, aes(x = motif_alt_id, y = TE_subfamily)) +
  geom_point(aes(size = coverage_percent, 
                 fill = log_enrichment),
             shape = 21,
             stroke = 0.1) +
  scale_size_continuous(name = "Proportion (%)", 
                        range = c(0.5, 2.5),
                        breaks = c(20, 40, 60, 80,100),
                        limits = c(19, 100)) +
  scale_fill_gradient(name = "Enrichment Intensity",
                      low = "#fcf3ef",
                      high = "#6A0B10") +
  facet_grid(. ~ TF_type, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(0.1, "cm"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 7, color = "black",face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    panel.grid.major = element_line(color = "#ebe7e6", size = 0.3),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 7, hjust = 0.5, color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 6, face = "bold", 
                              margin = margin(t = 1.5, b = 1.5, unit = "pt")),  
    strip.background = element_rect(fill = "white", color = "black",size = 0.3),
    strip.placement = "outside",  
    panel.spacing = unit(0, "lines"),
    panel.spacing.x = unit(0, "lines")  
  ) +
  labs(
    x = "Transcription factor", 
    y = "TE Subfamily"
  )

p_heatmap
# ggsave(p_heatmap,file="Fig3_TE_subfamily_all_fimo_results_coverage_percent_gt20_heatmap_add_anno.pdf",width = 12,height = 4)
ggsave(p_heatmap,file="top50_Fig3_TE_subfamily_all_fimo_results_coverage_percent_gt20_heatmap_add_anno.pdf",width = 7,height = 4)
heatmap_data$motif_alt_id <- factor(heatmap_data$motif_alt_id, levels = motif_order)
heatmap_data$TE_subfamily <- factor(heatmap_data$TE_subfamily, levels = rev(te_order))

tf_annotation <- heatmap_data %>%
  select(motif_alt_id, TF_type) %>%
  distinct() %>%
  arrange(motif_alt_id)  

tf_type_colors <- c(
  "KRAB" = "#9E9DCB",         
  "C2H2 ZF" = "#C1747B",       
  "NR" = "#E1C270", 
  "AP-2" = "#929F74",          
  "bHLH" = "#DEA368",         
  "Ets" = "#57B3C3",           
  "EBF1" = "#CF8DB4",          
  "IRF" = "#636363",           
  "Rel" = "#FD79A8",           
  "SMAD" = "#A0E7E5",          
  "STAT" = "#2A4492",          
  "NA" = "#BDBFC3"             
)

p_heatmap <- ggplot(heatmap_data, aes(x = motif_alt_id, y = TE_subfamily)) +
  geom_point(aes(size = coverage_percent, 
                 fill = log_enrichment),
             shape = 21,
             stroke = 0.1) +
  scale_size_continuous(name = "Proportion (%)", 
                        range = c(0.5, 3),
                        breaks = c(20, 40, 60, 80,100),
                        limits = c(19, 100)) +
  scale_fill_gradient(name = "Enrichment Intensity",
                      low = "#fcf3ef",
                      high = "#6A0B10") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(0.1, "cm"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 7, color = "black",face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    panel.grid.major = element_line(color = "#ebe7e6", size = 0.3),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 7, hjust = 0.5, color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 2, r = 5, b = 5, l = 5)  
  ) +
  labs(
    x = "Transcription factor", 
    y = "TE Subfamily"
  )

p_annotation <- ggplot(tf_annotation, aes(x = motif_alt_id, y = 1)) +
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

p_combined <- p_annotation / p_heatmap + 
  plot_layout(heights = c(0.05, 1))  
p_combined
# ggsave(p_combined,file="top50_Fig3_TE_subfamily_all_fimo_results_coverage_percent_gt20_heatmap_add_anno_v2.pdf",width = 7,height = 4.5)
ggsave(p_combined,file="Fig3_TE_subfamily_all_fimo_results_coverage_percent_gt20_heatmap_add_anno_v2.pdf",width = 12,height = 4.5)
########################################################################################################
library(tidyverse)
library(ggplot2)
te_stats <- bed_data %>%
  group_by(name) %>%
  summarise(
    total_loci = n(),                  
    total_length = sum(length),         
    .groups = 'drop'
  )

loci_lengths <- te_stats %>%
  select(name, total_length)
colnames(loci_lengths)<-c("loci_name","total_length")

# motif_position_data <- all_fimo_results %>%
#   left_join(loci_lengths, by = "loci_name") %>%
#   filter(!is.na(total_length)) %>%
#   mutate(
#     motif_center = (start + stop) / 2,
#     relative_position = motif_center / total_length,
#     relative_position = pmax(0, pmin(1, relative_position)),
#     position_bin = case_when(
#       relative_position <= 0.2 ~ "0-20%",
#       relative_position <= 0.4 ~ "20-40%", 
#       relative_position <= 0.6 ~ "40-60%",
#       relative_position <= 0.8 ~ "60-80%",
#       relative_position <= 1.0 ~ "80-100%",
#       TRUE ~ "Other"
#     ),
#     position_bin = factor(position_bin, 
#                           levels = c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%"))
#   )


motif_position_data <- all_fimo_results %>%
  left_join(loci_lengths, by = "loci_name") %>%
  filter(!is.na(total_length)) %>%
  mutate(
    motif_center = (start + stop) / 2,
    relative_position = motif_center / total_length,
    relative_position = pmax(0, pmin(1, relative_position)),
    position_bin = case_when(
      relative_position <= 0.1 ~ "0-10%",
      relative_position <= 0.2 ~ "10-20%",
      relative_position <= 0.3 ~ "20-30%", 
      relative_position <= 0.4 ~ "30-40%",
      relative_position <= 0.5 ~ "40-50%",
      relative_position <= 0.6 ~ "50-60%",
      relative_position <= 0.7 ~ "60-70%",
      relative_position <= 0.8 ~ "70-80%",
      relative_position <= 0.9 ~ "80-90%",
      relative_position <= 1.0 ~ "90-100%",
      TRUE ~ "Other"
    ),
    position_bin = factor(position_bin, 
                          levels = c("0-10%", "10-20%", "20-30%", "30-40%", "40-50%",
                                     "50-60%", "60-70%", "70-80%", "80-90%", "90-100%"))
  )

loci_position_counts <- motif_position_data %>%
  group_by(loci_name, TE_subfamily, position_bin) %>%
  summarise(
    motif_count_in_bin = n(),  
    .groups = 'drop'
  )

loci_total_motifs <- motif_position_data %>%
  group_by(loci_name, TE_subfamily) %>%
  summarise(
    total_motifs = n(),  
    .groups = 'drop'
  )

loci_position_proportions <- loci_position_counts %>%
  left_join(loci_total_motifs, by = c("loci_name", "TE_subfamily")) %>%
  mutate(
    motif_proportion = motif_count_in_bin / total_motifs  # 占比
  ) %>%
  complete(loci_name, position_bin, 
           fill = list(motif_count_in_bin = 0, motif_proportion = 0)) %>%
  group_by(loci_name) %>%
  fill(TE_subfamily, total_motifs, .direction = "downup") %>%
  ungroup()

custom_colors <- c("0-20%" = "#78005e",
                   "20-40%" = "#833c87",
                   "40-60%" = "#8662a1",
                   "60-80%" = "#858abb",
                   "80-100%" ="#b7c8df")
custom_colors = c("0-10%"="#4d004b", 
                  "10-20%"="#810f7c",
                  "20-30%"="#88419d", 
                  "30-40%"="#8c6bb1",
                  "40-50%"="#8c96c6",
                  "50-60%"="#bfd3e6",
                  "60-70%"="#d0e3f0", 
                  "70-80%"="#e0f0f8",  
                  "80-90%"="#f0f8fc",  
                  "90-100%"="#f8fcff")

# loci_position_proportions$position_bin<-factor(loci_position_proportions$position_bin,levels = c( "0-20%",  "20-40%",  "40-60%",  "60-80%", "80-100%" ))
# loci_position_proportions$position_bin <- factor(loci_position_proportions$position_bin,
# levels = c("80-100%", "60-80%", "40-60%", "20-40%", "0-20%"))
p <- ggplot(loci_position_proportions, aes(x = position_bin, y = motif_proportion, fill = position_bin)) +
  # geom_boxplot(alpha = 1, outlier.size = 0.2) +
  geom_boxplot(alpha = 1, outlier.size = 0.1, size = 0.3) +  
  facet_wrap(~ TE_subfamily, scales = "free_x", ncol = 6) + 
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "",
    x = "Sequence Location",
    y = "Proportion(%)",
    fill = "Sequence Location"
  ) +
  # theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 0),
    strip.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.position = "right"
  )

print(p)

# ggsave("TE_subfamily_fimo_results_sequence_location_proportion.pdf", p, width = 16, height = 12)
ggsave("TE_subfamily_fimo_results_sequence_location_proportion_10.pdf", p, width = 18, height = 12)

library(ggplot2)

custom_colors = c("0-10%"="#4d004b", 
                  "10-20%"="#810f7c",
                  "20-30%"="#88419d", 
                  "30-40%"="#8c6bb1",
                  "40-50%"="#8c96c6",
                  "50-60%"="#bfd3e6",
                  "60-70%"="#d0e3f0", 
                  "70-80%"="#e0f0f8",  
                  "80-90%"="#f0f8fc",  
                  "90-100%"="#f8fcff")


color_bar_data <- data.frame(
  segment = c("0-10%", "10-20%", "20-30%", "30-40%", "40-50%","50-60%", "60-70%", "70-80%", "80-90%", "90-100%"),
  x_start = c(0, 1, 2, 3, 4,5,6,7,8,9),
  x_end = c(1, 2, 3, 4, 5,6,7,8,9,10),
  y = 1
)

color_bar_data$segment <- factor(color_bar_data$segment,
                                 levels = c("0-10%", "10-20%", "20-30%", "30-40%", "40-50%","50-60%", "60-70%", "70-80%", "80-90%", "90-100%"))

p_color_bar <- ggplot(color_bar_data, aes(xmin = x_start, xmax = x_end, ymin = 0.5, ymax = 1.5, fill = segment)) +
  geom_rect(color = "white", size = 0.5) +
  scale_fill_manual(values = custom_colors) +
  geom_text(aes(x = (x_start + x_end)/2, y = 1, label = segment), 
            color = "black", size = 2.5, fontface = "bold") +
  scale_x_continuous(limits = c(0, 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0)) +
  labs(
    title = "",
    x = "Sequence Location",
    y = "TE Loci"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    # plot.margin = margin(10, 10, 10, 10),
    axis.text = element_blank(),           
    axis.ticks = element_blank(),          
    panel.grid = element_blank(),          
    panel.border = element_blank(),        
    # axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold",angle = 0),
  )

ggsave("TE_loci_sequence_location_10.pdf", p_color_bar, width = 8, height = 1)
