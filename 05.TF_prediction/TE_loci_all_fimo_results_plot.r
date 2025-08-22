library(patchwork)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(viridis)
library(data.table)
library(readxl)
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
# 
#   fimo_file <- file.path(dir, "fimo.tsv")
# 
#   if (file.exists(fimo_file)) {
# 
#     fimo_data <- read.delim(fimo_file, stringsAsFactors = FALSE)
# 
#     if (nrow(fimo_data) > 0) {
#       fimo_data$TE_subfamily <- te_name
#       all_fimo_results <- rbind(all_fimo_results, fimo_data)
#     }
#   } else {
#   }
# }

# save(all_fimo_results, file = "all_fimo_results.rdata")
load("all_fimo_results.rdata")

all_fimo_results <- all_fimo_results %>%
  mutate(motif_alt_id = ifelse(
    motif_alt_id == "" | is.na(motif_alt_id),  
    str_extract(motif_id, "^[^.]+"),          
    motif_alt_id                               
  ))
############################################################
#loci
############################################################
loci_info<-read.csv("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_loci/plots/Fig2_All_species_Up_DE-TEs_Loci_padj0.05log2FC1_v2_intergenic_top9.csv",header = T)

te_bed <-fread("/data2t_2/pathogen_TE_2025_New/01.Genomic_features_Gencode/hg38_TE_anno_custom_v20240110_0based_with_feature_summary_v2.bed")

setDT(te_bed)
setDT(loci_info)

loci_info[, c("chrom", "position") := {
  parts <- strsplit(Subfamily, "_")
  chrom_pos <- sapply(parts, function(x) x[length(x)])
  chrom_pos_split <- strsplit(chrom_pos, ":")
  list(
    chrom = sapply(chrom_pos_split, function(x) x[1]),
    position = as.integer(sapply(chrom_pos_split, function(x) x[2]))-1
  )
}]

result <- te_bed[loci_info, on = .(chrom, start = position), nomatch = 0]

loci_df <- result[, .(transcript_id,UniqPathogens)]

filtered_fimo_loci <- all_fimo_results %>%
  mutate(loci_name = gsub("\\([+-]\\)$", "", sequence_name)) %>%
  filter(loci_name %in% loci_df$transcript_id)


all_fimo_results<-filtered_fimo_loci
all_fimo_results$loci_name <- gsub("\\([+-]\\)$", "", all_fimo_results$sequence_name)
bed_lookup <- bed_data %>%
  select(name, TE_subfamily, length, chr, start, end) %>%
  rename(loci_name = name, loci_length = length)

fimo_with_length <- all_fimo_results %>%
  left_join(bed_lookup, by = "loci_name", suffix = c("", "_bed"))


fimo_clean <- fimo_with_length %>%
  filter(!is.na(loci_length)) %>%
  select(-TE_subfamily_bed) %>%
  rename(TE_subfamily_orig = TE_subfamily)

te_stats <- bed_data %>%
  group_by(name) %>%
  summarise(
    total_loci = n(),                    
    total_length = sum(length),         
    .groups = 'drop'
  )
colnames(te_stats)<-c("loci_name","total_loci","total_length")

motif_stats <- fimo_clean %>%
  group_by(loci_name, motif_alt_id) %>%
  summarise(
    motif_instances_per_loci = n(), 
    motif_instances = sum(motif_instances_per_loci),   
    unique_loci_with_motif = n_distinct(loci_name),    
    mean_score = mean(score),                   
    mean_pvalue = mean(p.value),                  
    total_loci_length_with_motif = sum(loci_length),  
    .groups = 'drop'
  ) 

motif_enrichment <- motif_stats %>%
  left_join(te_stats, by = "loci_name") %>%
  filter(!is.na(total_loci))  

motif_enrichment <- motif_enrichment %>%
  mutate(
    motif_density_per_kb = (motif_instances / total_length) * 1000,
    enrichment_intensity = motif_density_per_kb,
    normalized_count = motif_instances / sqrt(total_length),
    relative_enrichment = motif_density_per_kb / mean(motif_density_per_kb, na.rm = TRUE)
  )

te_motif_diversity <- motif_enrichment %>%
  group_by(loci_name) %>%
  summarise(
    total_loci = first(total_loci),
    motif_types = n(),
    max_enrichment = max(enrichment_intensity),
    .groups = 'drop'
  )

top_te_families <- te_motif_diversity %>%
  arrange(desc(max_enrichment)) %>%
  # head(30) %>% 
  pull(loci_name)

important_motifs <- motif_enrichment %>%
  filter(loci_name %in% top_te_families) %>%
  filter(enrichment_intensity > quantile(enrichment_intensity, 0.8, na.rm = TRUE)) %>%
  count(motif_alt_id, sort = TRUE) %>%
  filter(n >= 2) %>% 
  head(100) %>% 
  # head(50) %>%  
  pull(motif_alt_id)

heatmap_data <- motif_enrichment %>%
  filter(loci_name %in% top_te_families,
         motif_alt_id %in% important_motifs) %>%
  dplyr::select(loci_name, motif_alt_id, enrichment_intensity, 
                motif_instances, unique_loci_with_motif, total_loci) %>%
  mutate(
    coverage_percent = (unique_loci_with_motif / total_loci) * 100,
    log_enrichment = log10(enrichment_intensity + 1)
  )
# heatmap_data<-heatmap_data[heatmap_data$coverage_percent>=50,]
TE_loci_name<-result[,c("transcript_id","Subfamily")]
heatmap_data<-merge(heatmap_data,TE_loci_name,by.x = "loci_name",by.y = "transcript_id")
te_order <- heatmap_data %>%
  group_by(Subfamily) %>%
  summarise(mean_enrichment = sum(enrichment_intensity, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(mean_enrichment)) %>%
  pull(Subfamily)

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


# View matching statistics
table(heatmap_data$TF_type, useNA = "ifany")

heatmap_data$TF_type[heatmap_data$TF_type == "Nuclear receptor"] <- "NR"
heatmap_data$TF_type[heatmap_data$motif_alt_id == "STAT1::STAT2"] <- "STAT"
heatmap_data$TF_type[heatmap_data$motif_alt_id == "Pparg::Rxra"] <- "NR"
heatmap_data$TF_type[heatmap_data$motif_alt_id == "Zfp809"] <- "KRAB"

heatmap_data$motif_alt_id <- factor(heatmap_data$motif_alt_id, levels = motif_order)
heatmap_data$Subfamily <- factor(heatmap_data$Subfamily, levels = rev(te_order))
test <- heatmap_data[is.na(heatmap_data$TF_type), ]

# write_csv(heatmap_data,file="/data2t_2/pathogen_TE_2025_New/08.TF_prediction/Fig3_TE_loci_fimo_results_heatmap_unsort_no_pvalue.csv")
p_heatmap<- ggplot(heatmap_data, aes(x = motif_alt_id, y = Subfamily)) +
  geom_point(aes(fill = log_enrichment),shape = 21,stroke = 0.1,size = 1.8) +
  # scale_fill_gradient(name = "Enrichment Intensity\n(log10)",
  scale_fill_gradient(name = "Enrichment Intensity",
                      low = "#fcf3ef",
                      high = "#6A0B10") +
  theme_minimal() +
  theme(
    # Show XY axis lines
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(0.1, "cm"),
    
    # Text settings
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 7, color = "black",face = "bold"),
    
    # Legend settings
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    
    # Grid line settings
    panel.grid.major = element_line(color = "#ebe7e6", size = 0.3),
    panel.grid.minor = element_blank(),
    
    # Title settings
    plot.title = element_text(size = 7, hjust = 0.5, color = "black"),
    
    # Panel background
    panel.background = element_rect(fill = "white", color = NA),
  ) +
  labs(
    x = "Transcription factor", 
    y = "TE Loci"
  )

p_heatmap
# ggsave(p_heatmap,file="Fig3_TE_loci_fimo_results_heatmap_unsort_no_pvalue.pdf",width = 12,height = 4)

################################################################
# Group by TF_type display, remove gaps between facets
p_heatmap<- ggplot(heatmap_data, aes(x = motif_alt_id, y = Subfamily)) +
  geom_point(aes(fill = log_enrichment),shape = 21,stroke = 0.1,size = 1.8) +
  scale_fill_gradient(name = "Enrichment Intensity",
                      low = "#fcf3ef",
                      high = "#6A0B10") +
  # Add facet
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
    # Facet title
    strip.text = element_text(size = 6, face = "bold", 
                              margin = margin(t = 1.5, b = 1.5, unit = "pt")),  # Reduce top and bottom margins
    strip.background = element_rect(fill = "white", color = "black",size = 0.3),
    strip.placement = "outside",  # Place labels outside axis labels
    # Key modification: set panel.spacing to 0
    panel.spacing = unit(0, "lines"),
    panel.spacing.x = unit(0, "lines")  # Ensure x-direction spacing is 0
  ) +
  labs(
    x = "Transcription factor", 
    y = "TE Subfamily"
  )

p_heatmap
# ggsave(p_heatmap,file="top50_Fig3_TE_loci_fimo_results_heatmap_unsort_no_pvalue_add_anno.pdf",width = 7,height = 4)
# ggsave(p_heatmap,file="Fig3_TE_loci_fimo_results_heatmap_unsort_no_pvalue_add_anno.pdf",width = 12,height = 4.3)

# First set factor order
heatmap_data$motif_alt_id <- factor(heatmap_data$motif_alt_id, levels = motif_order)
heatmap_data$Subfamily <- factor(heatmap_data$Subfamily, levels = rev(te_order))

# Create TF type annotation data, maintain original order
tf_annotation <- heatmap_data %>%
  select(motif_alt_id, TF_type) %>%
  distinct() %>%
  arrange(motif_alt_id)  # Arrange by factor order

# Assign colors to different TF types

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
  "NA" = "#BDBFC3",
  "C2H2 ZF; AT hook"="#DBB0D4",
  "Homeodomain"="#A79798"
)
# Create main heatmap
p_heatmap<- ggplot(heatmap_data, aes(x = motif_alt_id, y = Subfamily)) +
  geom_point(aes(fill = log_enrichment),shape = 21,stroke = 0.1,size = 1.8) +
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
    plot.margin = margin(t = 2, r = 5, b = 5, l = 5)  # Increase top margin
  ) +
  labs(
    x = "Transcription factor", 
    y = "TE Subfamily"
  )

# Create top annotation bar
p_annotation <- ggplot(tf_annotation, aes(x = motif_alt_id, y = 1)) +
  geom_tile(aes(fill = TF_type), height = 1) +
  scale_fill_manual(values = tf_type_colors, 
                    name = "TF Type",
                    na.value = "grey80",
                    # Set legend to horizontal layout
                    guide = guide_legend(nrow = 1)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",  # Horizontal direction
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.x = unit(0.2, "cm"),  # Adjust spacing between legend items
    plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
  )

# Combine graphics
p_combined <- p_annotation / p_heatmap + 
  plot_layout(heights = c(0.05, 1))  # Adjust height ratio
p_combined
# ggsave(p_combined,file="top50_Fig3_TE_loci_fimo_results_heatmap_unsort_no_pvalue_add_anno_v2.pdf",width = 7,height = 4.5)
ggsave(p_combined,file="Fig3_TE_loci_fimo_results_heatmap_unsort_no_pvalue_add_anno_v2.pdf",width = 12,height = 4.5)




########################################################################
# TE sequence position motif enrichment analysis - loci level
library(tidyverse)
library(ggplot2)

# 1. Prepare length data
loci_lengths <- te_stats %>%
  select(loci_name, total_length)

# 2. Calculate position intervals for each motif
# motif_position_data <- all_fimo_results %>%
#   left_join(loci_lengths, by = "loci_name") %>%
#   filter(!is.na(total_length)) %>%
#   mutate(
#     # Calculate motif center position
#     motif_center = (start + stop) / 2,
#     # Calculate relative position (between 0-1)
#     relative_position = motif_center / total_length,
#     # Ensure relative position is within reasonable range
#     relative_position = pmax(0, pmin(1, relative_position)),
#     # Assign to 5 intervals
#     position_bin = case_when(
#       relative_position <= 0.2 ~ "0-20%",
#       relative_position <= 0.4 ~ "20-40%", 
#       relative_position <= 0.6 ~ "40-60%",
#       relative_position <= 0.8 ~ "60-80%",
#       relative_position <= 1.0 ~ "80-100%",
#       TRUE ~ "Other"
#     ),
#     # Set factor order
#     position_bin = factor(position_bin, 
#                           levels = c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%"))
#   )

motif_position_data <- all_fimo_results %>%
  left_join(loci_lengths, by = "loci_name") %>%
  filter(!is.na(total_length)) %>%
  mutate(
    # Calculate motif center position
    motif_center = (start + stop) / 2,
    # Calculate relative position (between 0-1)
    relative_position = motif_center / total_length,
    # Ensure relative position is within reasonable range
    relative_position = pmax(0, pmin(1, relative_position)),
    # Assign to 10 intervals
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
    # Set factor order
    position_bin = factor(position_bin, 
                          levels = c("0-10%", "10-20%", "20-30%", "30-40%", "40-50%",
                                     "50-60%", "60-70%", "70-80%", "80-90%", "90-100%"))
  )


# 3. Calculate motif count for each loci in each position interval
loci_position_counts <- motif_position_data %>%
  group_by(loci_name, TE_subfamily, position_bin) %>%
  summarise(
    motif_count_in_bin = n(),  # Motif count in this position interval
    .groups = 'drop'
  )

# 4. Calculate total motif count for each loci
loci_total_motifs <- motif_position_data %>%
  group_by(loci_name, TE_subfamily) %>%
  summarise(
    total_motifs = n(),  # Total motif count for this loci
    .groups = 'drop'
  )

# 5. Calculate motif proportion for each loci in each position interval
loci_position_proportions <- loci_position_counts %>%
  left_join(loci_total_motifs, by = c("loci_name", "TE_subfamily")) %>%
  mutate(
    motif_proportion = motif_count_in_bin / total_motifs  # Proportion
  ) %>%
  # Ensure each loci has data for 5 position intervals (fill 0 for intervals without motifs)
  complete(loci_name, position_bin, 
           fill = list(motif_count_in_bin = 0, motif_proportion = 0)) %>%
  # Refill TE_subfamily and total_motifs information
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
loci_position_proportions$position_bin <- factor(
  loci_position_proportions$position_bin,
  levels = rev(levels(loci_position_proportions$position_bin))
)


# 然后绘图
p_loci <- ggplot(loci_position_proportions, aes(x = loci_name, y = motif_proportion, fill = position_bin)) +
  geom_col(position = "stack", width = 0.8, color = "white", size = 0.2) +
  scale_fill_manual(values = custom_colors, name = "Sequence Location") +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"),
                     breaks = seq(0, 1, 0.2),
                     expand = c(0, 0),
                     limits = c(0, 1)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 10, face = "bold"),
    legend.position = "none",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", size = 0.3),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "",
    x = "TE Loci",
    y = "Proportion(%)",
  )
p_loci
# ggsave(p_loci,file="TE_loci_fimo_results_sequence_location_proportion_10.pdf",width = 6,height = 4)

# loci_position_proportions$position_bin<-factor(loci_position_proportions$position_bin,levels = c( "0-20%",  "20-40%",  "40-60%",  "60-80%", "80-100%" ))

loci_position_proportions$position_bin<-factor(loci_position_proportions$position_bin,levels = c("0-10%", "10-20%", "20-30%", "30-40%", "40-50%",
                                                                                                 "50-60%", "60-70%", "70-80%", "80-90%", "90-100%"))

p_distribution <- ggplot(loci_position_proportions, aes(x = position_bin, y = motif_proportion, fill = position_bin)) +
  geom_boxplot(alpha = 1, outlier.size = 0.1, size = 0.3) +  
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "",
    x = "Sequence Location",
    y = "Proportion(%)"
  ) +
  # theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)  # Add this line
  )

print(p_distribution)
library(patchwork)
p <- p_loci + p_distribution + 
  plot_layout(widths = c(2, 1))
p
ggsave(p,file="TE_loci_fimo_results_sequence_location_proportion_distribution_10.pdf",width = 8,height = 4)
