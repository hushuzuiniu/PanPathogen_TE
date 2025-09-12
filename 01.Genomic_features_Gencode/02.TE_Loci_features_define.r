library(data.table)
library(ggplot2)
library(stringr)
library(VennDiagram)
library(dplyr)
setwd("/data2t_2/hushu/01.Genomic_features_Gencode_v2/")
# bedtools intersect -a hg38_TE_anno_custom_v20240110_0based.bed -b all_genomic_regions_v2.bed -wa -wb > hg38_TE_anno_custom_v20240110_0based_with_feature_v2.bed
intersections <- fread("hg38_TE_anno_custom_v20240110_0based_with_feature_v2.bed", header=FALSE)

colnames(intersections) <- c(
  "TE_chrom", "TE_start", "TE_end", "TE_description", "TE_dot", "TE_strand",
  "Gene_chrom", "Gene_start", "Gene_end","Gene_width","Gene_strand","feature_type"
)

intersections[, te_id := paste(TE_chrom, TE_start, TE_end, TE_description, sep="_")]

result <- intersections[, .(
  chrom = first(TE_chrom),
  start = first(TE_start),
  end = first(TE_end),
  description = first(TE_description),
  dot = first(TE_dot),
  strand = first(TE_strand),
  feature_type_annotation = paste(unique(feature_type), collapse=",")
), by = te_id]

result[, te_id := NULL]

########################################################################################################
# Priority order: CDS > 5UTR > 3UTR > noncoding_exon > intron > intergenic
prioritize_feature <- function(features) {
  features_vec <- unlist(strsplit(features, ","))
  if (any(grepl("exon", features_vec))) {
    return("exon")
  } else if (any(grepl("intron", features_vec))) {
    return("intron")
  } else {
    return("intergenic")
  }
}

result[, prioritized_feature := sapply(feature_type_annotation, prioritize_feature)]

# Extract class_id and family_id from description
extract_info <- function(desc) {
  # Extract gene_id
  gene_match <- str_match(desc, "gene_id \"([^\"]+)\"")
  gene_id <- ifelse(!is.na(gene_match[1]), gene_match[2], NA)
  
  # Extract transcript_id
  transcript_match <- str_match(desc, "transcript_id \"([^\"]+)\"")
  transcript_id <- ifelse(!is.na(transcript_match[1]), transcript_match[2], NA)
  
  # Extract class_id
  class_match <- str_match(desc, "class_id \"([^\"]+)\"")
  class_id <- ifelse(!is.na(class_match[1]), class_match[2], NA)
  
  # Extract family_id
  family_match <- str_match(desc, "family_id \"([^\"]+)\"")
  family_id <- ifelse(!is.na(family_match[1]), family_match[2], NA)
  
  return(data.frame(gene_id=gene_id,transcript_id=transcript_id,class_id = class_id, family_id = family_id))
}

# Extract class_id and family_id
info_df <- as.data.frame(t(sapply(result$description, extract_info)))
result$gene_id<-info_df$gene_id
result$transcript_id<-info_df$transcript_id
result$class_id <- info_df$class_id
result$family_id <- info_df$family_id

fwrite(result, "hg38_TE_anno_custom_v20240110_0based_with_feature_summary_v2.bed", sep="\t", quote=FALSE)

result$gene_id <- as.character(result$gene_id)
result$transcript_id <- as.character(result$transcript_id)
result$class_id <- as.character(result$class_id)
result$prioritized_feature <- as.character(result$prioritized_feature)

plot_data <- as.data.frame(table(result$class_id, result$prioritized_feature))
colnames(plot_data) <- c("class_id", "prioritized_feature", "N")
plot_data <- plot_data[plot_data$class_id != "NA" & plot_data$N > 0, ]


plot_data$prioritized_feature <- factor(
  plot_data$prioritized_feature, 
  levels = c("exon", "intron", "intergenic")
  # levels = c("CDS", "5UTR", "3UTR", "noncoding_exon", "intron", "intergenic")
)

p <- ggplot(plot_data, aes(x = class_id, y = N, fill = prioritized_feature)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = N), position = position_dodge(width = 0.9), 
            vjust = -0.5, size =3) + 
  theme_classic() +
  labs(
    title = "Distribution of Gencode Features by TE Class",
    x = "TE Class",
    y = "Count",
    fill = "Feature"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linetype = "dotted", color = "gray80"),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.line = element_line(color = "black"),
  ) +
  scale_fill_manual(values = c(
    # "CDS" = "#8559A5",
    # "3UTR" = "#D52928",
    # "5UTR" = "#249E68",
    "exon" = "#8C574C",
    "intron" = "#1F78B4",
    "intergenic" = "#F57E20"
  )) +
  scale_x_discrete(expand = expansion(mult = c(0.13, 0.13))) +
  geom_vline(xintercept = seq(1.5, length(unique(plot_data$class_id))-0.5), 
             linetype = "dotted", color = "gray70")

print(p)
ggsave("FigureS1B_Gencode_features_distribution_TE_class_v2.pdf", p, width = 12, height = 7)
########################################################################################################
result<-fread("hg38_TE_anno_custom_v20240110_0based_with_feature_summary_v2.bed")
TE_Loci_info <- result[, c("gene_id", "transcript_id", "class_id", "family_id", "prioritized_feature")]
TE_Loci_info$region <- TE_Loci_info$prioritized_feature
TE_Loci_info$region[TE_Loci_info$prioritized_feature %in% 
                      c("3UTR", "5UTR", "CDS", "noncoding_exon")] <- "exon"

region_counts <- table(TE_Loci_info$region)
print(region_counts)
gene_region <- unique(TE_Loci_info[, c("gene_id", "region")])

gene_locations <- gene_region %>%
  group_by(gene_id) %>%
  summarize(
    in_exon = any(region == "exon"),
    in_intron = any(region == "intron"),
    in_intergenic = any(region == "intergenic")
  )

venn_counts <- gene_locations %>%
  summarize(
    exon_only = sum(in_exon & !in_intron & !in_intergenic),
    intron_only = sum(!in_exon & in_intron & !in_intergenic),
    intergenic_only = sum(!in_exon & !in_intron & in_intergenic),
    exon_intron = sum(in_exon & in_intron & !in_intergenic),
    exon_intergenic = sum(in_exon & !in_intron & in_intergenic),
    intron_intergenic = sum(!in_exon & in_intron & in_intergenic),
    all_three = sum(in_exon & in_intron & in_intergenic)
  )

print(venn_counts)
# # A tibble: 1 × 7
# exon_only intron_only intergenic_only exon_intron exon_intergenic intron_intergenic all_three
# <int>       <int>           <int>       <int>           <int>             <int>     <int>
#   0         0           1               0           2               27              1043
venn_colors <- c("#7FB2D5", "#F47F72", "#BFBCDA")
category_names <- c("exon", "intron", "intergenic")

exon_only <- 0
intron_only <- 0
intergenic_only <- 1
exon_intron <- 0
exon_intergenic <- 2
intron_intergenic <- 27
all_three <- 1043

area1 <- exon_only + exon_intron + exon_intergenic + all_three  # exon 
area2 <- intron_only + exon_intron + intron_intergenic + all_three  # intron 
area3 <- intergenic_only + exon_intergenic + intron_intergenic + all_three  # intergenic 
n12 <- exon_intron + all_three  # exon  intron 
n23 <- intron_intergenic + all_three  # intron  intergenic 
n13 <- exon_intergenic + all_three  # exon  intergenic 
n123 <- all_three  

venn.plot <- draw.triple.venn(
  area1 = area1,         # exon 
  area2 = area2,         # intron 
  area3 = area3,         # intergenic 
  n12 = n12,             # exon & intron
  n23 = n23,             # intron & intergenic
  n13 = n13,             # exon & intergenic
  n123 = n123,           # exon & intron & intergenic
  category = category_names,
  cat.cex = 1.5,
  cat.fontface = "bold", 
  cat.col = c("#7FB2D5", "#F47F72", "#BFBCDA"),
  cat.pos = c(0, 0, 180),
  cat.dist = c(0.05, 0.05, 0.05),
  
  fill = venn_colors,
  alpha = 1,
  cex = 1.2,
  fontface = "bold",
  lwd = 2,
  euler.d = FALSE, 
  scaled = FALSE, 
  main = "TE subfamilies",
  main.cex = 2,
  main.fontface = "bold"
)

grid::grid.newpage()
grid::grid.draw(venn.plot)
pdf("FigureS1C_TE_subfamil_venn_diagram_v2.pdf", width = 4, height = 4)
grid::grid.newpage()
grid::grid.draw(venn.plot)
dev.off()



















