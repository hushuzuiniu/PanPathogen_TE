library(GenomicFeatures)
library(GenomicRanges)

setwd("/data2t_2/pathogen_TE_2025_New/01.Genomic_features_Gencode/")
gtf_file <- "/data2t_2/pathogen_TE_2025_New/01.Genomic_features_Gencode/gencode.v26.primary_assembly.annotation.filtered.gtf"
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

chrom_sizes_url <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
chrom_sizes <- read.table(chrom_sizes_url, header = FALSE, stringsAsFactors = FALSE)
colnames(chrom_sizes) <- c("chrom", "size")

valid_chrom <- grep("^chr[0-9XYM]+$", chrom_sizes$chrom, value = TRUE)
chrom_sizes <- chrom_sizes[chrom_sizes$chrom %in% valid_chrom, ]

# 3.1 (Genes)
gene_regions <- genes(txdb)
gene_regions <- gene_regions[seqnames(gene_regions) %in% valid_chrom]
gene_regions <- reduce(gene_regions,ignore.strand=TRUE)
gene_df<-as.data.frame(gene_regions)
gene_length <- sum(width(gene_regions))

# 3.2 (Exons)
exon_regions <- exons(txdb)
exon_regions <- exon_regions[seqnames(exon_regions) %in% valid_chrom]
exon_regions <- reduce(exon_regions,ignore.strand=TRUE)
# exon_regions <- reduce(exon_regions)
exon_df<-as.data.frame(exon_regions)
exon_length <- sum(width(exon_regions))

# 3.3 (introns)
introns_regions <- setdiff(gene_regions, exon_regions)
introns_length <- sum(width(introns_regions))
introns_df<-as.data.frame(introns_regions)

genome_ranges <- GRanges(
  seqnames = chrom_sizes$chrom,
  ranges = IRanges(start = 1, end = chrom_sizes$size),
  strand = "*"
)
intergenic_regions<-setdiff(genome_ranges,gene_regions)
intergenic_length <- sum(width(intergenic_regions))
intergenic_df<-as.data.frame(intergenic_regions)
total_genome_length <- sum(chrom_sizes$size)

#########################################################################################################
add_type <- function(gr, type) {
  mcols(gr)$feature_type <- type
  return(gr)
}

intron_typed <- add_type(introns_regions, "intron")
intergenic_typed <- add_type(intergenic_regions, "intergenic")
exon_typed <- add_type(exon_regions, "exon")
all_regions <- c(intron_typed, exon_typed,intergenic_typed)


all_df <- data.frame(
  chrom = seqnames(all_regions),
  start = start(all_regions)-1,
  end = end(all_regions),
  width = width(all_regions),
  strand = strand(all_regions),
  feature_type = mcols(all_regions)$feature_type
)

# write.table(all_df, file = "all_genomic_regions_v2.bed", sep = "\t", quote = FALSE, row.names = FALSE)
#################################################################
library(ggplot2)
library(dplyr)
introns_length <- introns_length
intergenic_length <- intergenic_length
exon_length <- exon_length
total_genome_length <- total_genome_length

introns_percent <- introns_length/total_genome_length * 100
intergenic_percent <- intergenic_length/total_genome_length * 100
exon_percent<-exon_length/total_genome_length * 100

genomic_features <- data.frame(
  feature = c("exon", "intron", "intergenic"),
  percentage = c(exon_percent,introns_percent,intergenic_percent)
)

genomic_features$percentage <- round(genomic_features$percentage, 1)

genomic_features <- genomic_features %>%
  arrange(desc(feature)) %>%
  mutate(
    ypos = cumsum(percentage) - 0.5 * percentage,
    label = paste0(feature, "\n", percentage, "%")
  )


genomic_features <- genomic_features %>%
  arrange(percentage)

feature_colors <- c(
  "CDS" = "#8559A5",        
  "3'UTR" = "#D52928",      
  "5'UTR" = "#249E68",    
  "exon" = "#8C574C",
  "intron" = "#1F78B4",       
  "intergenic" = "#F57E20"    
)

p_exploded <- ggplot(genomic_features, aes(x = "", y = percentage, fill = feature)) +
  geom_bar(stat = "identity", width = 0.9) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = feature_colors) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  labs(
    title = "Genomic features in the Gencode",
    subtitle = "human transcriptome (Hg38 v43)"
  ) +
  geom_text(
    aes(y = ypos, label = label),
    size = 4
  )


print(p_exploded)
# ggsave("Genomic_features_Gencode_v2.pdf", p_exploded, width = 4, height = 4)

feature_colors <- c(
  "exon" = "#B6ABCC",
  "intron" = "#CCD9EC",       
  "intergenic" = "#6FA4AF"    
)
p_donut <- ggplot(genomic_features, aes(x = 2, y = percentage, fill = feature)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = feature_colors) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  labs(
    title = "Genomic features in the Gencode",
    subtitle = "human transcriptome (Hg38 v43)"
  ) +
  geom_text(
    aes(y = ypos, label = label),
    size = 4
  ) +
  xlim(0.5, 2.5) 
p_donut
ggsave("FigS1A_Genomic_features_Gencode_v2.pdf", p_donut, width = 4, height = 4)
