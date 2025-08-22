# HLA - PEAKS Studio DDA
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(VennDiagram)
library(limma)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggrepel)
library(reshape2)
library(dplyr)
library(tidyr)
library(data.table)
setwd("/data2t_2/pathogen_TE_2025_New/14.Whole_proteome/MHC/result-all-hla-191mb")
metadata <- data.frame(
  PEAKS_ID = c("Area Sample 1","Area Sample 2"),
  # Sample_ID = c("lu109375", "lu109376"),
  Condition = c("Mock","Infection"))

dia_db.peptides<-read_csv("db.peptides.csv")


area_columns <- paste0("Area Sample ", 1:2)
na_counts <- rowSums(is.na(dia_db.peptides[, area_columns]))

dia_db.peptides$Area_Sample_NA <- na_counts
create_te_status <- function(accession_string) {
  if (is.na(accession_string) || accession_string == "") {
    return("Non_TE")
  }
  accession_parts <- strsplit(accession_string, ":")[[1]]
  all_human <- all(grepl("\\|TE_HUMAN$", accession_parts))
  return(ifelse(all_human, "TE", "Non_TE"))
}

dia_db.peptides$TE_Status <- sapply(dia_db.peptides$Accession, create_te_status)
mock_samples <- metadata[metadata$Infection == "Mock", "PEAKS_ID"]
Infection_samples <- metadata[metadata$Infection == "Infection", "PEAKS_ID"]

mock_area_cols <- paste0("Area ", mock_samples)
Infection_area_cols <- paste0("Area ", Infection_samples)

dia_db.peptides$Mean_Mock_Area <- dia_db.peptides$`Area Sample 1`
dia_db.peptides$Mean_Infection_Area <- dia_db.peptides$`Area Sample 2`

te_bed <-fread("/data2t_2/hushu/01.Genomic_features_Gencode_v2/hg38_TE_anno_custom_v20240110_0based_with_feature_summary_v2.bed")

build_te_lookup_table <- function(te_bed_dt) {
  lookup_table <- te_bed_dt[, .(
    TE_Subfamily = paste(unique(gene_id[!is.na(gene_id)]), collapse = ";"),
    TE_Family = paste(unique(family_id[!is.na(family_id)]), collapse = ";"),
    TE_Class = paste(unique(class_id[!is.na(class_id)]), collapse = ";"),
    chrom = paste(unique(chrom[!is.na(chrom)]), collapse = ";"),
    start = paste(unique(start[!is.na(start)]), collapse = ";"),
    end = paste(unique(end[!is.na(end)]), collapse = ";"),
    prioritized_feature = paste(unique(prioritized_feature[!is.na(prioritized_feature)]), collapse = ";")
  ), by = transcript_id]
  
  setkey(lookup_table, transcript_id)
  return(lookup_table)
}

extract_te_info_fast <- function(accession_vector, lookup_table) {
  n <- length(accession_vector)
  
  result_df <- data.frame(
    TE_Name = character(n),
    TE_Subfamily = character(n),
    TE_Family = character(n),
    TE_Class = character(n),
    chrom = character(n),
    start = character(n),
    end = character(n),
    prioritized_feature = character(n),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(n)) {
    accession <- accession_vector[i]
    
    if (is.na(accession) || accession == "") {
      result_df[i, ] <- NA
      next
    }
    
    te_parts <- strsplit(accession, ":")[[1]]
    te_names <- gsub("\\|TE_HUMAN$", "", te_parts)
    te_main_names <- gsub("\\(.*$", "", te_names)
    
    result_df$TE_Name[i] <- paste(te_names, collapse = ";")
    matches <- lookup_table[transcript_id %in% te_main_names]
    
    if (nrow(matches) > 0) {
      result_df$TE_Subfamily[i] <- paste(unique(unlist(strsplit(matches$TE_Subfamily, ";"))), collapse = ";")
      result_df$TE_Family[i] <- paste(unique(unlist(strsplit(matches$TE_Family, ";"))), collapse = ";")
      result_df$TE_Class[i] <- paste(unique(unlist(strsplit(matches$TE_Class, ";"))), collapse = ";")
      result_df$chrom[i] <- paste(unique(unlist(strsplit(matches$chrom, ";"))), collapse = ";")
      result_df$start[i] <- paste(unique(unlist(strsplit(matches$start, ";"))), collapse = ";")
      result_df$end[i] <- paste(unique(unlist(strsplit(matches$end, ";"))), collapse = ";")
      result_df$prioritized_feature[i] <- paste(unique(unlist(strsplit(matches$prioritized_feature, ";"))), collapse = ";")
    } else {
      result_df$TE_Subfamily[i] <- paste(te_main_names, collapse = ";")
      result_df$TE_Family[i] <- NA
      result_df$TE_Class[i] <- NA
      result_df$chrom[i] <- NA
      result_df$start[i] <- NA
      result_df$end[i] <- NA
      result_df$prioritized_feature[i] <- NA
    }
  }
  
  return(result_df)
}


te_bed_dt <- as.data.table(te_bed)
lookup_table <- build_te_lookup_table(te_bed_dt)
te_info_df <- extract_te_info_fast(dia_db.peptides$Accession, lookup_table)

dia_db.peptides$TE_Name <- te_info_df$TE_Name
dia_db.peptides$TE_Subfamily <- te_info_df$TE_Subfamily
dia_db.peptides$TE_Family <- te_info_df$TE_Family
dia_db.peptides$TE_Class <- te_info_df$TE_Class
dia_db.peptides$chrom <- te_info_df$chrom
dia_db.peptides$start <- te_info_df$start
dia_db.peptides$end <- te_info_df$end
dia_db.peptides$prioritized_feature <- te_info_df$prioritized_feature

dia_db.peptides <- dia_db.peptides %>%
  mutate(
    Area.log2FC = log2(
      (coalesce(Mean_Infection_Area, 0) + 0.1) / 
        (coalesce(Mean_Mock_Area, 0) + 0.1)
    )
  )


dia_db.peptides_filtered<-dia_db.peptides[dia_db.peptides$Area_Sample_NA<=1,]
dia_db.peptides_filtered<-dia_db.peptides_filtered[dia_db.peptides_filtered$TE_Status=="TE",]
dia_db.peptides_filtered<-dia_db.peptides_filtered[dia_db.peptides_filtered$Length>=7,]

dia_db.peptides_filtered$Area.log2FC[dia_db.peptides_filtered$Area.log2FC >= 3] <- 3
dia_db.peptides_filtered$Area.log2FC[dia_db.peptides_filtered$Area.log2FC <= -3] <- -3
dia_db.peptides_filtered.hist <- dia_db.peptides_filtered %>% mutate(peptide.logFC = cut(Area.log2FC, seq(from=-4, to=4, by=0.5)))

colors <- c(
  "DNA" = "#f8d196",
  "LINE" = "#d7a9cb", 
  "LTR" = "#6fa4af",
  "SINE" = "#8290bb",
  "Retroposon" = "#277899",
  "DNA;LTR" = "#e5bb94",
  "LINE;DNA" = "#a4c0a2"
)


p1 <- ggplot(dia_db.peptides_filtered.hist, aes(x = peptide.logFC, fill = TE_Class)) + 
  geom_bar(width = 0.9) +  
  scale_fill_manual(name="TE Class", values=colors) + 
  theme_classic() + 
  scale_x_discrete(drop = FALSE,
                   breaks = c("(-3.5,-3]","(0,0.5]", "(0.5,1]","(1.5,2]", "(2.5,3]"),
                   labels = c("Mock only","0.5", "1","2", "Infection only")) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), 
                     limits = c(0, 6)) +  
  xlab('Peptide log2 Fold Change') + 
  ylab('Number of Peptides') + 
  geom_vline(xintercept = c(3, 13), colour="grey", linetype=2, alpha=1) + 
  theme(
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5, 
                               colour = "black"),  
    axis.text.y = element_text(size = 8, colour = "black"),  
    axis.title.x = element_text(size = 10, colour = "black", 
                                margin = margin(t = 10, r = 0, b = 0, l = 0)),  
    axis.title.y = element_text(size = 10, colour = "black"),  
    plot.margin = margin(t=5, r=5, b=5, l=5, unit="mm"),
    plot.title = element_text(size = 10)
  ) +
  ggtitle('Recurrent Up-TE Peptide Infection/Mock')

p1

ggsave(p1,file="MHC_Recurrent_Up_TE_peptide_distribution.pdf",width = 5,height = 3)
write.table(dia_db.peptides, file = "dda_db.peptides_add_anno.txt",quote = FALSE, row.names = FALSE, sep = "\t")
write.table(dia_db.peptides_filtered, file = "dda_db.peptides_filtered_add_anno.txt",quote = FALSE, row.names = FALSE, sep = "\t")


length_TE_Class_table <- dia_db.peptides_filtered.hist %>%
  group_by(Length, TE_Class) %>%
  summarise(count = n(), .groups = 'drop')

p <- ggplot(length_TE_Class_table, aes(x = factor(Length), y = count, fill = TE_Class)) +
  geom_bar(stat = "identity", position = "stack", size = 0.1, width = 0.8) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(expand = c(0.5, 0.05))+
  labs(x = "Peptide Length",
       y = "Count",
       fill = "TE Class") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 7),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.position = "none",
    plot.margin = margin(t=5, r=5, b=5, l=5, unit="mm")
  ) +  
  ggtitle('Recurrent Up-TE Peptide Length Distribution by TE Class') +
  scale_y_continuous(breaks = seq(0, max(aggregate(count ~ Length, length_TE_Class_table, sum)$count), by = 2))
p

ggsave(p,file="MHC_Recurrent_Up-TE_peptide_length_distribution.pdf",width = 2,height = 2)

