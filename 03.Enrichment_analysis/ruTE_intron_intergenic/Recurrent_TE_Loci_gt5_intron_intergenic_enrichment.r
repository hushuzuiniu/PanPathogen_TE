library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)  
library(enrichplot)
library(ggplot2)
library(DOSE)
library(ggnewscale)
library(stringr)
library(gground)
library(ggprism)
library(tidyverse)
library(patchwork)
library(gridExtra)
library(dplyr)
##################################################################
setwd("/data2t_2/pathogen_TE_2025_New/03.Enrichment_analysis/intron_intergenic/")
recurrent_up_intron_intergenic_TE_Loci<-read.csv("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_loci/plots/Recurrent_Up_TE_Loci_intergenic_intron_UniqPathogens_gt1.csv")
table(recurrent_up_intron_intergenic_TE_Loci$UniqPathogens)
gt5_Up_TE_Loci<-recurrent_up_intron_intergenic_TE_Loci[recurrent_up_intron_intergenic_TE_Loci$UniqPathogens>=5,]
te_bed <-fread("/data2t_2/pathogen_TE_2025_New/01.Genomic_features_Gencode/hg38_TE_anno_custom_v20240110_0based_with_feature_summary_v2.bed")

setDT(te_bed)
setDT(gt5_Up_TE_Loci)

gt5_Up_TE_Loci[, c("chrom", "position") := {
  parts <- strsplit(Subfamily, "_")
  chrom_pos <- sapply(parts, function(x) x[length(x)])  
  chrom_pos_split <- strsplit(chrom_pos, ":")
  list(
    chrom = sapply(chrom_pos_split, function(x) x[1]),
    position = as.integer(sapply(chrom_pos_split, function(x) x[2]))-1
  )
}]

result <- te_bed[gt5_Up_TE_Loci, on = .(chrom, start = position), nomatch = 0]

final_result <- result[, .(chrom, start, end, transcript_id, UniqPathogens, Subfamily, 
                           class_id, family_id, feature_type_annotation)]
write.table(final_result, file = "/data2t_2/pathogen_TE_2025_New/03.Enrichment_analysis/intron_intergenic/Recurrent_TE_Loci_gt5_intergenic_intron_n5481.bed", row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t")

##################bedtools linux########################################################################
cd /data2t_2/pathogen_TE_2025_New/03.Enrichment_analysis/intron_intergenic
bedtools sort -g /data2t_2/hushu/00.ref/hg38.p13.fa.fai -i Recurrent_TE_Loci_gt5_intergenic_intron_n5481.bed > Recurrent_TE_Loci_gt5_intergenic_intron_n5481.genome.sorted.bed
grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y)\s' /data2t_2/pathogen_TE_2025_New/01.Genomic_features_Gencode/genes_regions.bed > /data2t_2/pathogen_TE_2025_New/01.Genomic_features_Gencode/genes_regions.standard_chr.bed
bedtools sort -g /data2t_2/hushu/00.ref/hg38.p13.fa.fai -i /data2t_2/pathogen_TE_2025_New/01.Genomic_features_Gencode/genes_regions.standard_chr.bed > /data2t_2/pathogen_TE_2025_New/01.Genomic_features_Gencode/genes_regions.standard_chr.sorted.bed

bedtools closest -a Recurrent_TE_Loci_gt5_intergenic_intron_n5481.genome.sorted.bed -b /data2t_2/pathogen_TE_2025_New/01.Genomic_features_Gencode/genes_regions.standard_chr.sorted.bed -d -t first -k 1 > Recurrent_TE_Loci_gt5_intergenic_intron_n5481_closest_genes.bed
############################################################################################################
TE_closest_genes<-read.csv("Recurrent_TE_Loci_gt5_intergenic_intron_n5481_closest_genes.bed",sep = "\t",header = F)
colnames(TE_closest_genes) <- c(
  "chr1", "start", "end", "TE_id", "UniqPathogens","TE_location" ,"class", "family","feature",
  "chr2", "gene_start", "gene_end", "gene_width", "strand",  "gene_info","distance"
)
################################################################
gene_list <- unique(TE_closest_genes$gene_info)
gene_gtf<-fread("/data2t_2/hushu/00.ref/hg38.p13.gene.anno.gtf")
gene_gtf_filtered <- gene_gtf[gene_gtf$V3 == "gene", ]
extract_id_name <- function(attribute_string) {
  gene_id_match <- str_match(attribute_string, 'gene_id "([^"]+)"')[, 2]
  gene_name_match <- str_match(attribute_string, 'gene_name "([^"]+)"')[, 2]
  return(data.frame(
    gene_id = gene_id_match,
    gene_name = gene_name_match,
    stringsAsFactors = FALSE
  ))
}
id_name_list <- lapply(gene_gtf_filtered$V9, extract_id_name)
id_name_df <- do.call(rbind, id_name_list)
id_name_df <- unique(id_name_df)
head(id_name_df)
gene_list<-as.data.frame(gene_list)
result <- merge(gene_list, id_name_df, by.x = "gene_list", by.y = "gene_id", all.x = TRUE)

gene_ids <- bitr(result$gene_name, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
entrez_ids <- gene_ids$ENTREZID

go_bp <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",  
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

go_mf <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",  
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

go_cc <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",  
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

kegg_result <- enrichKEGG(gene = entrez_ids,
                          organism = "hsa",  
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
kegg_result <- setReadable(kegg_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

save_enrichment_results <- function(enrichment_obj, file_name) {
  if (is.null(enrichment_obj) || nrow(enrichment_obj@result) == 0) {
    cat(paste0("No significant enrichment found for ", file_name, "\n"))
    return(NULL)
  }
  write.csv(enrichment_obj@result, paste0(file_name, ".csv"), row.names = FALSE)
  return(enrichment_obj@result)
}
go_bp_df <- save_enrichment_results(go_bp, "Recurrent_TE_Loci_recurrent_up_intron_intergenic_TE_Loci_closest_genes_GOBP_enrichment")
go_mf_df <- save_enrichment_results(go_mf, "Recurrent_TE_Loci_recurrent_up_intron_intergenic_TE_Loci_closest_genes_GOMF_enrichment")
go_cc_df <- save_enrichment_results(go_cc, "Recurrent_TE_Loci_recurrent_up_intron_intergenic_TE_Loci_closest_genes_GOCC_enrichment")
kegg_df <- save_enrichment_results(kegg_result, "Recurrent_TE_Loci_recurrent_up_intron_intergenic_TE_Loci_closest_genes_KEGG_enrichment")


prepare_go_data <- function(go_bp, go_cc, go_mf, top_n = 10) {
  if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
    bp_data <- go_bp@result %>%
      filter(p.adjust < 0.05) %>%
      arrange(p.adjust) %>%
      head(top_n) %>%
      mutate(Category = "BP")
  } else {
    bp_data <- NULL
  }
  
  if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
    cc_data <- go_cc@result %>%
      filter(p.adjust < 0.05) %>%
      arrange(p.adjust) %>%
      head(top_n) %>%
      mutate(Category = "CC")
  } else {
    cc_data <- NULL
  }
  
  if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
    mf_data <- go_mf@result %>%
      filter(p.adjust < 0.05) %>%
      arrange(p.adjust) %>%
      head(top_n) %>%
      mutate(Category = "MF")
  } else {
    mf_data <- NULL
  }
  
  combined_data <- bind_rows(bp_data, cc_data, mf_data)
  if (nrow(combined_data) == 0) {
    return(NULL)
  }
  combined_data$neg_log10_padj <- -log10(combined_data$p.adjust)
  
  combined_data$GeneRatio_numeric <- sapply(combined_data$GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
  
  combined_data$Category <- factor(combined_data$Category, levels = c("BP", "CC", "MF"))
  combined_data <- combined_data %>%
    group_by(Category) %>%
    mutate(Description = factor(Description, levels = rev(Description))) %>%
    ungroup()
  
  return(combined_data)
}

combined_go_data <- prepare_go_data(go_bp, go_cc, go_mf, top_n = 10)
########################################################################
process_go_result <- function(go_result, ontology, gene_ids) {
  if (is.null(go_result) || nrow(go_result@result) == 0) {
    return(NULL)
  }
  
  result <- go_result@result
  result$ONTOLOGY <- ontology
  id_to_symbol <- setNames(gene_ids$SYMBOL, gene_ids$ENTREZID)
  result$geneID <- sapply(strsplit(result$geneID, "/"), function(ids) {
    valid_ids <- ids[ids %in% names(id_to_symbol)]
    symbols <- id_to_symbol[valid_ids]
    paste(symbols, collapse = ", ")
  })
  
  return(result)
}

process_kegg_result <- function(kegg_result, gene_ids) {
  if (is.null(kegg_result) || nrow(kegg_result@result) == 0) {
    return(NULL)
  }
  result <- kegg_result@result
  id_to_symbol <- setNames(gene_ids$SYMBOL, gene_ids$ENTREZID)
  result$geneID <- sapply(strsplit(result$geneID, "/"), function(ids) {
    valid_ids <- ids[ids %in% names(id_to_symbol)]
    symbols <- id_to_symbol[valid_ids]
    paste(symbols, collapse = ", ")
  })
  
  return(result)
}

bp_data <- process_go_result(go_bp, "BP", gene_ids)
cc_data <- process_go_result(go_cc, "CC", gene_ids)
mf_data <- process_go_result(go_mf, "MF", gene_ids)
kegg_data <- process_kegg_result(kegg_result, gene_ids)

GO <- rbind(bp_data, cc_data, mf_data)
KEGG <- kegg_data

if (is.null(KEGG) || nrow(KEGG) == 0) {
  warning("")
  KEGG <- data.frame(
    ID = character(0),
    Description = character(0),
    GeneRatio = character(0),
    BgRatio = character(0),
    pvalue = numeric(0),
    p.adjust = numeric(0),
    qvalue = numeric(0),
    Count = integer(0),
    geneID = character(0)
  )
}

use_pathway <- GO %>%
  filter(p.adjust < 0.05) %>%
  group_by(ONTOLOGY) %>%
  arrange(p.adjust, desc(Count)) %>%
  slice_head(n = 5) %>%
  rbind(
    if(nrow(KEGG) > 0) {
      KEGG %>%
        filter(p.adjust < 0.05) %>%
        arrange(p.adjust, desc(Count)) %>%
        slice_head(n = 5) %>%
        mutate(ONTOLOGY = 'KEGG')
    }
  ) %>%
  ungroup() %>%
  mutate(ONTOLOGY = factor(ONTOLOGY,
                           levels = rev(c('BP','CC','MF','KEGG')))) %>%
  dplyr::arrange(ONTOLOGY, p.adjust) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  tibble::rowid_to_column('index')
simplify_gene_list <- function(gene_string, max_genes = 25) {
  genes <- unlist(strsplit(gene_string, ", "))
  if (length(genes) <= max_genes) {
    return(gene_string)
  } else {
    return(paste0(paste(genes[1:max_genes], collapse = ", "), " etc"))
  }
}

use_pathway$geneID_simplified <- sapply(use_pathway$geneID, simplify_gene_list)

width <- 0.5
xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1

rect.data <- group_by(use_pathway, ONTOLOGY) %>%
  reframe(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -3* width,
    xmax = -2* width,
    ymax = cumsum(n),
    ymin = lag(ymax, default =0) +0.6,
    ymax = ymax +0.4
  )


pal <- c("BP" = "#368abf",  "CC" = "#cac54d","MF" = "#40a647","KEGG" = "#c7232c" )

enrichment_plot <- use_pathway %>%
  ggplot(aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
  geom_round_col(
    aes(y = Description), width = 0.7, alpha = 0.8) +
  geom_text(
    aes(x = 0.1, label = Description),
    hjust = 0, size = 4) +
  geom_text(
    aes(x = 0.3, label = geneID_simplified, colour = ONTOLOGY),
    hjust = 0, vjust = 2.8, size = 2,show.legend = FALSE,fontface = 'italic') +
  geom_point(
    aes(x = -width, size = Count, fill = ONTOLOGY),shape = 21,stroke = 0.5) +
  geom_text(
    aes(x = -width, label = Count),size = 3) + scale_size_continuous(name = 'Count', range = c(2, 8)) +
  geom_round_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,fill = ONTOLOGY),
    data = rect.data,radius = unit(1,'mm'),inherit.aes =FALSE) +
  geom_text(
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
    data = rect.data,inherit.aes = FALSE,fontface = "bold",size = 4.5,color = "white",angle = 90) +
  geom_segment(
    aes(x = 0, y = 0, xend = 0, yend = max(index) + 1),linewidth = 1,inherit.aes = FALSE) +
  geom_segment(
    aes(x = 0, y = 0, xend = xaxis_max, yend = 0),linewidth = 1,inherit.aes = FALSE) +
  labs(y = NULL) +
  scale_fill_manual(name = 'Category', values = pal) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(
    breaks = seq(0, xaxis_max, 5),
    limits = c(-3 * width, xaxis_max),  
    expand = expansion(c(0, 0))
  ) +
  theme_prism() +
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right",
    legend.title = element_text(),
    # plot.title = element_text(hjust = 0.5, size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
enrichment_plot

########################################################################
########################################################################
process_enrichment <- function(enrich_result, top_n = 20) {
  result_df <- as.data.frame(enrich_result)
  significant_paths <- result_df %>% 
    filter(p.adjust <= 0.05) %>%
    arrange(pvalue) %>%
    head(top_n)
  significant_paths <- significant_paths[order(significant_paths$pvalue, decreasing = TRUE), ]
  
  significant_paths$neg_log10_p <- -log10(significant_paths$p.adjust)
  return(significant_paths)
}

plot_enrichment_barplot <- function(data, title, bar_color = "lightpink", x_lab = "-log10(pvalue)", count_lab = "Count", show_legend = TRUE) {
  
  if(!"Count" %in% colnames(data)) {
    stop(" 'Count' ")
  }
  
  data <- data[order(data$neg_log10_p), ]
  
  data$Description <- factor(data$Description, levels = data$Description)
  
  max_neg_log10 <- max(data$neg_log10_p, na.rm = TRUE)
  max_count <- max(data$Count, na.rm = TRUE)
  min_count <- 0  
  
  if (max_neg_log10 <= 10) {
    pvalue_breaks <- seq(0, ceiling(max_neg_log10) + 1, by = 2)
  } else if (max_neg_log10 <= 20) {
    pvalue_breaks <- seq(0, ceiling(max_neg_log10) + 2, by = 5)
  } else {
    pvalue_breaks <- seq(0, ceiling(max_neg_log10) + 5, by = 10)
  }
  pvalue_max <- max(pvalue_breaks)
  
  count_breaks <- pretty(c(0, max_count), n = 5)
  count_max <- max(count_breaks)
  count_min <- 0
  
  p <- ggplot(data, aes(y = Description)) +
    geom_bar(aes(x = neg_log10_p), 
             stat = "identity", 
             fill = bar_color, 
             alpha = 0.6,
             width = 0.8) +
    
    geom_path(aes(x = Count * pvalue_max / count_max, y = as.numeric(Description)), 
              color = "black", 
              size = 0.5) +
    geom_point(aes(x = Count * pvalue_max / count_max), 
               color = "black", 
               size = 1) +
    
    scale_x_continuous(
      name = x_lab,
      position = "top",
      breaks = pvalue_breaks,
      limits = c(0, pvalue_max),
      expand = expansion(mult = c(0, 0.05)),
      sec.axis = sec_axis(
        trans = ~ . * count_max / pvalue_max,
        name = count_lab,
        breaks = count_breaks
      )
    ) +
    
    scale_y_discrete(expand = expansion(mult = c(0.02, 0.02))) +
    
    labs(
      title = title,
      y = ""
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(size = 10, hjust = 1, color = "black"),
      axis.text.x.top = element_text(size = 10, color = "black"),
      axis.text.x.bottom = element_text(size = 10, color = "black"),
      axis.title.x.top = element_text(size = 12, color = "black", margin = margin(b = 10)),
      axis.title.x.bottom = element_text(size = 12, color = "black", margin = margin(t = 10)),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", size = 0.5),
      panel.grid.minor.x = element_blank(),
      axis.ticks.x.top = element_line(color = "black"),
      axis.ticks.x.bottom = element_line(color = "black"),
      axis.ticks.y = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.line.x.top = element_line(color = "black", size = 0.5),
      axis.line.x.bottom = element_line(color = "black", size = 0.5),
      axis.line.y.left = element_line(color = "black", size = 0.5),
      legend.position = "none"
    )
  
  if(show_legend) {
    legend_x <- pvalue_max * 0.7
    legend_y_top <- length(levels(data$Description))-6
    p <- p + 
      annotate("rect", 
               xmin = legend_x, xmax = legend_x + pvalue_max * 0.08,
               ymin = legend_y_top - 0.8, ymax = legend_y_top - 0.4,
               fill = bar_color, alpha = 0.7) +
      annotate("text", 
               x = legend_x + pvalue_max * 0.12, 
               y = legend_y_top - 0.6,
               label = x_lab, hjust = 0, size = 3.5) +
      
      annotate("point", 
               x = legend_x + pvalue_max * 0.04, 
               y = legend_y_top - 1.4,
               color = "black", size = 1) +
      annotate("segment", 
               x = legend_x, xend = legend_x + pvalue_max * 0.08,
               y = legend_y_top - 1.4, yend = legend_y_top - 1.4,
               color = "black", size = 0.5) +
      annotate("text", 
               x = legend_x + pvalue_max * 0.12, 
               y = legend_y_top - 1.4,
               label = count_lab, hjust = 0, size = 3.5)
  }
  
  return(p)
}



gobp_data <- process_enrichment(go_bp)
gobp_plot <- plot_enrichment_barplot(data = gobp_data,title = "GOBP enrichment",bar_color = "#7d578c")
print(gobp_plot)
gocc_data <- process_enrichment(go_cc)
gocc_plot <- plot_enrichment_barplot(data = gocc_data,title = "GOCC enrichment",bar_color = "#7d578c")
print(gocc_plot)

gomf_data <- process_enrichment(go_mf)
gomf_plot <- plot_enrichment_barplot(data = gomf_data,title = "GOMF enrichment",bar_color = "#7d578c")
print(gomf_plot)
kegg_data <- process_enrichment(kegg_result)
kegg_plot <- plot_enrichment_barplot(data = kegg_data,title = "KEGG enrichment",bar_color = "#4e62a5")
print(kegg_plot)
ggsave("Fig3s_Recurrent_TE_Loci_recurrent_up_intron_intergenic_TE_Loci_closest_genes_GOBP_enrichment.pdf", gobp_plot, width = 6, height = 4)
ggsave("Fig3s_Recurrent_TE_Loci_recurrent_up_intron_intergenic_TE_Loci_closest_genes_GOCC_enrichment.pdf", gocc_plot, width = 8, height = 2)
ggsave("Fig3s_Recurrent_TE_Loci_recurrent_up_intron_intergenic_TE_Loci_closest_genes_GOMF_enrichment.pdf", gomf_plot, width = 8, height = 4)
ggsave("Fig3s_Recurrent_TE_Loci_recurrent_up_intron_intergenic_TE_Loci_closest_genes_KEGG_enrichment.pdf", kegg_plot, width = 7, height = 4)
