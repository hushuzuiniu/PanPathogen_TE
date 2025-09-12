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
setwd("/data2t_2/pathogen_TE_2025_New/03.Enrichment_analysis/")

recurrent_genes<-read.csv("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_Gene/plots/All_species_split_DE-Gene_padj0.05log2FC1.csv",header = T)
recurrent_genes<- recurrent_genes[,c("GeneID","UniqPathogens")]
recurrent_genes<-recurrent_genes[recurrent_genes$UniqPathogens>=5,]
gene_list <- unique(recurrent_genes$GeneID)
## 2893
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

go_bp_df <- save_enrichment_results(go_bp, "Recurrent_Gene_gt5_n2893/Recurrent_Gene_gt5_n2893_GOBP_enrichment")
go_mf_df <- save_enrichment_results(go_mf, "Recurrent_Gene_gt5_n2893/Recurrent_Gene_gt5_n2893_GOMF_enrichment")
go_cc_df <- save_enrichment_results(go_cc, "Recurrent_Gene_gt5_n2893/Recurrent_Gene_gt5_n2893_GOCC_enrichment")
kegg_df <- save_enrichment_results(kegg_result, "Recurrent_Gene_gt5_n2893/Recurrent_Gene_gt5_n2893_KEGG_enrichment")

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
    pvalue_breaks <- seq(0, ceiling(max_neg_log10) + 4, by = 5)
  } else if (max_neg_log10 <= 50) {
    pvalue_breaks <- seq(0, ceiling(max_neg_log10) + 5, by = 5)
  } else {
    pvalue_breaks <- seq(0, ceiling(max_neg_log10) + 10, by = 20)
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


ggsave("Recurrent_Gene_gt5_n2893/Recurrent_Gene_gt5_n2893_GOBP_enrichment.pdf", gobp_plot, width = 6, height = 4)
ggsave("Recurrent_Gene_gt5_n2893/Recurrent_Gene_gt5_n2893_GOCC_enrichment.pdf", gocc_plot, width = 6, height = 4)
ggsave("Recurrent_Gene_gt5_n2893/Recurrent_Gene_gt5_n2893_GOMF_enrichment.pdf", gomf_plot, width = 8, height = 4)
ggsave("Recurrent_Gene_gt5_n2893/Recurrent_Gene_gt5_n2893_KEGG_enrichment.pdf", kegg_plot, width = 7, height = 4)
