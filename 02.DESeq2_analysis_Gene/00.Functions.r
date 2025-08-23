TE_Gene_analysis <- function(
    metadata_file,   # e.g. "raw_data/new_all_sample_info_A549_n248.xlsx"
    raw_count_file,  # e.g. "raw_data/new_all_sample_readscounts_matrix_combined_n4842_TE_subfamily.txt"
    outDir,          # e.g. "DE_results/"
    dataset,         # e.g. "Virus_IAV_11"
    cond1,           # e.g. "Pre_0m"
    cond2            # e.g. "Post_360m"
) {
  #-------- deal with the problem: two pathogens in one dataset ----------------------------
  metadata_input <- read.xlsx(metadata_file)
  subset_info <- metadata_input[metadata_input$Dataset == dataset,]
  
  #--------------------- DE analysis (TE)-----------------------------------------------------------------------
  
  raw_count <- read.csv(raw_count_file, sep="\t", row.names = 1)
  outDir <- outDir
  subset_df_TE <- raw_count[, subset_info$Sample_Name]
  
  # ------------- filtered TE------------------------------
  subset_df_TE_1 <- subset_df_TE[rowSums(subset_df_TE) > 2, ]  #n=1073
  #---------get sample info form count table ---------------------------
  count <- subset_df_TE_1  #n=1066
  Condition <- subset_info$Condition
  
  sample <- data.frame(Sample = subset_info$Sample_Name,
                       Condition = Condition)
  
  sample <- sample[sample$Condition %in% c(cond1, cond2),]
  count <- count[, sample$Sample]
  ######################################
  colnames(sample)[2] <- "Condition"
  rownames(sample) <- NULL
  sample <- column_to_rownames(sample, var = "Sample")
  sample$Condition <- factor(sample$Condition)
  sample <- sample %>%
    mutate(order = match(rownames(sample), colnames(count))) %>%
    arrange(order) %>%
    select(-order) 
  
  print(head(count))
  print(head(sample))
  
  #------------------------DESeq2-------------------------------------------
  dataset_dds1 <- DESeqDataSetFromMatrix(countData = count, colData = sample, design = ~ Condition)
  summary(dataset_dds1)
  
  #size factors
  dataset_dds2 <- estimateSizeFactors(dataset_dds1)
  summary(dataset_dds2)
  dataset_dds3 <- DESeq(dataset_dds2, quiet = F)
  summary(dataset_dds3)
  
  dataset_count_norm <- counts(dataset_dds3, normalized = T)
  summary(dataset_count_norm)
  dataset_count_norm1 <- as.data.frame(dataset_count_norm)
  
  
  dataset_res <- results(dataset_dds3, contrast = c("Condition", cond2, cond1))
  dataset_res_ord <- dataset_res[order(dataset_res$padj),]
  print(head(dataset_res_ord)) 
  summary(dataset_res_ord) 
  
  # plotCounts(dataset_dds_norm, gene="HERVI-int:ERV1:LTR", intgroup="Condition") #
  # plotCounts(dataset_dds_norm, gene="AluSp_c5104420", intgroup="Condition") #
  
  res_temp <- as.data.frame(dataset_res_ord) #n=
  res_final <- na.omit(res_temp)  
  
  res_final_ord <- res_final[order(res_final$padj),]
  write.csv(res_final_ord, paste0(outDir, "/", dataset, "_", cond2, "_vs_", cond1, "_all_res_Gene.csv"), na = "NA")
  
  
  res_sig <- res_final_ord[which(res_final_ord$padj < 0.05 & abs(res_final_ord$log2FoldChange) > 0.5),]  #n=xxxxx
  
  if (nrow(res_sig) == 0) {
    print("This dataset has 0 significant DE Gene (padj<0.05, log2FC>0.5) !")
  } else {
    
    res_sig[which(res_sig$log2FoldChange > 0), 'sig'] <- 'Up'
    res_sig[which(res_sig$log2FoldChange < 0), 'sig'] <- 'Down'
    
    res_sig2 <- as.data.frame(res_sig)
    print(head(res_sig2))
    
    res_sig3 <- res_sig2[order(res_sig2$padj),]
    
    # res_TE <- filter(res_sig3, str_detect(rownames(res_sig3), ":"))
    res_TE<-res_sig3
    print("padj<0.05, log2FC>0.5条件下DE TE/TE_Loci/Gene数量:")
    print(nrow(res_sig3))
    write.csv(res_TE, paste0(outDir, "/", dataset, "_", cond2, "_vs_", cond1, "_sig_DE-Gene_padj0.05log2FC0.5.csv"),  na = "NA", row.names = T)
  }
    res_sig <- res_final_ord[which(res_final_ord$padj < 0.05 & abs(res_final_ord$log2FoldChange) > 1),]  #n=xxxxx
  
  if (nrow(res_sig) == 0) {
    print("This dataset has 0 significant DE TE  (padj<0.05, log2FC>1) !")
  } else {
    
    res_sig[which(res_sig$log2FoldChange > 0), 'sig'] <- 'Up'
    res_sig[which(res_sig$log2FoldChange < 0), 'sig'] <- 'Down'
    
    res_sig2 <- as.data.frame(res_sig)
    print(head(res_sig2))
    
    res_sig3 <- res_sig2[order(res_sig2$padj),]
    
    # res_Gene <- filter(res_sig3, !str_detect(rownames(res_sig3), ":"))
    # res_TE <- filter(res_sig3, str_detect(rownames(res_sig3), ":"))
    res_TE<-res_sig3
    print("padj<0.05, log2FC>1条件下DE TE/TE_Loci/Gene数量:")
    print(nrow(res_sig3))
    write.csv(res_TE, paste0(outDir, "/", dataset, "_", cond2, "_vs_", cond1, "_sig_DE-Gene_padj0.05log2FC1.csv"),  na = "NA", row.names = T)
  }
}

union_DE_genes <- function(
    prefix,  #eg."Virus_IAV",           
    suffix,  #eg."sig_DE-Gene_padj0.05log2FC0.5.csv" 
    outDir   #eg."DE_results/A549" 
) {
  raw_count <- read.csv("~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",sep = "\t")
  all_genes <- raw_count$GeneID
  
  result_df <- data.frame(GeneID = all_genes, stringsAsFactors = FALSE)
  result_df[[prefix]] <- "/"

  files <- list.files(pattern = paste0("^", prefix, ".*_", suffix))
  
  all_de_genes <- data.frame(GeneID = character(),
                             Regulation = character(),
                             stringsAsFactors = FALSE)
  
  for (file in files) {
    de_data <- read.csv(file, row.names = 1)
    
    genes <- rownames(de_data)
    regulation <- de_data$sig
    
    de_genes <- data.frame(
      GeneID = genes,
      Regulation = regulation,
      stringsAsFactors = FALSE
    )
    
    all_de_genes <- rbind(all_de_genes, de_genes)
  }
  
  de_summary <- all_de_genes %>%
    dplyr::group_by(GeneID) %>%
    dplyr::summarize(Regulation = dplyr::case_when(
      all(Regulation == "Up") ~ "Up",
      all(Regulation == "Down") ~ "Down",
      TRUE ~ "Mixed"  
    ))
  
  for (i in 1:nrow(de_summary)) {
    gene <- de_summary$GeneID[i]
    reg <- de_summary$Regulation[i]
    
    idx <- which(result_df$GeneID == gene)
    if (length(idx) > 0) {
      result_df[[prefix]][idx] <- reg
    }
  }
  
  out_file_name <- paste0("All_", prefix, "_", suffix)
  out_file <- file.path(outDir, out_file_name)
  write.csv(result_df, out_file, row.names = FALSE, quote = FALSE)
}
