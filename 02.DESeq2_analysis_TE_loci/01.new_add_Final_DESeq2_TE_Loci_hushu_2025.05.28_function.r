library(DESeq2)
library(tidyverse)
library(dplyr)
library(openxlsx)

source("~/pathogen_TE_2025/02.DESeq2_analysis/00.Functions.r")
setwd("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_loci/")

TE_DE_analysis(
  metadata_file  = "/data2t_2/hushu/02.DESeq2_analysis/raw_data/RNA_seq_metadata_new_add_250525.xlsx",
  raw_count_file = "/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data/Virus_RSV_3_readscounts_matrix_TE_Loci.txt",
  outDir         = "./new_add_DE_results/",
  dataset        = "Virus_RSV_3",
  cond1          = "Pre_0m",
  cond2          = "Post_960m"
)

TE_DE_analysis(
  metadata_file  = "/data2t_2/hushu/02.DESeq2_analysis/raw_data/RNA_seq_metadata_new_add_250525.xlsx",
  raw_count_file = "/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data/Virus_RSV_3_readscounts_matrix_TE_Loci.txt",
  outDir         = "./new_add_DE_results/",
  dataset        = "Virus_RSV_3",
  cond1          = "Pre_0m",
  cond2          = "Post_1440m"
)


##Virus_SARSCoV2_1
TE_DE_analysis(
  metadata_file  = "/data2t_2/hushu/02.DESeq2_analysis/raw_data/RNA_seq_metadata_new_add_250525.xlsx",
  raw_count_file = "/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data/Virus_SARSCoV2_1_readscounts_matrix_TE_Loci.txt",
  outDir         = "./new_add_DE_results/",
  dataset        = "Virus_SARSCoV2_1",
  cond1          = "Pre_120m",
  cond2          = "Post_120m"
)
TE_DE_analysis(
  metadata_file  = "/data2t_2/hushu/02.DESeq2_analysis/raw_data/RNA_seq_metadata_new_add_250525.xlsx",
  raw_count_file = "/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data/Virus_SARSCoV2_1_readscounts_matrix_TE_Loci.txt",
  outDir         = "./new_add_DE_results/",
  dataset        = "Virus_SARSCoV2_1",
  cond1          = "Pre_240m",
  cond2          = "Post_240m"
)

TE_DE_analysis(
  metadata_file  = "/data2t_2/hushu/02.DESeq2_analysis/raw_data/RNA_seq_metadata_new_add_250525.xlsx",
  raw_count_file = "/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data/Virus_SARSCoV2_1_readscounts_matrix_TE_Loci.txt",
  outDir         = "./new_add_DE_results/",
  dataset        = "Virus_SARSCoV2_1",
  cond1          = "Pre_360m",
  cond2          = "Post_360m"
)
TE_DE_analysis(
  metadata_file  = "/data2t_2/hushu/02.DESeq2_analysis/raw_data/RNA_seq_metadata_new_add_250525.xlsx",
  raw_count_file = "/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data/Virus_SARSCoV2_1_readscounts_matrix_TE_Loci.txt",
  outDir         = "./new_add_DE_results/",
  dataset        = "Virus_SARSCoV2_1",
  cond1          = "Pre_540m",
  cond2          = "Post_540m"
)
TE_DE_analysis(
  metadata_file  = "/data2t_2/hushu/02.DESeq2_analysis/raw_data/RNA_seq_metadata_new_add_250525.xlsx",
  raw_count_file = "/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data/Virus_SARSCoV2_1_readscounts_matrix_TE_Loci.txt",
  outDir         = "./new_add_DE_results/",
  dataset        = "Virus_SARSCoV2_1",
  cond1          = "Pre_1120m",
  cond2          = "Post_1120m"
)
TE_DE_analysis(
  metadata_file  = "/data2t_2/hushu/02.DESeq2_analysis/raw_data/RNA_seq_metadata_new_add_250525.xlsx",
  raw_count_file = "/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data/Virus_SARSCoV2_1_readscounts_matrix_TE_Loci.txt",
  outDir         = "./new_add_DE_results/",
  dataset        = "Virus_SARSCoV2_1",
  cond1          = "Pre_2240m",
  cond2          = "Post_2240m"
)

##Virus_SARSCoV2_2
TE_DE_analysis(
  metadata_file  = "/data2t_2/hushu/02.DESeq2_analysis/raw_data/RNA_seq_metadata_new_add_250525.xlsx",
  raw_count_file = "/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data/Virus_SARSCoV2_2_readscounts_matrix_TE_Loci.txt",
  outDir         = "./new_add_DE_results/",
  dataset        = "Virus_SARSCoV2_2",
  cond1          = "Pre_1440m",
  cond2          = "Post_1440m"
)
TE_DE_analysis(
  metadata_file  = "/data2t_2/hushu/02.DESeq2_analysis/raw_data/RNA_seq_metadata_new_add_250525.xlsx",
  raw_count_file = "/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data/Virus_SARSCoV2_2_readscounts_matrix_TE_Loci.txt",
  outDir         = "./new_add_DE_results/",
  dataset        = "Virus_SARSCoV2_2",
  cond1          = "Pre_2880m",
  cond2          = "Post_2880m"
)


##Bacteria_MTB_1
TE_DE_analysis(
  metadata_file  = "/data2t_2/hushu/02.DESeq2_analysis/raw_data/RNA_seq_metadata_new_add_250525.xlsx",
  raw_count_file = "/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data/Bacteria_MTB_1_readscounts_matrix_TE_Loci.txt",
  outDir         = "./new_add_DE_results/",
  dataset        = "Bacteria_MTB_1",
  cond1          = "Pre_1080m",
  cond2          = "Post_1080m"
)


##Virus_HCV_1 Huh-7
TE_DE_analysis(
  metadata_file  = "/data2t_2/hushu/02.DESeq2_analysis/raw_data/RNA_seq_metadata_new_add_250525.xlsx",
  raw_count_file = "/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data/Virus_HCV_1_readscounts_matrix_TE_Loci.txt",
  outDir         = "./new_add_DE_results/",
  dataset        = "Virus_HCV_1",
  cond1          = "Pre_4320m_Huh-7",
  cond2          = "Post_4320m_Huh-7"
)

TE_DE_analysis(
  metadata_file  = "/data2t_2/hushu/02.DESeq2_analysis/raw_data/RNA_seq_metadata_new_add_250525.xlsx",
  raw_count_file = "/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data/Virus_HCV_1_readscounts_matrix_TE_Loci.txt",
  outDir         = "./new_add_DE_results/",
  dataset        = "Virus_HCV_1",
  cond1          = "Pre_4320m_Huh-7",
  cond2          = "Post_7200m_Huh-7"
)

##Virus_HCV_1 Huh-7.5.1
TE_DE_analysis(
  metadata_file  = "/data2t_2/hushu/02.DESeq2_analysis/raw_data/RNA_seq_metadata_new_add_250525.xlsx",
  raw_count_file = "/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data/Virus_HCV_1_readscounts_matrix_TE_Loci.txt",
  outDir         = "./new_add_DE_results/",
  dataset        = "Virus_HCV_1",
  cond1          = "Pre_4320m_Huh-7.5.1",
  cond2          = "Post_4320m_Huh-7.5.1"
)

TE_DE_analysis(
  metadata_file  = "/data2t_2/hushu/02.DESeq2_analysis/raw_data/RNA_seq_metadata_new_add_250525.xlsx",
  raw_count_file = "/data2t_2/hushu/02.DESeq2_analysis_TE_loci/new_add_raw_data/Virus_HCV_1_readscounts_matrix_TE_Loci.txt",
  outDir         = "./new_add_DE_results/",
  dataset        = "Virus_HCV_1",
  cond1          = "Pre_4320m_Huh-7.5.1",
  cond2          = "Post_7200m_Huh-7.5.1"
)


###合并同一个species的差异基因
###A549 Virus_IAV
setwd("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_loci/DE_results/A549/")

union_DE_genes_TE_loci(
  prefix = "Virus_IAV",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./"
)

union_DE_genes_TE_loci(
  prefix = "Virus_SARSCoV2",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./"
)

union_DE_genes_TE_loci(
  prefix = "Virus_RSV",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./"
)

union_DE_genes_TE_loci(
  prefix = "Bacteria_Spn",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./"
)

##A549每个dataset分开得到up或down或/的结果
setwd("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_loci/DE_results/A549/")
###A549 Virus_IAV_11
union_DE_genes_TE_loci(
  prefix = "Virus_IAV_11",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./A549_split/"
)

union_DE_genes_TE_loci(
  prefix = "Virus_SARSCoV2_1",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./A549_split/"
)
union_DE_genes_TE_loci(
  prefix = "Virus_SARSCoV2_2",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./A549_split/"
)

union_DE_genes_TE_loci(
  prefix = "Virus_SARSCoV2_3",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./A549_split/"
)



###合并同一个species的差异基因
###Dendritic Bacteria_Bpe
setwd("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_loci/DE_results/Dendritic/")

union_DE_genes_TE_loci(
  prefix = "Bacteria_Bpe",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./"
)

union_DE_genes_TE_loci(
  prefix = "Bacteria_Sau",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./"
)

union_DE_genes_TE_loci(
  prefix = "Protozoa_Pf",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./"
)

union_DE_genes_TE_loci(
  prefix = "Virus_HBV",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./"
)

union_DE_genes_TE_loci(
  prefix = "Bacteria_MTB",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./"
)

setwd("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_loci/DE_results/Dendritic/")
union_DE_genes_TE_loci(
  prefix = "Bacteria_Bpe_2",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./Dendritic_split/"
)

union_DE_genes_TE_loci(
  prefix = "Bacteria_Sau_2",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./Dendritic_split/"
)
union_DE_genes_TE_loci(
  prefix = "Protozoa_Pf_6",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./Dendritic_split/"
)

union_DE_genes_TE_loci(
  prefix = "Virus_HBV_1",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./Dendritic_split/"
)

union_DE_genes_TE_loci(
  prefix = "Bacteria_MTB_1",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./Dendritic_split/"
)


setwd("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_loci/DE_results/hSAEC/")
union_DE_genes_TE_loci(
  prefix = "Virus_RSV",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./"
)


union_DE_genes_TE_loci(
  prefix = "Virus_RSV_3",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./hSAEC_split/"
)


setwd("/data2t_2/pathogen_TE_2025_New/02.DESeq2_analysis_TE_loci/DE_results/Huh-7/")
union_DE_genes_TE_loci(
  prefix = "Virus_HCV",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./"
)


union_DE_genes_TE_loci(
  prefix = "Virus_HCV",
  suffix = "sig_DE-TEs_padj0.05log2FC1.csv",
  outDir = "./Huh-7_split/"
)


