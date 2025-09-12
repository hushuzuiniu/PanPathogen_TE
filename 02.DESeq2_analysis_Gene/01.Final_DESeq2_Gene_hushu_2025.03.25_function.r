library(DESeq2)
library(tidyverse)
library(dplyr)
library(openxlsx)
source("~/pathogen_TE_2025/02.DESeq2_analysis_Gene/00.Functions.r")
##################################A549##############################################
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene")
##A549 celltype Virus_IAV_11 dataset Post_360m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Virus_IAV_11",
  cond1          = "Pre_0m",
  cond2          = "Post_360m"
)

###A549 celltype Virus_IAV_11 dataset Post_720m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Virus_IAV_11",
  cond1          = "Pre_0m",
  cond2          = "Post_720m"
)

###A549 celltype Virus_IAV_3 dataset Post_360m vs Pre_360m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Virus_IAV_3",
  cond1          = "Pre_360m",
  cond2          = "Post_360m"
)

###A549 celltype Virus_RSV_1 dataset Post_1440m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Virus_RSV_1",
  cond1          = "Pre_0m",
  cond2          = "Post_1440m"
)

###A549 celltype Virus_RSV_1 dataset Post_2880m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Virus_RSV_1",
  cond1          = "Pre_0m",
  cond2          = "Post_2880m"
)

###A549 celltype Virus_IAV_6 dataset Post_240m vs Pre_240m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Virus_IAV_6",
  cond1          = "Pre_240m",
  cond2          = "Post_240m"
)

###A549 celltype Virus_IAV_6 dataset Post_480m vs Pre_480m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Virus_IAV_6",
  cond1          = "Pre_480m",
  cond2          = "Post_480m"
)

###A549 celltype Virus_IAV_6 dataset Post_720m vs Pre_720m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Virus_IAV_6",
  cond1          = "Pre_720m",
  cond2          = "Post_720m"
)
###A549 celltype Virus_SARSCoV2_10 dataset Post_120m_VIC vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Virus_SARSCoV2_10",
  cond1          = "Pre_0m",
  cond2          = "Post_120m_VIC"
)

###A549 celltype Virus_SARSCoV2_10 dataset Post_480m_VIC vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Virus_SARSCoV2_10",
  cond1          = "Pre_0m",
  cond2          = "Post_480m_VIC"
)

###A549 celltype Virus_SARSCoV2_10 dataset Post_1440m_VIC vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Virus_SARSCoV2_10",
  cond1          = "Pre_0m",
  cond2          = "Post_1440m_VIC"
)

###A549 celltype Virus_SARSCoV2_10 dataset Post_120m_Alpha vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Virus_SARSCoV2_10",
  cond1          = "Pre_0m",
  cond2          = "Post_120m_Alpha"
)

###A549 celltype Virus_SARSCoV2_10 dataset Post_480m_Alpha vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Virus_SARSCoV2_10",
  cond1          = "Pre_0m",
  cond2          = "Post_480m_Alpha"
)


###A549 celltype Virus_SARSCoV2_10 dataset Post_1440m_Alpha vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Virus_SARSCoV2_10",
  cond1          = "Pre_0m",
  cond2          = "Post_1440m_Alpha"
)


###A549 celltype Virus_SARSCoV2_9 dataset Post_1440m vs Pre_1440m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Virus_SARSCoV2_9",
  cond1          = "Pre_1440m",
  cond2          = "Post_1440m"
)

###A549 celltype Virus_SARSCoV2_9 dataset Post_2880m vs Pre_2880m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Virus_SARSCoV2_9",
  cond1          = "Pre_2880m",
  cond2          = "Post_2880m"
)

###A549 celltype Virus_IAV_9 dataset Post_1440m vs Pre_1440m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Virus_IAV_9",
  cond1          = "Pre_1440m",
  cond2          = "Post_1440m"
)

###A549 celltype Virus_IAV_2 dataset Post_60m vs Pre_60m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Virus_IAV_2",
  cond1          = "Pre_60m",
  cond2          = "Post_60m"
)

###A549 celltype Bacteria_Spn_2 dataset Post_30m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Bacteria_Spn_2",
  cond1          = "Pre_0m",
  cond2          = "Post_30m"
)

###A549 celltype Bacteria_Spn_2 dataset Post_60m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Bacteria_Spn_2",
  cond1          = "Pre_0m",
  cond2          = "Post_60m"
)

###A549 celltype Bacteria_Spn_2 dataset Post_120m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_A549_n248.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/A549",
  dataset        = "Bacteria_Spn_2",
  cond1          = "Pre_0m",
  cond2          = "Post_120m"
)

#############################################################################
###A549 Virus_IAV
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene/DE_results/A549/")
union_DE_genes(
  prefix = "Virus_IAV",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)

union_DE_genes(
  prefix = "Virus_IAV",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)

###A549 Bacteria_Spn
union_DE_genes(
  prefix = "Bacteria_Spn",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)

union_DE_genes(
  prefix = "Bacteria_Spn",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)

###A549 Virus_RSV
union_DE_genes(
  prefix = "Virus_RSV",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)

union_DE_genes(
  prefix = "Virus_RSV",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)


###A549 Virus_SARSCoV2
union_DE_genes(
  prefix = "Virus_SARSCoV2",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)

union_DE_genes(
  prefix = "Virus_SARSCoV2",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)

########################################################################
##A549每个dataset分开得到up或down或/的结果
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene/DE_results/A549/")
###A549 Virus_IAV_11
union_DE_genes(
  prefix = "Virus_IAV_11",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./A549_split/"
)

union_DE_genes(
  prefix = "Virus_IAV_11",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./A549_split/"
)

###A549 Virus_IAV_3
union_DE_genes(
  prefix = "Virus_IAV_3",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./A549_split/"
)

union_DE_genes(
  prefix = "Virus_IAV_3",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./A549_split/"
)

###A549 Virus_IAV_6
union_DE_genes(
  prefix = "Virus_IAV_6",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./A549_split/"
)

union_DE_genes(
  prefix = "Virus_IAV_6",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./A549_split/"
)

###A549 Virus_IAV_2
union_DE_genes(
  prefix = "Virus_IAV_2",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./A549_split/"
)

union_DE_genes(
  prefix = "Virus_IAV_2",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./A549_split/"
)

###A549 Virus_IAV_9
union_DE_genes(
  prefix = "Virus_IAV_9",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./A549_split/"
)

union_DE_genes(
  prefix = "Virus_IAV_9",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./A549_split/"
)

###A549 Bacteria_Spn_2
union_DE_genes(
  prefix = "Bacteria_Spn_2",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./A549_split/"
)

union_DE_genes(
  prefix = "Bacteria_Spn_2",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./A549_split/"
)

###A549 Virus_RSV_1
union_DE_genes(
  prefix = "Virus_RSV_1",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./A549_split/"
)

union_DE_genes(
  prefix = "Virus_RSV_1",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./A549_split/"
)


###A549 Virus_SARSCoV2_9
union_DE_genes(
  prefix = "Virus_SARSCoV2_9",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./A549_split/"
)

union_DE_genes(
  prefix = "Virus_SARSCoV2_9",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./A549_split/"
)

###A549 Virus_SARSCoV2_10
union_DE_genes(
  prefix = "Virus_SARSCoV2_10",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./A549_split/"
)

union_DE_genes(
  prefix = "Virus_SARSCoV2_10",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./A549_split/"
)
##################################Macrophages##############################################
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene")
##Macrophages celltype Bacteria_Lmo_1 dataset Post_120m vs Pre_120m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Bacteria_Lmo_1",
  cond1          = "Pre_120m",
  cond2          = "Post_120m"
)
##Macrophages celltype Bacteria_Lmo_4 dataset Post_300m vs Pre_300m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Bacteria_Lmo_4",
  cond1          = "Pre_300m",
  cond2          = "Post_300m"
)
##Macrophages celltype Bacteria_MTB_4 dataset Post_2880m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Bacteria_MTB_4",
  cond1          = "Pre_0m",
  cond2          = "Post_2880m"
)
##Macrophages celltype Bacteria_MTB_6 dataset Post_120m vs Pre_120m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Bacteria_MTB_6",
  cond1          = "Pre_120m",
  cond2          = "Post_120m"
)
##Macrophages celltype Bacteria_MTB_6 dataset Post_2880m vs Pre_2880m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Bacteria_MTB_6",
  cond1          = "Pre_2880m",
  cond2          = "Post_2880m"
)

##Macrophages celltype Bacteria_MTB_8 dataset Post_1440m vs Pre_1440m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Bacteria_MTB_8",
  cond1          = "Pre_1440m",
  cond2          = "Post_1440m"
)

##Macrophages celltype Bacteria_NTHi_1 dataset Post_360m vs Pre_360m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Bacteria_NTHi_1",
  cond1          = "Pre_360m",
  cond2          = "Post_360m"
)
##Macrophages celltype Bacteria_NTHi_1 dataset Post_1440m vs Pre_1440m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Bacteria_NTHi_1",
  cond1          = "Pre_1440m",
  cond2          = "Post_1440m"
)
##Macrophages celltype Bacteria_SalEn_1 dataset Post_240m vs Pre_240m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Bacteria_SalEn_1",
  cond1          = "Pre_240m",
  cond2          = "Post_240m"
)
##Macrophages celltype Bacteria_SalTy_4 dataset Post_120m vs Pre_120m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Bacteria_SalTy_4",
  cond1          = "Pre_120m",
  cond2          = "Post_120m"
)

##Macrophages celltype Fungi_Af_4 dataset Post_120m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Fungi_Af_4",
  cond1          = "Pre_0m",
  cond2          = "Post_120m"
)
##Macrophages celltype Fungi_Af_4 dataset Post_360m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Fungi_Af_4",
  cond1          = "Pre_0m",
  cond2          = "Post_360m"
)

##Macrophages celltype Fungi_Af_6 dataset Post_60m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Fungi_Af_6",
  cond1          = "Pre_0m",
  cond2          = "Post_60m"
)
##Macrophages celltype Fungi_Af_6 dataset Post_360m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Fungi_Af_6",
  cond1          = "Pre_0m",
  cond2          = "Post_360m"
)
##Macrophages celltype Fungi_Af_8 dataset Post_120m vs Pre_120m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Fungi_Af_8",
  cond1          = "Pre_120m",
  cond2          = "Post_120m"
)
##Macrophages celltype Fungi_Af_8 dataset Post_240m vs Pre_240m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Fungi_Af_8",
  cond1          = "Pre_240m",
  cond2          = "Post_240m"
)
##Macrophages celltype Fungi_Af_8 dataset Post_360m vs Pre_360m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Fungi_Af_8",
  cond1          = "Pre_360m",
  cond2          = "Post_360m"
)
##Macrophages celltype Fungi_Af_8 dataset Post_480m vs Pre_480m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Fungi_Af_8",
  cond1          = "Pre_480m",
  cond2          = "Post_480m"
)
##Macrophages celltype Virus_IAV_9 dataset Post_1440m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Virus_IAV_9",
  cond1          = "Pre_0m",
  cond2          = "Post_1440m"
)
##Macrophages celltype Virus_SARSCoV2_1 dataset Post_240m vs Pre_240m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Virus_SARSCoV2_1",
  cond1          = "Pre_240m",
  cond2          = "Post_240m"
)

##Macrophages celltype Virus_SARSCoV2_1 dataset Post_720m vs Pre_720m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Virus_SARSCoV2_1",
  cond1          = "Pre_720m",
  cond2          = "Post_720m"
)
##Macrophages celltype Virus_SARSCoV2_1 dataset Post_2160m vs Pre_2160m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Macrophages_n986.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Macrophages/",
  dataset        = "Virus_SARSCoV2_1",
  cond1          = "Pre_2160m",
  cond2          = "Post_2160m"
)

##################################################################
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene/DE_results/Macrophages/")
###Macrophages Bacteria_Lmo
union_DE_genes(
  prefix = "Bacteria_Lmo",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)

union_DE_genes(
  prefix = "Bacteria_Lmo",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)

###Macrophages Bacteria_MTB
union_DE_genes(
  prefix = "Bacteria_MTB",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)

union_DE_genes(
  prefix = "Bacteria_MTB",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)

###Macrophages Bacteria_NTHi
union_DE_genes(
  prefix = "Bacteria_NTHi",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)

union_DE_genes(
  prefix = "Bacteria_NTHi",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)

###Macrophages Bacteria_SalEn
union_DE_genes(
  prefix = "Bacteria_SalEn",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)

union_DE_genes(
  prefix = "Bacteria_SalEn",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)

###Macrophages Bacteria_SalTy
union_DE_genes(
  prefix = "Bacteria_SalTy",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)

union_DE_genes(
  prefix = "Bacteria_SalTy",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)

###Macrophages Fungi_Af
union_DE_genes(
  prefix = "Fungi_Af",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)

union_DE_genes(
  prefix = "Fungi_Af",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)

###Macrophages Virus_IAV
union_DE_genes(
  prefix = "Virus_IAV",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)

union_DE_genes(
  prefix = "Virus_IAV",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)
###Macrophages Virus_SARSCoV2
union_DE_genes(
  prefix = "Virus_SARSCoV2",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)

union_DE_genes(
  prefix = "Virus_SARSCoV2",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)


####Macrophages split
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene/DE_results/Macrophages/")
###Macrophages Bacteria_Lmo_1
union_DE_genes(
  prefix = "Bacteria_Lmo_1",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Macrophages_split/"
)

union_DE_genes(
  prefix = "Bacteria_Lmo_1",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Macrophages_split/"
)
###Macrophages Bacteria_Lmo_4
union_DE_genes(
  prefix = "Bacteria_Lmo_4",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Macrophages_split/"
)

union_DE_genes(
  prefix = "Bacteria_Lmo_4",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Macrophages_split/"
)

###Macrophages Bacteria_MTB_4
union_DE_genes(
  prefix = "Bacteria_MTB_4",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Macrophages_split/"
)

union_DE_genes(
  prefix = "Bacteria_MTB_4",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Macrophages_split/"
)

###Macrophages Bacteria_MTB_6
union_DE_genes(
  prefix = "Bacteria_MTB_6",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Macrophages_split/"
)

union_DE_genes(
  prefix = "Bacteria_MTB_6",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Macrophages_split/"
)
###Macrophages Bacteria_MTB_8
union_DE_genes(
  prefix = "Bacteria_MTB_8",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Macrophages_split/"
)

union_DE_genes(
  prefix = "Bacteria_MTB_8",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Macrophages_split/"
)
###Macrophages Bacteria_NTHi_1
union_DE_genes(
  prefix = "Bacteria_NTHi_1",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Macrophages_split/"
)

union_DE_genes(
  prefix = "Bacteria_NTHi_1",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Macrophages_split/"
)

###Macrophages Bacteria_SalEn_1
union_DE_genes(
  prefix = "Bacteria_SalEn_1",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Macrophages_split/"
)

union_DE_genes(
  prefix = "Bacteria_SalEn_1",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Macrophages_split/"
)

###Macrophages Bacteria_SalTy_4
union_DE_genes(
  prefix = "Bacteria_SalTy_4",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Macrophages_split/"
)

union_DE_genes(
  prefix = "Bacteria_SalTy_4",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Macrophages_split/"
)

###Macrophages Fungi_Af_4
union_DE_genes(
  prefix = "Fungi_Af_4",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Macrophages_split/"
)

union_DE_genes(
  prefix = "Fungi_Af_4",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Macrophages_split/"
)

###Macrophages Fungi_Af_6
union_DE_genes(
  prefix = "Fungi_Af_6",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Macrophages_split/"
)

union_DE_genes(
  prefix = "Fungi_Af_6",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Macrophages_split/"
)
###Macrophages Fungi_Af_8
union_DE_genes(
  prefix = "Fungi_Af_8",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Macrophages_split/"
)

union_DE_genes(
  prefix = "Fungi_Af_8",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Macrophages_split/"
)

###Macrophages Virus_IAV_9
union_DE_genes(
  prefix = "Virus_IAV_9",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Macrophages_split/"
)

union_DE_genes(
  prefix = "Virus_IAV_9",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Macrophages_split/"
)

###Macrophages Virus_SARSCoV2_1
union_DE_genes(
  prefix = "Virus_SARSCoV2_1",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Macrophages_split/"
)

union_DE_genes(
  prefix = "Virus_SARSCoV2_1",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Macrophages_split/"
)



##################################B_cell##############################################
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene")
##B_cell celltype Virus_EBV_2 dataset Post_120m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_B_cell_n85.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/B_cell/",
  dataset        = "Virus_EBV_2",
  cond1          = "Pre_0m",
  cond2          = "Post_2880m"
)
##B_cell celltype Virus_EBV_2 dataset Post_5760m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_B_cell_n85.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/B_cell/",
  dataset        = "Virus_EBV_2",
  cond1          = "Pre_0m",
  cond2          = "Post_5760m"
)
##B_cell celltype Virus_EBV_2 dataset Post_10080m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_B_cell_n85.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/B_cell/",
  dataset        = "Virus_EBV_2",
  cond1          = "Pre_0m",
  cond2          = "Post_10080m"
)

##B_cell celltype Virus_EBV_2 dataset Post_20160m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_B_cell_n85.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/B_cell/",
  dataset        = "Virus_EBV_2",
  cond1          = "Pre_0m",
  cond2          = "Post_20160m"
)

##B_cell celltype Virus_EBV_2 dataset Post_30240m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_B_cell_n85.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/B_cell/",
  dataset        = "Virus_EBV_2",
  cond1          = "Pre_0m",
  cond2          = "Post_30240m"
)
##B_cell celltype Virus_EBV_2 dataset Post_40320m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_B_cell_n85.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/B_cell/",
  dataset        = "Virus_EBV_2",
  cond1          = "Pre_0m",
  cond2          = "Post_40320m"
)
##B_cell celltype Virus_EBV_4 dataset Post_1440m vs Pre_1440m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_B_cell_n85.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/B_cell/",
  dataset        = "Virus_EBV_4",
  cond1          = "Pre_1440m",
  cond2          = "Post_1440m"
)

##B_cell celltype Virus_EBV_5 dataset Post_2880m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_B_cell_n85.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/B_cell/",
  dataset        = "Virus_EBV_5",
  cond1          = "Pre_0m",
  cond2          = "Post_2880m"
)

##B_cell celltype Virus_EBV_5 dataset Post_10080m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_B_cell_n85.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/B_cell/",
  dataset        = "Virus_EBV_5",
  cond1          = "Pre_0m",
  cond2          = "Post_10080m"
)

##B_cell celltype Virus_EBV_5 dataset Post_30240m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_B_cell_n85.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/B_cell/",
  dataset        = "Virus_EBV_5",
  cond1          = "Pre_0m",
  cond2          = "Post_30240m"
)

##B_cell celltype Protozoa_Pf_5 dataset Post_2880m vs Pre_2880m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_B_cell_n85.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/B_cell/",
  dataset        = "Protozoa_Pf_5",
  cond1          = "Pre_2880m",
  cond2          = "Post_2880m"
)


##B_cell celltype Virus_HBV_1 dataset Post_120m vs Pre_120m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_B_cell_n85.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/B_cell/",
  dataset        = "Virus_HBV_1",
  cond1          = "Pre_120m",
  cond2          = "Post_120m"
)

#########################
###B_cell Virus_HBV
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene/DE_results/B_cell/")
union_DE_genes(
  prefix = "Virus_HBV",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)
union_DE_genes(
  prefix = "Virus_HBV",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)

###B_cell Virus_EBV
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene/DE_results/B_cell/")
union_DE_genes(
  prefix = "Virus_EBV",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)
union_DE_genes(
  prefix = "Virus_EBV",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)

###B_cell Protozoa_Pf
union_DE_genes(
  prefix = "Protozoa_Pf",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)
union_DE_genes(
  prefix = "Protozoa_Pf",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)


#########################
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene/DE_results/B_cell/")
###B_cell split Virus_HBV_1
union_DE_genes(
  prefix = "Virus_HBV_1",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./B_cell_split/"
)
union_DE_genes(
  prefix = "Virus_HBV_1",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./B_cell_split/"
)

###B_cell split Virus_EBV_4
union_DE_genes(
  prefix = "Virus_EBV_4",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./B_cell_split/"
)
union_DE_genes(
  prefix = "Virus_EBV_4",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./B_cell_split/"
)

###B_cell split Virus_EBV_5
union_DE_genes(
  prefix = "Virus_EBV_5",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./B_cell_split/"
)
union_DE_genes(
  prefix = "Virus_EBV_5",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./B_cell_split/"
)
###B_cell split Virus_EBV_2
union_DE_genes(
  prefix = "Virus_EBV_2",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./B_cell_split/"
)
union_DE_genes(
  prefix = "Virus_EBV_2",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./B_cell_split/"
)
###B_cell split Protozoa_Pf_5
union_DE_genes(
  prefix = "Protozoa_Pf_5",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./B_cell_split/"
)
union_DE_genes(
  prefix = "Protozoa_Pf_5",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./B_cell_split/"
)


##################################T_cell##############################################
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene")
##T_cell celltype Virus_HBV_1 dataset post_120m vs pre_120m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_T_cell_n108.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/T_cell/",
  dataset        = "Virus_HBV_1",
  cond1          = "pre_120m",
  cond2          = "post_120m"
)

##T_cell celltype Virus_HIV1_10 dataset Post_1440m vs Pre_1440m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_T_cell_n108.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/T_cell/",
  dataset        = "Virus_HIV1_10",
  cond1          = "Pre_1440m",
  cond2          = "Post_1440m"
)

##T_cell celltype Virus_HIV1_11 dataset Post_1440m vs Pre_1440m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_T_cell_n108.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/T_cell/",
  dataset        = "Virus_HIV1_11",
  cond1          = "Pre_1440m",
  cond2          = "Post_1440m"
)

##T_cell celltype Virus_HIV1_14 dataset Post_720m vs Pre_720m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_T_cell_n108.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/T_cell/",
  dataset        = "Virus_HIV1_14",
  cond1          = "Pre_720m",
  cond2          = "Post_720m"
)

##T_cell celltype Virus_HIV1_14 dataset Post_1440m vs Pre_1440m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_T_cell_n108.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/T_cell/",
  dataset        = "Virus_HIV1_14",
  cond1          = "Pre_1440m",
  cond2          = "Post_1440m"
)


#########################
###T_cell Virus_HBV
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene/DE_results/T_cell/")
union_DE_genes(
  prefix = "Virus_HBV",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)
union_DE_genes(
  prefix = "Virus_HBV",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)

union_DE_genes(
  prefix = "Virus_HIV1",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)
union_DE_genes(
  prefix = "Virus_HIV1",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)

#########################
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene/DE_results/T_cell/")
###T_cell split Virus_HBV_1
union_DE_genes(
  prefix = "Virus_HBV_1",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./T_cell_split/"
)
union_DE_genes(
  prefix = "Virus_HBV_1",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./T_cell_split/"
)


###T_cell split Virus_HIV1_10
union_DE_genes(
  prefix = "Virus_HIV1_10",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./T_cell_split/"
)
union_DE_genes(
  prefix = "Virus_HIV1_10",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./T_cell_split/"
)

###T_cell split Virus_HIV1_11
union_DE_genes(
  prefix = "Virus_HIV1_11",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./T_cell_split/"
)
union_DE_genes(
  prefix = "Virus_HIV1_11",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./T_cell_split/"
)
###T_cell split Virus_HIV1_14
union_DE_genes(
  prefix = "Virus_HIV1_14",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./T_cell_split/"
)
union_DE_genes(
  prefix = "Virus_HIV1_14",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./T_cell_split/"
)


##################################Liver##############################################
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene")
##Liver celltype Virus_HBV_2 dataset Post_21600m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Liver_n176.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Liver/",
  dataset        = "Virus_HBV_2",
  cond1          = "Pre_0m",
  cond2          = "Post_21600m"
)

##Liver celltype Virus_HBV_5 dataset Post_2880m vs Pre_2880m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Liver_n176.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Liver/",
  dataset        = "Virus_HBV_5",
  cond1          = "Pre_2880m",
  cond2          = "Post_2880m"
)

##Liver celltype Virus_HBV_5 dataset Post_10080m vs Pre_10080m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Liver_n176.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Liver/",
  dataset        = "Virus_HBV_5",
  cond1          = "Pre_10080m",
  cond2          = "Post_10080m"
)

##Liver celltype Virus_HBV_5 dataset Post_40320m vs Pre_40320m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Liver_n176.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Liver/",
  dataset        = "Virus_HBV_5",
  cond1          = "Pre_40320m",
  cond2          = "Post_40320m"
)

#########################
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene/DE_results/Liver/")
###Liver Virus_HBV
union_DE_genes(
  prefix = "Virus_HBV",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)
union_DE_genes(
  prefix = "Virus_HBV",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)

#########################
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene/DE_results/Liver/")
###Liver split Virus_HBV_5
union_DE_genes(
  prefix = "Virus_HBV_5",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Liver_split/"
)
union_DE_genes(
  prefix = "Virus_HBV_5",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Liver_split/"
)

###Liver split Virus_HBV_2
union_DE_genes(
  prefix = "Virus_HBV_2",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Liver_split/"
)
union_DE_genes(
  prefix = "Virus_HBV_2",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Liver_split/"
)

##################################PBMCs##############################################
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene")
##PBMCs celltype Fungi_Af_9 dataset Post_96m vs Pre_96m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_PBMCs_n243.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/PBMCs/",
  dataset        = "Fungi_Af_9",
  cond1          = "Pre_96m",
  cond2          = "Post_96m"
)

##PBMCs celltype Fungi_Af_9 dataset Post_1440m vs Pre_1440m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_PBMCs_n243.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/PBMCs/",
  dataset        = "Fungi_Af_9",
  cond1          = "Pre_1440m",
  cond2          = "Post_1440m"
)

##PBMCs celltype Virus_HBV_1 dataset Post_120m vs Pre_120m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_PBMCs_n243.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/PBMCs/",
  dataset        = "Virus_HBV_1",
  cond1          = "Pre_120m",
  cond2          = "Post_120m"
)

##PBMCs celltype Protozoa_Pf_2 dataset Post_0m vs Pre_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_PBMCs_n243.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/PBMCs/",
  dataset        = "Protozoa_Pf_2",
  cond1          = "Pre_0m",
  cond2          = "Post_0m"
)

##PBMCs celltype Protozoa_Pf_2 dataset Post_129600m vs Pre_129600m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_PBMCs_n243.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/PBMCs/",
  dataset        = "Protozoa_Pf_2",
  cond1          = "Pre_129600m",
  cond2          = "Post_129600m"
)
##PBMCs celltype Protozoa_Pf_2 dataset Post_216000m vs Pre_216000m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_PBMCs_n243.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/PBMCs/",
  dataset        = "Protozoa_Pf_2",
  cond1          = "Pre_216000m",
  cond2          = "Post_216000m"
)

# ##PBMCs celltype Protozoa_Pf_2 dataset Post_216000m vs Pre_0m
# TE_Gene_analysis(
#   metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_PBMCs_n243.xlsx",
#   raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
#   outDir         = "DE_results/PBMCs/",
#   dataset        = "Protozoa_Pf_2",
#   cond1          = "Pre_0m",
#   cond2          = "Post_216000m"
# )


#########################
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene/DE_results/PBMCs/")
###PBMCs Fungi_Af
union_DE_genes(
  prefix = "Fungi_Af",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)
union_DE_genes(
  prefix = "Fungi_Af",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)

###PBMCs Protozoa_Pf
union_DE_genes(
  prefix = "Protozoa_Pf",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)
union_DE_genes(
  prefix = "Protozoa_Pf",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)
###PBMCs Virus_HBV
union_DE_genes(
  prefix = "Virus_HBV",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)
union_DE_genes(
  prefix = "Virus_HBV",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)

#########################
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene/DE_results/PBMCs/")
###PBMCs Fungi_Af_9
union_DE_genes(
  prefix = "Fungi_Af_9",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./PBMCs_split/"
)
union_DE_genes(
  prefix = "Fungi_Af_9",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./PBMCs_split/"
)

###PBMCs Protozoa_Pf_2
union_DE_genes(
  prefix = "Protozoa_Pf_2",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./PBMCs_split/"
)
union_DE_genes(
  prefix = "Protozoa_Pf_2",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./PBMCs_split/"
)

###PBMCs Virus_HBV_1
union_DE_genes(
  prefix = "Virus_HBV_1",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./PBMCs_split/"
)
union_DE_genes(
  prefix = "Virus_HBV_1",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./PBMCs_split/"
)


##################################Monocyte##############################################
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene")
##Monocyte celltype Bacteria_Lmo_6 dataset Post_Neonate_120m vs Pre_Neonate_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Monocyte_n283.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Monocyte/",
  dataset        = "Bacteria_Lmo_6",
  cond1          = "Pre_Neonate_0m",
  cond2          = "Post_Neonate_120m"
)

##Monocyte celltype Bacteria_Lmo_6 dataset Post_Neonate_360m vs Pre_Neonate_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Monocyte_n283.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Monocyte/",
  dataset        = "Bacteria_Lmo_6",
  cond1          = "Pre_Neonate_0m",
  cond2          = "Post_Neonate_360m"
)

##Monocyte celltype Bacteria_Lmo_6 dataset Post_Adult_120m vs Pre_Adult_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Monocyte_n283.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Monocyte/",
  dataset        = "Bacteria_Lmo_6",
  cond1          = "Pre_Adult_0m",
  cond2          = "Post_Adult_120m"
)

##Monocyte celltype Bacteria_Lmo_6 dataset Post_Adult_360m vs Pre_Adult_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Monocyte_n283.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Monocyte/",
  dataset        = "Bacteria_Lmo_6",
  cond1          = "Pre_Adult_0m",
  cond2          = "Post_Adult_360m"
)

##Monocyte celltype Bacteria_Lmo_6 dataset Post_OlderAdult_120m vs Pre_OlderAdult_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Monocyte_n283.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Monocyte/",
  dataset        = "Bacteria_Lmo_6",
  cond1          = "Pre_OlderAdult_0m",
  cond2          = "Post_OlderAdult_120m"
)

##Monocyte celltype Bacteria_Lmo_6 dataset Post_OlderAdult_360m vs Pre_OlderAdult_0m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Monocyte_n283.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Monocyte/",
  dataset        = "Bacteria_Lmo_6",
  cond1          = "Pre_OlderAdult_0m",
  cond2          = "Post_OlderAdult_360m"
)

##Monocyte celltype Virus_RSV_5 dataset Post_1440m_RSVA vs Pre_1440m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Monocyte_n283.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Monocyte/",
  dataset        = "Virus_RSV_5",
  cond1          = "Pre_1440m",
  cond2          = "Post_1440m_RSVA"
)

##Monocyte celltype Virus_RSV_5 dataset Post_1440m_RSVB vs Pre_1440m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Monocyte_n283.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Monocyte/",
  dataset        = "Virus_RSV_5",
  cond1          = "Pre_1440m",
  cond2          = "Post_1440m_RSVB"
)

#########################
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene/DE_results/Monocyte/")
###Monocyte Bacteria_Lmo
union_DE_genes(
  prefix = "Bacteria_Lmo",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)
union_DE_genes(
  prefix = "Bacteria_Lmo",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)
###Monocyte Virus_RSV
union_DE_genes(
  prefix = "Virus_RSV",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)
union_DE_genes(
  prefix = "Virus_RSV",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)


#########################
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene/DE_results/Monocyte/")
###Monocyte Bacteria_Lmo_6
union_DE_genes(
  prefix = "Bacteria_Lmo_6",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Monocyte_split/"
)
union_DE_genes(
  prefix = "Bacteria_Lmo_6",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Monocyte_split/"
)
###Monocyte Virus_RSV_5
union_DE_genes(
  prefix = "Virus_RSV_5",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Monocyte_split/"
)
union_DE_genes(
  prefix = "Virus_RSV_5",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Monocyte_split/"
)


##################################Dendritic##############################################
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene")
##Dendritic celltype Bacteria_Bpe_2 dataset Post_180m_strainA vs Pre_180m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Dendritic_n289.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Dendritic/",
  dataset        = "Bacteria_Bpe_2",
  cond1          = "Pre_180m",
  cond2          = "Post_180m_strainA"
)
##Dendritic celltype Bacteria_Bpe_2 dataset Post_180m_strainB vs Pre_180m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Dendritic_n289.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Dendritic/",
  dataset        = "Bacteria_Bpe_2",
  cond1          = "Pre_180m",
  cond2          = "Post_180m_strainB"
)

##Dendritic celltype Bacteria_Bpe_2 dataset Post_180m_strainD vs Pre_180m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Dendritic_n289.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Dendritic/",
  dataset        = "Bacteria_Bpe_2",
  cond1          = "Pre_180m",
  cond2          = "Post_180m_strainD"
)

##Dendritic celltype Bacteria_Bpe_2 dataset Post_180m_strainF vs Pre_180m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Dendritic_n289.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Dendritic/",
  dataset        = "Bacteria_Bpe_2",
  cond1          = "Pre_180m",
  cond2          = "Post_180m_strainF"
)

##Dendritic celltype Bacteria_Bpe_2 dataset Post_180m_strainG vs Pre_180m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Dendritic_n289.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Dendritic/",
  dataset        = "Bacteria_Bpe_2",
  cond1          = "Pre_180m",
  cond2          = "Post_180m_strainG"
)
##Dendritic celltype Bacteria_Bpe_2 dataset Post_180m_strainH vs Pre_180m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Dendritic_n289.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Dendritic/",
  dataset        = "Bacteria_Bpe_2",
  cond1          = "Pre_180m",
  cond2          = "Post_180m_strainH"
)
##Dendritic celltype Bacteria_Bpe_2 dataset Post_180m_strainJ vs Pre_180m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Dendritic_n289.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Dendritic/",
  dataset        = "Bacteria_Bpe_2",
  cond1          = "Pre_180m",
  cond2          = "Post_180m_strainJ"
)
##Dendritic celltype Bacteria_Bpe_2 dataset Post_180m_strainK vs Pre_180m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Dendritic_n289.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Dendritic/",
  dataset        = "Bacteria_Bpe_2",
  cond1          = "Pre_180m",
  cond2          = "Post_180m_strainK"
)

##Dendritic celltype Bacteria_Bpe_2 dataset Post_180m_strainL vs Pre_180m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Dendritic_n289.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Dendritic/",
  dataset        = "Bacteria_Bpe_2",
  cond1          = "Pre_180m",
  cond2          = "Post_180m_strainL"
)

##Dendritic celltype Bacteria_Sau_2 dataset Post_1440m vs Pre_1440m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Dendritic_n289.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Dendritic/",
  dataset        = "Bacteria_Sau_2",
  cond1          = "Pre_1440m",
  cond2          = "Post_1440m"
)

##Dendritic celltype Protozoa_Pf_6 dataset Post_360m vs Pre_360m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Dendritic_n289.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Dendritic/",
  dataset        = "Protozoa_Pf_6",
  cond1          = "Pre_360m",
  cond2          = "Post_360m"
)

##Dendritic celltype Virus_HBV_1 dataset Post_120m vs Pre_120m
TE_Gene_analysis(
  metadata_file  = "~/pathogen_TE_2025/02.DESeq2_analysis/raw_data/new_all_sample_info_Dendritic_n289.xlsx",
  raw_count_file = "~/pathogen_TE_2025/02.DESeq2_analysis_Gene/raw_data/new_all_sample_readscounts_matrix_combined_n4842_Gene.txt",
  outDir         = "DE_results/Dendritic/",
  dataset        = "Virus_HBV_1",
  cond1          = "Pre_120m",
  cond2          = "Post_120m"
)

#########################
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene/DE_results/Dendritic/")
###Dendritic Bacteria_Bpe
union_DE_genes(
  prefix = "Bacteria_Bpe",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)
union_DE_genes(
  prefix = "Bacteria_Bpe",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)
###Dendritic Bacteria_Sau
union_DE_genes(
  prefix = "Bacteria_Sau",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)
union_DE_genes(
  prefix = "Bacteria_Sau",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)

###Dendritic Protozoa_Pf
union_DE_genes(
  prefix = "Protozoa_Pf",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)
union_DE_genes(
  prefix = "Protozoa_Pf",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)

###Dendritic Virus_HBV
union_DE_genes(
  prefix = "Virus_HBV",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./"
)
union_DE_genes(
  prefix = "Virus_HBV",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./"
)
#########################
setwd("~/pathogen_TE_2025/02.DESeq2_analysis_Gene/DE_results/Dendritic/")
###Dendritic split Bacteria_Bpe_2
union_DE_genes(
  prefix = "Bacteria_Bpe_2",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Dendritic_split/"
)
union_DE_genes(
  prefix = "Bacteria_Bpe_2",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Dendritic_split/"
)
###Dendritic split Bacteria_Sau_2
union_DE_genes(
  prefix = "Bacteria_Sau_2",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Dendritic_split/"
)
union_DE_genes(
  prefix = "Bacteria_Sau_2",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Dendritic_split/"
)

###Dendritic split Protozoa_Pf_6
union_DE_genes(
  prefix = "Protozoa_Pf_6",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Dendritic_split/"
)
union_DE_genes(
  prefix = "Protozoa_Pf_6",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Dendritic_split/"
)
###Dendritic split Virus_HBV_1
union_DE_genes(
  prefix = "Virus_HBV_1",
  suffix = "sig_DE-Gene_padj0.05log2FC0.5.csv",
  outDir = "./Dendritic_split/"
)
union_DE_genes(
  prefix = "Virus_HBV_1",
  suffix = "sig_DE-Gene_padj0.05log2FC1.csv",
  outDir = "./Dendritic_split/"
)

