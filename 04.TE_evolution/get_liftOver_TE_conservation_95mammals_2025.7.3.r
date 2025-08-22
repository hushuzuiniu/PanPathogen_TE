library(tidyverse)
library(tidyr)
library(RColorBrewer)


setwd("C:/Users/28257/CuiLab Dropbox/WEI Xiaoman/pathogen_TE_2025_New/09.TE_evolution/01.liftover")

df <- read.table("./input1.combined_TE_liftOver_results_n95_mammals.txt") #n=11210
colnames(df) <- c("TE_Copy","TE","Type","Species")

map <- filter(df, df$Type=="Mapped") #n=5605
map2 <- subset(map, select = -c(Type))  #n=5605
map2_wide <- pivot_wider(map2, names_from = "TE", values_from = "TE_Copy") #n=95

name <- read.csv("./input2.mammals_uniq_n95_name_list_NEW.csv", fileEncoding = "UTF-8") 

df_final <- merge(map2_wide, name, by.x = "Species", by.y = "SpeciesAbrreviation") #n=64
df_final_ord <- df_final[order(df_final$AluJr, decreasing = T),]

df_final_ord_new <- cbind(df_final_ord[, c("LatinName","CommonName","ChineseName","Species","ChainFile")],
                          df_final_ord[, !(names( df_final_ord) %in% c("LatinName","CommonName","ChineseName","Species","ChainFile"))])


copy_df <- openxlsx::read.xlsx("./input3.rmsk_TE_CopyNumber_Subfamily_Family_Class_n1073_v20240110.xlsx")
human_te_copy_all <- copy_df[copy_df$TE_Subfamily %in% colnames(df_final_ord_new)[-c(1:5)],]
human_te_copy <- human_te_copy_all[,c(1:2)]

human_te_copy_wide <- pivot_wider(human_te_copy , names_from = TE_Subfamily, values_from = CopyNumber)
human_te_copy_wide_ord <- human_te_copy_wide[,colnames(df_final_ord_new)[-c(1:5)]]
prefix_info <- data.frame(
  LatinName    = "Homo_sapiens",
  CommonName   = "Human",
  ChineseName  = "人类",
  Species      = "Hg38",
  ChainFile    = "NA"
)
df_human <- cbind(prefix_info, human_te_copy_wide_ord)
liftover_all_new <- rbind(df_human,df_final_ord_new)

# save Chinese characters using write_excel_csv instead of write.csv
write_excel_csv(liftover_all_new, "./output1.TE_liftOver_results_n95_mammals_plus_human.csv") 


#############   plot orthologous heatmap ###########################

selected_species <- read.table("./input4.selected_23_mammal_species_ordered_by_tree.txt")
res_selected <- liftover_all_new[liftover_all_new$LatinName %in% c(selected_species$V1),]

res_df <- res_selected[,c(4, 6:length(res_selected))]
rownames(res_df) <- NULL
res_df2 <- column_to_rownames(res_df, "Species")
res_df2 <- as.matrix(res_df2)

new_df <- sweep(res_df2, 2, res_df2[1, ], FUN = "/")
res_selected_ord <- res_selected[match(selected_species$V1, res_selected$LatinName),]
new_df_ord <- new_df[res_selected_ord$Species, ]

dfam_old <- openxlsx::read.xlsx("./dfam_results/TE_loci_and_subfamily_n59_Species_Taxa_from_Dfam.xlsx")
our_te_n55 <- read.table("./Recurrent_TE_loci_and_subfamily_uniq_n55.txt")

our_te_n55_dfam_part <- merge(our_te_n55,dfam_old, by.x = "V1", by.y = "TE", all.x = T)

openxlsx::write.xlsx(our_te_n55_dfam_part, "./dfam_results/our_te_n55_dfam_part_results.xlsx")

taxa <- read.csv("./TE_n55_Species_Taxa_from_Dfam_final.csv")
taxa_df<- taxa[,c("TE","Genome")]
taxa_df_new <- separate_rows(taxa_df, Genome, sep = ",")

stars_matrix <- matrix("", nrow = nrow(new_df_ord), ncol = ncol(new_df_ord))
rownames(stars_matrix) <- rownames(new_df_ord)
colnames(stars_matrix) <- colnames(new_df_ord)

for (i in 1:nrow(taxa_df_new)) {
  sp <- taxa_df_new$Genome[i]
  te <- taxa_df_new$TE[i]
  if (sp %in% rownames(new_df_ord_dec) && te %in% colnames(new_df_ord_dec)) {
    stars_matrix[sp, te] <- "*"
  }
}



# col_sums <- colSums(new_df_ord, na.rm = T)
# 
# # decreasing
# col_sums_dec <- sort(col_sums, decreasing = TRUE) 
# TE_copy_num <- data.frame(TE = names(col_sums_dec),
#                           Sum = as.numeric(col_sums_dec),
#                           row.names = NULL)
# new_df_ord_dec <- new_df_ord[, TE_copy_num$TE]

taxa_df_count <- taxa_df_new %>%
  count(TE, name = "TE_count")%>%
  arrange(desc(TE_count))

ordered_colnames <- taxa_df_count$TE

new_df_ord_decTE <- new_df_ord[, ordered_colnames]
stars_matrix_decTE <- stars_matrix[,ordered_colnames]


te_anno_final_3_ordered <- te_anno_final_3[ordered_colnames, , drop = FALSE]

te_table <- openxlsx::read.xlsx("./input3.rmsk_TE_CopyNumber_Subfamily_Family_Class_n1073_v20240110.xlsx")

te_anno <- data.frame(TE = colnames(new_df_ord_decTE))

te_anno_new <- merge(te_anno, te_table, by.x = "TE", by.y = "TE_Subfamily", all.x = T)
te_anno_ord<- te_anno_new[match(colnames(new_df_ord_decTE),te_anno_new$TE),]



te_subfamily <- read.table("./TE_uniq_process/Recurrent_TEs_subfamily_intergenic_n34_uniq_n34.txt")
te_loci <- read.table("./TE_uniq_process/Recurrent_TE_Loci_intergenic_n39_uniq_n25.txt")


te_anno_final <- te_anno_ord %>%
  mutate(TE_Level = case_when(
    TE %in% te_subfamily$V1 & TE %in% te_loci$V1 ~ "Recurrent_TE_Subfamily&Loci",   # n=4
    TE %in% te_subfamily$V1 & !TE %in% te_loci$V1 ~ "Recurrent_TE_Subfamily", # n=34
    !TE %in% te_subfamily$V1 & TE %in% te_loci$V1 ~ "Recurrent_TE_Loci", #n=21
    TRUE ~ NA_character_  
  ))


te_anno_final_2 <- te_anno_final[,c(1,3,4,5)]
rownames(te_anno_final_2) <- NULL
te_anno_final_3 <- column_to_rownames(te_anno_final_2, var = "TE")
te_anno_final_3 <- te_anno_final_3[,c(3,2,1)]



table(te_anno_final_3$TE_Class) # n=4:DNA LINE  LTR SINE
table(te_anno_final_3$TE_Family) # n=12



# annotation_colors <- list(
#   TE_Class = c("DNA" = "#377EB8", "LINE" = "#FFED6F", "LTR" = "#E41A1C","SINE" = "#B2DF8A"),
#   TE_Family =  c("Alu" = "#b3e19b","CR1" = "#FFDA76","ERV1" = "#a2d2e7","ERVK" = "#cdb6da", "ERVL" = "#de149e","ERVL-MaLR" = "#ff9d9f",
#                  "Gypsy" = "#e99b78","hAT" = "#ff8831", "hAT-Charlie" = "#6fb3a8", "L1" = "#67a8cd","L2" = "#50aa4b","MIR" = "#706D54",
#                  "TcMar-Tigger" = "#9a7fbd"),
#   TE_Level = c("Recurrent_TE_Subfamily&Loci"="#7f7625",
#                "Recurrent_TE_Loci"="#f7f2c6",
#                "Recurrent_TE_Subfamily"="#a0c49e")
#   )
# 

annotation_colors <- list(
  TE_Class = c("DNA" = "#f8d196", "LINE" = "#d7a9cb", "LTR" = "#6fa4af", "SINE" = "#8290bb"),
  TE_Family = c(
    "Alu"= "#619c60", "ERV1"= "#a2d2e7", "ERVL"="#dc4aa8",  "ERVL-MaLR"="#ff9d9f",  "Gypsy"= "#e99b78", 
    "hAT"="#dba9a8", "hAT-Blackjack"= "#4d6d7a","hAT-Charlie"="#ff8831","hAT-Tip100"= "#9a7fbd","L1"= "#ca9600", 
    "MIR"= "#706D54",  "TcMar-Tigger"= "#9a7fbd"),
    TE_Level = c("Recurrent_TE_Subfamily&Loci"="#7f7625",
                 "Recurrent_TE_Loci"="#f7f2c6",
                 "Recurrent_TE_Subfamily"="#a0c49e")
  )
  


# V1: TE order by decreasing of taxa sum_number 

library(ComplexHeatmap)
library(circlize)
library(grid)



cat("Start plotting...\n")
pdf("./final_TE_liftOver_Dfam_Heatmap_v1_ordTE_by_TaxaSumDec.pdf", width = 16, height = 10)

top_anno <- HeatmapAnnotation(
  df = te_anno_final_3,                
  col = annotation_colors,             
  annotation_legend_param = list(      
    title_gp = gpar(fontsize = 16,fontface = "bold"),    
    labels_gp = gpar(fontsize = 14)  
  ),
  annotation_name_gp = gpar(col = "black", fontsize = 14) 
)
ht <- Heatmap(
  matrix = new_df_ord_decTE,
  name = "% orthologous",  
  col = colorRamp2(seq(0, 1, length.out = 50), colorRampPalette(c("white", "#2c6fa5"))(50)),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 17),
  column_names_gp = gpar(fontsize = 17),
  column_names_rot = 90,
  border = TRUE,
  top_annotation = top_anno,
  cell_fun = function(j, i, x, y, width, height, fill) {
    txt <- stars_matrix_decTE[i, j]
    if (!is.na(txt) && txt != "") {
      grid.text(txt, x = x, y = y, gp = gpar(fontsize = 14, col = "grey30"))
    }
  },
  heatmap_legend_param = list(
    title = "% orthologous",
    title_gp = gpar(fontsize = 16,fontface = "bold"),
    labels_gp = gpar(fontsize = 14),
    legend_height = unit(0.2, "npc"),
    legend_width = unit(1, "cm")
  )
)

draw(ht,
     heatmap_legend_side = "right",      
     annotation_legend_side = "right"   
)

dev.off()
cat("Plot saved!\n")










# other color test ----------------------------









cat("Start plotting...\n")
pdf("./final_TE_liftOver_Dfam_Heatmap_v1_ordTE_by_TaxaSumDec_purple3.pdf", width = 17, height = 10)


top_anno <- HeatmapAnnotation(
  df = te_anno_final_3,                
  col = annotation_colors,             
  annotation_legend_param = list(      
    title_gp = gpar(fontsize = 16,fontface = "bold"),    
    labels_gp = gpar(fontsize = 14)   
  ),
  annotation_name_gp = gpar(col = "black", fontsize = 14) 
)

ht <- Heatmap(
  matrix = new_df_ord_decTE,
  name = "% orthologous",  
  #col = colorRamp2(seq(0, 1, length.out = 50), colorRampPalette(c("white", "#2c6fa5"))(50)),
  #col = colorRampPalette(c("white", "orange", "darkred"))(100),
  #col = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
  #col = colorRampPalette(c("white", "#b79de0","#8a4f9e"))(100),
  #col = colorRampPalette(c("white", "#d39dc2","#8a4f9e"))(100),
  col = colorRampPalette(c("white", "#e8aebc", "#d96998",  "#8a4f9e"))(100),
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 17),
  column_names_gp = gpar(fontsize = 17),
  column_names_rot = 90,
  border = TRUE,
  top_annotation = top_anno,
  cell_fun = function(j, i, x, y, width, height, fill) {
    txt <- stars_matrix_decTE[i, j]
    if (!is.na(txt) && txt != "") {
      grid.text(txt, x = x, y = y, gp = gpar(fontsize = 14, col = "grey30"))
    }
  },
  heatmap_legend_param = list(
    title = "% orthologous",
    title_gp = gpar(fontsize = 16,fontface = "bold"),
    labels_gp = gpar(fontsize = 14),
    legend_height = unit(0.2, "npc"),
    legend_width = unit(1, "cm")
  )
)

draw(ht,
     heatmap_legend_side = "right",      
     annotation_legend_side = "right"   
)

dev.off()
cat("Plot saved!\n")





cat("Start plotting...\n")
pdf("./final_TE_liftOver_Dfam_Heatmap_v1_ordTE_by_TaxaSumDec_red4.pdf", width = 17, height = 10)

top_anno <- HeatmapAnnotation(
  df = te_anno_final_3,                
  col = annotation_colors,             
  annotation_legend_param = list(      
    title_gp = gpar(fontsize = 16,fontface = "bold"),   
    labels_gp = gpar(fontsize = 14)   
  ),
  annotation_name_gp = gpar(col = "black", fontsize = 14) 
)

ht <- Heatmap(
  matrix = new_df_ord_decTE,
  name = "% orthologous",  
  col = colorRampPalette(c("white", "#ebcce2","#c87d98"))(50),
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 17),
  column_names_gp = gpar(fontsize = 17),
  column_names_rot = 90,
  border = TRUE,
  top_annotation = top_anno,
  cell_fun = function(j, i, x, y, width, height, fill) {
    txt <- stars_matrix_decTE[i, j]
    if (!is.na(txt) && txt != "") {
      grid.text(txt, x = x, y = y, gp = gpar(fontsize = 14, col = "grey30"))
    }
  },
  heatmap_legend_param = list(
    title = "% orthologous",
    title_gp = gpar(fontsize = 16,fontface = "bold"),
    labels_gp = gpar(fontsize = 14),
    legend_height = unit(0.2, "npc"),
    legend_width = unit(1, "cm")
  )
)

draw(ht,
     heatmap_legend_side = "right",     
     annotation_legend_side = "right" 
)

dev.off()
cat("Plot saved!\n")




cat("Start plotting...\n")
pdf("./final_TE_liftOver_Dfam_Heatmap_v1_ordTE_by_TaxaSumDec_green1.pdf", width = 17, height = 10)

top_anno <- HeatmapAnnotation(
  df = te_anno_final_3,                
  col = annotation_colors,           
  annotation_legend_param = list(      
    title_gp = gpar(fontsize = 16,fontface = "bold"),    
    labels_gp = gpar(fontsize = 14)  
  ),
  annotation_name_gp = gpar(col = "black", fontsize = 14) 
)
ht <- Heatmap(
  matrix = new_df_ord_decTE,
  name = "% orthologous",  
  
  col = colorRampPalette(c("white", "#669877"))(50),
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 17),
  column_names_gp = gpar(fontsize = 17),
  column_names_rot = 90,
  border = TRUE,
  top_annotation = top_anno,
  cell_fun = function(j, i, x, y, width, height, fill) {
    txt <- stars_matrix_decTE[i, j]
    if (!is.na(txt) && txt != "") {
      grid.text(txt, x = x, y = y, gp = gpar(fontsize = 14, col = "grey30"))
    }
  },
  heatmap_legend_param = list(
    title = "% orthologous",
    title_gp = gpar(fontsize = 16,fontface = "bold"),
    labels_gp = gpar(fontsize = 14),
    legend_height = unit(0.2, "npc"),
    legend_width = unit(1, "cm")
  )
)

draw(ht,
     heatmap_legend_side = "right",      
     annotation_legend_side = "right"   
)

dev.off()
cat("Plot saved!\n")







cat("Start plotting...\n")
pdf("./final_TE_liftOver_Dfam_Heatmap_v1_ordTE_by_TaxaSumDec_blue1.pdf", width = 17, height = 10)

top_anno <- HeatmapAnnotation(
  df = te_anno_final_3,               
  col = annotation_colors,             
  annotation_legend_param = list(      
    title_gp = gpar(fontsize = 16,fontface = "bold"),  
    labels_gp = gpar(fontsize = 14)   
  ),
  annotation_name_gp = gpar(col = "black", fontsize = 14) 
)
ht <- Heatmap(
  matrix = new_df_ord_decTE,
  name = "% orthologous",  
  
  col = colorRampPalette(c("white", "#4f6b9f"))(50),

  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 17),
  column_names_gp = gpar(fontsize = 17),
  column_names_rot = 90,
  border = TRUE,
  top_annotation = top_anno,
  cell_fun = function(j, i, x, y, width, height, fill) {
    txt <- stars_matrix_decTE[i, j]
    if (!is.na(txt) && txt != "") {
      grid.text(txt, x = x, y = y, gp = gpar(fontsize = 14, col = "grey30"))
    }
  },
  heatmap_legend_param = list(
    title = "% orthologous",
    title_gp = gpar(fontsize = 16,fontface = "bold"),
    labels_gp = gpar(fontsize = 14),
    legend_height = unit(0.2, "npc"),
    legend_width = unit(1, "cm")
  )
)

draw(ht,
     heatmap_legend_side = "right",      
     annotation_legend_side = "right"   
)

dev.off()
cat("Plot saved!\n")

