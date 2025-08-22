setwd("C:/Users/28257/CuiLab Dropbox/WEI Xiaoman/pathogen_TE_2025_New/10.TE_GeneHancer/02.All_Recurrent_Up_TE_vs_GennHancerRE/")
setwd('/Users/hushu/CuiLab Dropbox/Hu Shu/pathogen_TE_2025_New/10.TE_GeneHancer/02.All_Recurrent_Up_TE_vs_GennHancerRE')


############# 1.  fisher test ######################################################

df_gt0bp <- read.csv("./TE_vs_RE_gt0bp_results_v2.csv", row.names = 1)
df_gt0.5TE <-  read.csv("./TE_vs_RE_gt0.5TE_results_v2.csv", row.names = 1)
df_gt1TE <-  read.csv("./TE_vs_RE_gt1TE_results_v2.csv", row.names = 1)

mat_gt0bp <- as.matrix(df_gt0bp[,-3])
mat_gt0.5TE <- as.matrix(df_gt0.5TE[,-3])
mat_gt1TE <- as.matrix(df_gt1TE[,-3])

fisher.test(mat_gt0bp)
fisher.test(mat_gt0.5TE)
fisher.test(mat_gt1TE)

############# 2. plot ######################################################

library(ggplot2)

# df <- data.frame(
#   Test = c("Any", "50%", "100%"),
#   OR = c(1.298635, 1.39496, 1.462387),
#   CI_lower = c(1.250346, 1.340602, 1.402349),
#   CI_upper = c(1.348410, 1.451088, 1.524539),
#   pval = c("<2.2e-16", "<2.2e-16","<2.2e-16")
# )

df <- data.frame(
  Test = c("Any", "50%", "100%"),
  OR = c(1.313179, 1.408556, 1.48338),
  CI_lower = c(1.271177, 1.359974, 1.428468),
  CI_upper = c(1.356342, 1.458573, 1.540080),
  pval = c("<2.2e-16", "<2.2e-16","<2.2e-16")
)

df$Test <- factor(df$Test, levels = c("Any", "50%", "100%"))


colors <- c("Any" = "#7fb4e0",       
            "50%" = "#3884d4",  
            "100%" = "#083d73") 

p2 <- ggplot(df, aes(x = Test, y = OR, color = Test)) +
  geom_point(size = 4) +
  scale_color_manual(values = colors) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.15,  size=1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", size=0.8) +
  geom_text(aes(label = paste0("p", pval)), vjust =-4, size = 4.2) +
  ylab("Enrichment Ratio\n(recurrent up-TE/other TE)") +
  xlab("Minimum Fraction of TE overlap with Regulatory Element") +
  #xlab("Type 0f TE overlap with Regulatory elemnets") +
  ggtitle("Enrichment of Recurrent Up-TEs\nin Regulatory Elements") +
  scale_y_continuous(limits = c(0.9, 1.6), breaks = seq(0.9, 1.6, by = 0.1)) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(size = 0.7),              
    axis.ticks = element_line(size = 1.2),             
    axis.title = element_text(size = 12), 
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),  
    axis.text.x = element_text(hjust = 1),
    legend.position = "none"  
  ) 

p2

ggsave("ruTEs_enrichment_fisher_test_forest_plot_v2.pdf", p2, width = 5, height = 4.5)
# ggsave("ruTEs_enrichment_fisher_test_forest_plot_no0.pdf", p2, width = 5, height = 4.5)














