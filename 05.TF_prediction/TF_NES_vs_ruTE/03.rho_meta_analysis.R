# ============================================================
# TF-TE 相关性元分析
# Meta-analysis of correlation coefficients
# ============================================================

library(tidyverse)
library(metafor)  # 专业元分析包 v4.9-29
setwd('/Users/hushu/Dropbox (个人)/Hu Shu/pathogen_TE_2025_New/08.TF_prediction/TF_prediction_add_bg_v3/TF_NES_vs_ruTE/TF_activity_vs_ruTE_correlation_UpOnly/')
# 读取数据
df <- read_csv('/Users/hushu/Dropbox (个人)/Hu Shu/pathogen_TE_2025_New/08.TF_prediction/TF_prediction_add_bg_v3/TF_NES_vs_ruTE/TF_activity_vs_ruTE_correlation_UpOnly/per_dataset_TF_TE_correlation.csv')
# ============================================================
# 使用 metafor 包
# ============================================================
results_metafor <- df %>%
  group_by(TF, TE) %>%
  group_modify(~ {
    # escalc: 计算 Fisher z 和方差
    es <- escalc(measure = "ZCOR", ri = .x$Spearman_rho, ni = .x$N_samples)
    # rma: 随机效应模型
    model <- rma(yi, vi, data = es, method = "REML")
    # 反变换回相关系数尺度
    pred <- predict(model, transf = transf.ztor)
    tibble(
      `#up-TE datasets` = nrow(es),
      `#rho>0 & p<0.05` = sum(.x$Spearman_rho > 0 & .x$P_value < 0.05),
      Total_samples = sum(.x$N_samples),
      Meta_rho = pred$pred,
      CI_lower = pred$ci.lb,
      CI_upper = pred$ci.ub,
      Meta_p = model$pval,
      I2_percent = model$I2
    )
  }) %>%
  arrange(desc(Meta_rho))

print(results_metafor)
# 保存结果
write_csv(results_metafor, "meta_analysis_results.csv")

# ============================================================
# 森林图
# ============================================================
plot_forest <- function(tf, te) {
  sub_df <- df %>% filter(TF == tf, TE == te)
  
  es <- escalc(
    measure = "ZCOR",
    ri = sub_df$Spearman_rho,
    ni = sub_df$N_samples,
    slab = sub_df$Dataset
  )
  
  model <- rma(yi, vi, data = es, method = "REML")
  
  forest(model, 
         atransf = transf.ztor,  # 显示为相关系数
         xlab = "Correlation (ρ)",
         main = paste(tf, "-", te))
}

plot_forest("RELA", "LTR27E_intergenic")
plot_forest("RELA", "HERVIP10FH-int_intergenic")
