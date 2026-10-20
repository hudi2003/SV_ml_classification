#!/usr/bin/env Rscript
# 单独绘制 FST 混淆矩阵
# 用法：在已有 00_functions.R 和数据的目录下运行
# Rscript plot_fst_confusion.R

library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)

# 你自己的函数库（需包含 stratified_folds, rank_this_fold, fit_predict）
source("00_functions.R")

# ---------- 1. 读取数据 ----------
geno_mat   <- readRDS("geno_mat.rds")
breed_vec  <- as.factor(readRDS("breed_vec.rds"))
snp_names  <- readRDS("snp_names.rds")
breeds     <- sort(levels(breed_vec))

# ---------- 2. 加载嵌套 CV 结果，只为选择 FST 下的最优分类器 ----------
res <- fread("nested_cv_full_summary.csv")
res_fst <- res[metric == "FST"]

# 选择 >=95% 准确率的最小标记数组合；若无则取全局最高准确率
fst_best <- res_fst[mean_acc >= 0.95]
if (nrow(fst_best) == 0) {
  fst_best <- res_fst[which.max(mean_acc)]
} else {
  fst_best <- fst_best[which.min(n_markers)]
}

CONF_METRIC     <- "FST"
CONF_CLASSIFIER <- fst_best$classifier[1]
CONF_PANEL_SIZES <- c(100, 200, 300, 400, 500, 600, 1000, 1200)  # 根据需要调整

cat(sprintf("\n绘制 FST 混淆矩阵：分类器 = %s，面板大小 = %s\n",
            CONF_CLASSIFIER, paste(CONF_PANEL_SIZES, collapse = ",")))

# ---------- 3. 生成分层 CV 折 ----------
set.seed(2024)  # 固定随机种子，保证可重复
n_folds <- 5
folds <- stratified_folds(breed_vec, n_folds)

# ---------- 4. 按面板大小循环绘制混淆矩阵 ----------
plot_list <- list()

for (n in CONF_PANEL_SIZES) {
  cat(sprintf("\n=== 面板大小: %d 个 SNP ===\n", n))
  
  all_true <- c()
  all_pred <- c()
  
  for (k in 1:n_folds) {
    train_idx <- which(folds != k)
    test_idx  <- which(folds == k)
    
    train_geno  <- geno_mat[train_idx, , drop = FALSE]
    train_breed <- droplevels(breed_vec[train_idx])
    test_geno   <- geno_mat[test_idx, , drop = FALSE]
    test_breed  <- breed_vec[test_idx]
    
    # FST 排序
    ranking_list <- rank_this_fold(train_geno, train_breed, snp_names)
    ranked_snps  <- ranking_list[[CONF_METRIC]]   # 这里就是 FST
    top_snps     <- head(ranked_snps, n)
    
    tr_x <- train_geno[, top_snps, drop = FALSE]
    te_x <- test_geno[,  top_snps, drop = FALSE]
    
    pred <- fit_predict(CONF_CLASSIFIER, tr_x, train_breed, te_x)
    
    all_true <- c(all_true, as.character(test_breed))
    all_pred <- c(all_pred, as.character(pred))
  }
  
  overall_acc <- mean(all_true == all_pred)
  cat(sprintf("总体准确率: %.4f\n", overall_acc))
  
  # 构建混淆矩阵（比例）
  conf_mat  <- table(True = factor(all_true, levels = breeds),
                     Pred = factor(all_pred, levels = breeds))
  conf_prop <- sweep(conf_mat, 1, rowSums(conf_mat), "/")
  conf_df   <- as.data.frame(conf_prop)
  colnames(conf_df) <- c("True", "Pred", "Proportion")
  
  conf_df$label     <- sprintf("%.2f", conf_df$Proportion)
  conf_df$txt_color <- ifelse(conf_df$Proportion > 0.5, "white", "gray20")
  
  p <- ggplot(conf_df, aes(x = Pred, y = True, fill = Proportion)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = label, color = txt_color), size = 2.6, show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_gradient2(low = "white", mid = "#85B7EB", high = "#0C447C",
                         midpoint = 0.5, limits = c(0, 1),
                         name = "Proportion") +
    labs(x = "Predicted breed", y = "True breed",
         title = sprintf("%d markers (FST + %s)", n, CONF_CLASSIFIER),
         subtitle = sprintf("Accuracy: %.1f%%", overall_acc * 100)) +
    theme_classic(base_size = 20) +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y     = element_text(size = 14),
          plot.title      = element_text(size = 18, face = "bold"),
          plot.subtitle   = element_text(size = 14),
          legend.key.size = unit(1, "cm"),
          legend.title    = element_text(size = 16),
          legend.text     = element_text(size = 10))
  
  plot_list[[as.character(n)]] <- p
}

# ---------- 5. 保存图片 ----------
for (n in CONF_PANEL_SIZES) {
  p <- plot_list[[as.character(n)]]
  ggsave(sprintf("fig_confusion_%dSNPs_FST.pdf", n), p, width = 9, height = 8)
  ggsave(sprintf("fig_confusion_%dSNPs_FST.png", n), p, width = 9, height = 8, dpi = 300)
  cat(sprintf("已保存: fig_confusion_%dSNPs_FST.pdf / .png\n", n))
}

cat("\n所有 FST 混淆矩阵图绘制完成。\n")
