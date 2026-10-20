#!/usr/bin/env Rscript
# ================================================================
# 脚本：02_rigorous_nested_benchmark.R （内部并行版 + Stacking）
# 功能：彻底消除信息泄漏 + 自动利用多核加速 + 概率级Stacking集成
# 用法：jsub -n 8 -q normal ... 'Rscript 02_rigorous_nested_benchmark.R'
# ================================================================

library(data.table)
library(dplyr)
library(randomForest)
library(e1071)
library(MASS)
library(class)
library(rpart)
library(xgboost)
library(nnet)
library(glmnet)

HAS_LGB <- requireNamespace("lightgbm", quietly = TRUE)
if (HAS_LGB) library(lightgbm)

library(foreach)
library(doParallel)

source("00_functions.R")   # 统一从这里加载所有函数，包括 Stacking

QUICK_TEST <- FALSE     # 先用少量面板大小试跑，确认流程没问题后改 FALSE
ADD_RF_MDA <- FALSE

# ---- 1. 加载数据 ----
cat(format(Sys.time(), "%H:%M:%S"), "加载数据...\n")
geno_mat  <- readRDS("geno_mat.rds")
breed_vec <- as.factor(readRDS("breed_vec.rds"))
snp_names <- readRDS("snp_names.rds")
cat(sprintf("样本数: %d | 位点数: %d | 品种数: %d\n",
            nrow(geno_mat), ncol(geno_mat), nlevels(breed_vec)))

panel_sizes_full <- c(10, 20, 50, 100, 200, 300, 400, 500, 600, 700, 800,
                       900, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000)
panel_sizes <- if (QUICK_TEST) c(50, 200, 500) else panel_sizes_full

classifier_list <- c("RF", "SVM_linear", "SVM_RBF", "NB", "LDA",
                      "KNN", "DecisionTree", "XGBoost", "MLP", "Lasso")
if (HAS_LGB) classifier_list <- c(classifier_list, "LightGBM")

metric_list <- c("FST", "Delta", "AED", "PIC")
if (ADD_RF_MDA) metric_list <- c(metric_list, "RF_MDA")

n_folds <- if (QUICK_TEST) 3 else 5

folds <- stratified_folds(breed_vec, n_folds)
cat("各折各品种样本数检查:\n")
print(table(breed_vec, folds))

# ---- 并行环境初始化 ----
# n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
n_cores <- 32

if (is.na(n_cores) || n_cores < 1) n_cores <- 1
cat(sprintf("使用 %d 个核心进行并行计算\n", n_cores))
registerDoParallel(cores = n_cores)

# ---- 主循环 ----
results  <- data.frame()
detailed <- data.frame()

for (k in 1:n_folds) {
  cat(sprintf("\n========== Fold %d/%d ==========\n", k, n_folds))
  cat(format(Sys.time(), "%H:%M:%S"), "重新计算该折训练集的标记排序...\n")

  train_idx <- which(folds != k)
  test_idx  <- which(folds == k)

  train_geno  <- geno_mat[train_idx, , drop = FALSE]
  train_breed <- droplevels(breed_vec[train_idx])
  test_geno   <- geno_mat[test_idx, , drop = FALSE]
  test_breed  <- breed_vec[test_idx]

  fst_train   <- compute_fst_wc(train_geno, train_breed)
  dap         <- compute_delta_aed_pic(train_geno, train_breed)

  ranking_this_fold <- list(
    FST   = snp_names[order(-fst_train)],
    Delta = snp_names[order(-dap$Delta)],
    AED   = snp_names[order(-dap$AED)],
    PIC   = snp_names[order(-dap$PIC)]
  )

  if (ADD_RF_MDA) {
    cat(format(Sys.time(), "%H:%M:%S"), "训练RF计算该折的MDA重要性...\n")
    rf_imp_model <- randomForest(x = train_geno, y = train_breed,
                                  ntree = 500, importance = TRUE)
    mda <- importance(rf_imp_model, type = 1)[, "MeanDecreaseAccuracy"]
    ranking_this_fold$RF_MDA <- colnames(train_geno)[order(-mda)]
  }

  combos <- expand.grid(metric = metric_list,
                        n_markers = panel_sizes,
                        classifier = classifier_list,
                        stringsAsFactors = FALSE)
  n_max <- length(ranking_this_fold[[1]])
  combos <- combos[combos$n_markers <= n_max, ]
  
  # 新增这一行：打乱行顺序，避免Stacking_proba这类慢任务集中分配到少数几个核
  set.seed(42 + k)  # 每折用不同种子，但保证可复现
  combos <- combos[sample(nrow(combos)), ]

  cat(sprintf("该折待评估 %d 个组合，并行计算中...\n", nrow(combos)))

  fold_results <- foreach(i = seq_len(nrow(combos)),
                          .combine = rbind,
                          .packages = c("randomForest","e1071","MASS","class",
                                        "rpart","xgboost","nnet","glmnet"),
                          .export = c("fit_predict",  
                                      "stratified_folds",
                                      "train_geno","train_breed",
                                      "test_geno","test_breed",
                                      "ranking_this_fold",
                                      "k"),
			  .options.multicore = list(preschedule = FALSE),
                          .errorhandling = "remove") %dopar% {
    metric  <- combos$metric[i]
    n       <- combos$n_markers[i]
    cl      <- combos$classifier[i]

    if (cl == "LightGBM") {
      requireNamespace("lightgbm", quietly = TRUE)
      library(lightgbm)
    }

    top_snps <- ranking_this_fold[[metric]][1:n]
    tr_x <- train_geno[, top_snps, drop = FALSE]
    te_x <- test_geno[,  top_snps, drop = FALSE]

    t0 <- Sys.time()
    pred <- tryCatch(
      fit_predict(cl, tr_x, train_breed, te_x),
      error = function(e) rep(NA_character_, length(test_breed))
    )
    t1 <- Sys.time()
    acc <- mean(pred == as.character(test_breed), na.rm = TRUE)
    elapsed <- as.numeric(difftime(t1, t0, units = "secs"))

    data.frame(fold = k, metric = metric, n_markers = n,
               classifier = cl, acc = acc, time_sec = elapsed,
               stringsAsFactors = FALSE)
  }

  detailed <- rbind(detailed, fold_results)
  cat(sprintf("[Fold %d] 完成，已收集 %d 条记录\n", k, nrow(fold_results)))
}

stopImplicitCluster()

results <- detailed %>%
  dplyr::group_by(metric, n_markers, classifier) %>%
  dplyr::summarise(mean_acc = mean(acc, na.rm = TRUE),
                    sd_acc   = sd(acc, na.rm = TRUE),
                    mean_time = mean(time_sec, na.rm = TRUE),
                    .groups = "drop")

fwrite(results,  "nested_cv_full_summary.csv")
fwrite(detailed, "nested_cv_full_detail.csv")
cat(format(Sys.time(), "%H:%M:%S"), "严谨嵌套并行版（含Stacking）全部完成。\n")
cat("汇总结果: nested_cv_full_summary.csv\n")
cat("每折明细: nested_cv_full_detail.csv\n")
