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

source("00_functions.R")

QUICK_TEST <- FALSE
ADD_RF_MDA <- FALSE

cat(format(Sys.time(), "%H:%M:%S"), "Loading data...\n")
geno_mat  <- readRDS("geno_mat.rds")
breed_vec <- as.factor(readRDS("breed_vec.rds"))
snp_names <- readRDS("snp_names.rds")
cat(sprintf("Samples: %d | sites: %d | breeds: %d\n",
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
cat("Samples per fold and breed:\n")
print(table(breed_vec, folds))

n_cores <- 32

if (is.na(n_cores) || n_cores < 1) n_cores <- 1
cat(sprintf("Using %d cores\n", n_cores))
registerDoParallel(cores = n_cores)

results  <- data.frame()
detailed <- data.frame()

for (k in 1:n_folds) {
  cat(sprintf("\n========== Fold %d/%d ==========\n", k, n_folds))
  cat(format(Sys.time(), "%H:%M:%S"), "Recomputing marker ranking on the training set...\n")

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
    cat(format(Sys.time(), "%H:%M:%S"), "Training RF for MDA importance...\n")
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

  set.seed(42 + k)
  combos <- combos[sample(nrow(combos)), ]

  cat(sprintf("Evaluating %d combinations...\n", nrow(combos)))

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
  cat(sprintf("[Fold %d] done, %d records\n", k, nrow(fold_results)))
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
cat(format(Sys.time(), "%H:%M:%S"), "Nested benchmark finished.\n")
cat("Summary: nested_cv_full_summary.csv\n")
cat("Per-fold detail: nested_cv_full_detail.csv\n")
