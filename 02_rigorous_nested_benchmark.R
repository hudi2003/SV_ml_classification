# ---------------------------------------------------------------------------
# Nested cross-validation benchmark with inner hyper-parameter search
#
# Outer loop : 5-fold stratified CV (seed 42) on the 339 samples.
# Inner loop : 5-fold stratified CV *inside each outer training fold*; a grid
#              search selects the hyper-parameters, the model is then refit on
#              the whole outer training fold and applied once to the outer
#              test fold.
# Marker ranking (FST / Delta / AED / PIC) is recomputed on the outer training
# fold, and again inside every inner training split for the tuning step, so
# neither the outer test samples nor the inner validation samples influence
# the ranking or the hyper-parameters.
#
# Set NESTED_QUICK=1 for a small smoke test, NESTED_CORES=<n> for workers.
# ---------------------------------------------------------------------------
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
library(foreach)
library(doParallel)

source("00_functions.R")

QUICK_TEST <- identical(Sys.getenv("NESTED_QUICK", "0"), "1")
N_CORES    <- as.integer(Sys.getenv("NESTED_CORES", "32"))
INNER_K    <- as.integer(Sys.getenv("NESTED_INNER_K", "5"))
SEED_OUTER <- 42L

cat(format(Sys.time(), "%H:%M:%S"), "Loading data...\n")
geno_mat  <- readRDS("geno_mat.rds")
breed_vec <- as.factor(readRDS("breed_vec.rds"))
snp_names <- readRDS("snp_names.rds")

# Exclude break-end (BND) SVs from the ranking/classification panels, in line
# with the manuscript (BNDs are retained in the callset but excluded from
# population analyses). BND IDs carry a bracketed ALT.
bnd_idx <- grepl("\\[|\\]", snp_names)
if (any(bnd_idx)) {
  geno_mat  <- geno_mat[, !bnd_idx, drop = FALSE]
  snp_names <- snp_names[!bnd_idx]
}
cat(sprintf("Samples: %d | sites: %d | breeds: %d\n",
            nrow(geno_mat), ncol(geno_mat), nlevels(breed_vec)))

panel_sizes_full <- c(10, 20, 50, 100, 200, 300, 400, 500, 600, 700, 800,
                       900, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000)
panel_sizes <- if (QUICK_TEST) c(20, 100) else panel_sizes_full

classifier_list <- c("RF", "SVM_linear", "SVM_RBF", "NB", "LDA",
                      "KNN", "DecisionTree", "XGBoost", "MLP", "Lasso")
metric_list <- c("FST", "Delta", "AED", "PIC", "Random")

n_folds <- if (QUICK_TEST) 2L else 5L

folds <- stratified_folds(breed_vec, n_folds, seed = SEED_OUTER)
cat("Samples per fold and breed:\n")
print(table(breed_vec, folds))

if (is.na(N_CORES) || N_CORES < 1) N_CORES <- 1L
cat(sprintf("Using %d cores | inner folds: %d | quick test: %s\n",
            N_CORES, INNER_K, QUICK_TEST))
registerDoParallel(cores = N_CORES)

detailed <- data.frame()
selected <- data.frame()
t_start  <- Sys.time()

for (k in 1:n_folds) {
  cat(sprintf("\n========== Outer fold %d/%d ==========\n", k, n_folds))
  cat(format(Sys.time(), "%H:%M:%S"), "Recomputing marker ranking on the outer training set...\n")

  train_idx <- which(folds != k)
  test_idx  <- which(folds == k)

  train_geno  <- geno_mat[train_idx, , drop = FALSE]
  train_breed <- droplevels(breed_vec[train_idx])
  test_geno   <- geno_mat[test_idx, , drop = FALSE]
  test_breed  <- breed_vec[test_idx]

  ranking_this_fold <- rank_this_fold(train_geno, train_breed, snp_names,
                                      rng_seed = 100000L + k)

  cat(format(Sys.time(), "%H:%M:%S"),
      sprintf("Building %d inner folds and inner rankings...\n", INNER_K))
  inner_setup <- make_inner_setup(train_geno, train_breed, k = INNER_K,
                                  seed = 2024L + k, snp_names = snp_names)

  combos <- expand.grid(metric = metric_list,
                        n_markers = panel_sizes,
                        classifier = classifier_list,
                        stringsAsFactors = FALSE)
  combos <- combos[combos$n_markers <= length(ranking_this_fold[[1]]), ]

  set.seed(42 + k)
  combos <- combos[sample(nrow(combos)), ]

  cat(sprintf("Evaluating %d combinations...\n", nrow(combos)))

  fold_results <- foreach(i = seq_len(nrow(combos)),
                          .combine = rbind,
                          .packages = c("randomForest","e1071","MASS","class",
                                        "rpart","xgboost","nnet","glmnet"),
                          .export = c("fit_predict", "inner_cv_select",
                                      "nested_fit_predict", "hyper_grid",
                                      "params_to_string", "%||%",
                                      "train_geno", "train_breed",
                                      "test_geno", "test_breed",
                                      "ranking_this_fold", "inner_setup",
                                      "metric_list", "k"),
                          .options.multicore = list(preschedule = FALSE),
                          .errorhandling = "remove") %dopar% {
    metric <- combos$metric[i]
    n      <- combos$n_markers[i]
    cl     <- combos$classifier[i]

    seed_i <- 42L + 10000L * k +
              100L * match(metric, metric_list) +
              10L * n + match(cl, c("RF","SVM_linear","SVM_RBF","NB","LDA",
                                    "KNN","DecisionTree","XGBoost","MLP","Lasso"))

    t0 <- Sys.time()
    fit <- tryCatch(
      nested_fit_predict(cl, metric, n, train_geno, train_breed, test_geno,
                         ranking_this_fold, inner_setup, seed = seed_i),
      error = function(e) list(pred = rep(NA_character_, length(test_breed)),
                               params = list(), inner_acc = NA_real_,
                               inner_sd = NA_real_, searched = NA,
                               grid_size = NA_integer_,
                               error = conditionMessage(e))
    )
    t1 <- Sys.time()

    acc <- mean(fit$pred == as.character(test_breed), na.rm = TRUE)
    err <- if (!is.null(fit$error)) fit$error else NA_character_

    data.frame(fold = k, metric = metric, n_markers = n, classifier = cl,
               acc = acc, time_sec = as.numeric(difftime(t1, t0, units = "secs")),
               params = params_to_string(fit$params),
               inner_acc = fit$inner_acc, inner_sd = fit$inner_sd,
               inner_grid_size = fit$grid_size,
               status = ifelse(is.na(err), "ok", "error"),
               error = err,
               stringsAsFactors = FALSE)
  }

  if (nrow(fold_results) != nrow(combos)) {
    cat(sprintf("WARNING: %d/%d combinations returned a result in fold %d\n",
                nrow(fold_results), nrow(combos), k))
  }
  n_err <- sum(fold_results$status != "ok", na.rm = TRUE)
  if (n_err > 0) cat(sprintf("WARNING: %d combinations failed in fold %d\n", n_err, k))

  detailed <- rbind(detailed, fold_results)
  cat(sprintf("[Fold %d] done, %d records, elapsed %.1f min\n", k, nrow(fold_results),
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

  fwrite(detailed, "nested_cv_full_detail.csv")
  rm(fold_results); gc()
}

stopImplicitCluster()

detailed$params[detailed$params == ""] <- "fixed"

results <- detailed %>%
  dplyr::group_by(metric, n_markers, classifier) %>%
  dplyr::summarise(mean_acc  = mean(acc, na.rm = TRUE),
                   sd_acc    = sd(acc, na.rm = TRUE),
                   mean_time = mean(time_sec, na.rm = TRUE),
                   n_ok      = sum(status == "ok"),
                   .groups = "drop")

sel_tab <- detailed %>%
  dplyr::filter(status == "ok") %>%
  dplyr::select(fold, metric, n_markers, classifier, params,
                inner_acc, inner_sd, inner_grid_size)

fwrite(results,  "nested_cv_full_summary.csv")
fwrite(detailed, "nested_cv_full_detail.csv")
fwrite(sel_tab,  "nested_cv_selected_params.csv")

writeLines(c(
  sprintf("run_time = %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("outer_folds = %d", n_folds),
  sprintf("outer_seed = %d", SEED_OUTER),
  sprintf("inner_folds = %d", INNER_K),
  sprintf("panel_sizes = %s", paste(panel_sizes, collapse = ",")),
  sprintf("metrics = %s", paste(metric_list, collapse = ",")),
  sprintf("classifiers = %s", paste(classifier_list, collapse = ",")),
  sprintf("RF_grid = ntree{300,500,1000} x mtry{floor(sqrt(p)),floor(p/5),floor(p/2)}"),
  sprintf("SVM_linear_grid = cost{0.01,0.1,1,10}"),
  sprintf("cores = %d", N_CORES),
  sprintf("quick_test = %s", QUICK_TEST)
), "nested_cv_run_info.txt")

cat(format(Sys.time(), "%H:%M:%S"), "Nested benchmark finished.\n")
cat(sprintf("Total elapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat("Summary: nested_cv_full_summary.csv | Detail: nested_cv_full_detail.csv\n")
cat("Selected hyper-parameters: nested_cv_selected_params.csv\n")
