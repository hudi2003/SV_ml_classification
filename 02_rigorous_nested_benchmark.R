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

# --- parallel-decomposition support -----------------------------------------
# NESTED_FOLDS        "3" or "1,2"  -> run only these outer folds
# NESTED_CLASSIFIERS  "RF,SVM_linear" -> run only these classifiers
# NESTED_TAG          "f3_A"        -> suffix for every output file
# Folds and inner splits are always constructed with the same seeds, so the
# results of a split run are identical to the corresponding rows of a full run.
CLS_ENV <- Sys.getenv("NESTED_CLASSIFIERS", "")
if (nzchar(CLS_ENV)) {
  wanted <- trimws(strsplit(CLS_ENV, ",")[[1]])
  stopifnot(all(wanted %in% classifier_list))
  classifier_list <- classifier_list[classifier_list %in% wanted]
}
TAG <- Sys.getenv("NESTED_TAG", "full")
DETAIL_FILE  <- sprintf("nested_cv_detail_%s.csv", TAG)
SUMMARY_FILE <- sprintf("nested_cv_summary_%s.csv", TAG)
PARAMS_FILE  <- sprintf("nested_cv_params_%s.csv", TAG)
INFO_FILE    <- sprintf("nested_cv_runinfo_%s.txt", TAG)
metric_list <- c("FST", "Delta", "AED", "PIC", "Random")

n_folds <- if (QUICK_TEST) 2L else 5L
FOLD_ENV  <- Sys.getenv("NESTED_FOLDS", "")
fold_list <- if (nzchar(FOLD_ENV)) as.integer(strsplit(FOLD_ENV, ",")[[1]]) else
               seq_len(n_folds)
stopifnot(all(fold_list %in% seq_len(n_folds)))

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

# Resume support: combinations already present in the detail table are skipped,
# so a long run can be restarted without repeating finished work.
if (file.exists(DETAIL_FILE)) {
  detailed <- fread(DETAIL_FILE)
  cat(sprintf("Resuming from %s: %d records already done\n", DETAIL_FILE,
              nrow(detailed)))
} else {
  detailed <- data.frame()
}
CHUNK <- as.integer(Sys.getenv("NESTED_CHUNK", "25"))
CLS_ORDER <- c("RF", "SVM_linear", "SVM_RBF", "NB", "LDA",
               "KNN", "DecisionTree", "XGBoost", "MLP", "Lasso")

for (k in fold_list) {
  cat(sprintf("\n========== Outer fold %d (of %d) ==========\n", k, n_folds))
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
  total_fold <- nrow(combos)

  if (nrow(detailed)) {
    done_key <- paste(detailed$fold, detailed$metric, detailed$n_markers,
                      detailed$classifier, sep = "|")
    my_key <- paste(k, combos$metric, combos$n_markers, combos$classifier, sep = "|")
    combos <- combos[!(my_key %in% done_key), , drop = FALSE]
    cat(sprintf("Skipping combinations already in the detail table (%d left)\n",
                nrow(combos)))
  }

  cat(sprintf("Evaluating %d combinations in chunks of %d...\n", nrow(combos), CHUNK))
  fold_start <- Sys.time()
  chunks <- split(seq_len(nrow(combos)), ceiling(seq_len(nrow(combos)) / CHUNK))

  for (ch in chunks) {
    cc <- combos[ch, , drop = FALSE]
    # Run the expensive combinations first inside a chunk so that the slow jobs
    # (MLP / RF / XGBoost on large panels) do not form a long single-core tail.
    # Per-combination seeds are fixed, so the ordering changes no result.
    cost_w <- c(RF = 6, SVM_linear = 1, SVM_RBF = 1.5, NB = 0.5, LDA = 0.2,
                KNN = 0.2, DecisionTree = 1, XGBoost = 2, MLP = 30, Lasso = 1)
    cc <- cc[order(-cc$n_markers * cost_w[cc$classifier]), , drop = FALSE]
    chunk_results <- foreach(j = seq_len(nrow(cc)),
                             .combine = rbind,
                             .packages = c("randomForest","e1071","MASS","class",
                                           "rpart","xgboost","nnet","glmnet"),
                             .export = c("fit_predict", "inner_cv_select",
                                         "nested_fit_predict", "hyper_grid",
                                         "params_to_string", "%||%",
                                         "train_geno", "train_breed",
                                         "test_geno", "test_breed",
                                         "ranking_this_fold", "inner_setup",
                                         "metric_list", "CLS_ORDER", "k", "cc"),
                             .options.multicore = list(preschedule = FALSE),
                             .errorhandling = "remove") %dopar% {
      metric <- cc$metric[j]
      n      <- cc$n_markers[j]
      cl     <- cc$classifier[j]

      seed_i <- 42L + 10000L * k +
                100L * match(metric, metric_list) +
                10L * n + match(cl, CLS_ORDER)

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

    if (nrow(chunk_results) != nrow(cc)) {
      cat(sprintf("WARNING: %d/%d combinations returned a result\n",
                  nrow(chunk_results), nrow(cc)))
    }
    detailed <- rbind(detailed, chunk_results)
    fwrite(detailed, DETAIL_FILE)
    done_fold <- sum(detailed$fold == k)
    el <- as.numeric(difftime(Sys.time(), fold_start, units = "mins"))
    cat(sprintf("  [fold %d] %d/%d combos | %.1f min elapsed | ~%.1f min left in fold\n",
                k, done_fold, total_fold, el,
                if (done_fold > 0) el * (total_fold - done_fold) / done_fold else NA))
    rm(chunk_results); gc()
  }

  cat(sprintf("[Fold %d] done, elapsed %.1f min total\n", k,
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
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

fwrite(results,  SUMMARY_FILE)
fwrite(detailed, DETAIL_FILE)
fwrite(sel_tab,  PARAMS_FILE)

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
  sprintf("SVM_RBF_grid = cost{0.1,1,10} x gamma{1/2p,1/p,2/p}"),
  sprintf("NB_grid = laplace{0,1}"),
  sprintf("KNN_grid = k{1,3,5,7,9}"),
  sprintf("DecisionTree_grid = cp{0.001,0.005,0.01,0.05}"),
  sprintf("XGBoost_grid = nrounds{100,200} x max_depth{3,6} x eta{0.3}"),
  sprintf("MLP_grid = size{5,10,20} x decay{5e-4,5e-3}"),
  sprintf("LDA_grid = none (no hyper-parameters); Lasso = internal cv.glmnet"),
  sprintf("cores = %d", N_CORES),
  sprintf("quick_test = %s", QUICK_TEST)
), INFO_FILE)

cat(format(Sys.time(), "%H:%M:%S"), "Nested benchmark finished.\n")
cat(sprintf("Total elapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat(sprintf("Summary: %s | Detail: %s\n", SUMMARY_FILE, DETAIL_FILE))
cat(sprintf("Selected hyper-parameters: %s\n", PARAMS_FILE))
