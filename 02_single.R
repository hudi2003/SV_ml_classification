# Run ONE (fold, metric, panel size, classifier) cell of the nested benchmark.
# Intended for a 1-core job array element: Rscript 02_single.R <task_index>
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(randomForest); library(e1071)
  library(MASS); library(class); library(rpart); library(xgboost); library(nnet)
  library(glmnet); library(foreach); library(doParallel)
})
source("00_functions.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript 02_single.R <task_index>")
idx <- as.integer(args[1])

tasks <- fread("tasks.csv")
if (is.na(idx) || idx < 1 || idx > nrow(tasks)) stop("task index out of range")
tk <- tasks[idx]

OUTDIR <- Sys.getenv("NESTED_OUTDIR", "parts")
dir.create(OUTDIR, showWarnings = FALSE)
out_file <- sprintf("%s/row_%05d.csv", OUTDIR, idx)
if (file.exists(out_file)) { cat("task", idx, "already done\n"); quit(save = "no", status = 0) }

CLS_ORDER   <- c("RF", "SVM_linear", "SVM_RBF", "NB", "LDA",
                 "KNN", "DecisionTree", "XGBoost", "MLP", "Lasso")
METRIC_ORDER <- c("FST", "Delta", "AED", "PIC", "Random")
N_OUTER <- 5L
INNER_K <- as.integer(Sys.getenv("NESTED_INNER_K", "5"))

geno_mat  <- readRDS("geno_mat.rds")
breed_vec <- as.factor(readRDS("breed_vec.rds"))
snp_names <- readRDS("snp_names.rds")
bnd <- grepl("\\[|\\]", snp_names)
if (any(bnd)) {
  geno_mat  <- geno_mat[, !bnd, drop = FALSE]
  snp_names <- snp_names[!bnd]
}

folds <- stratified_folds(breed_vec, N_OUTER, seed = 42L)
k <- as.integer(tk$fold)
train_idx <- which(folds != k); test_idx <- which(folds == k)
train_geno  <- geno_mat[train_idx, , drop = FALSE]
train_breed <- droplevels(breed_vec[train_idx])
test_geno   <- geno_mat[test_idx, , drop = FALSE]
test_breed  <- breed_vec[test_idx]

ranking_this_fold <- rank_this_fold(train_geno, train_breed, snp_names,
                                    rng_seed = 100000L + k)
method <- tk$classifier
grid <- hyper_grid(method, as.integer(tk$n_markers))
inner_setup <- if (is.null(grid) && !identical(method, "LDA")) NULL else
  make_inner_setup(train_geno, train_breed, k = INNER_K, seed = 2024L + k,
                   snp_names = snp_names)

seed_i <- 42L + 10000L * k +
          100L * match(tk$metric, METRIC_ORDER) +
          10L * as.integer(tk$n_markers) +
          match(method, CLS_ORDER)

t0 <- Sys.time()
fit <- tryCatch(
  nested_fit_predict(method, tk$metric, as.integer(tk$n_markers),
                     train_geno, train_breed, test_geno,
                     ranking_this_fold, inner_setup, seed = seed_i),
  error = function(e) list(pred = rep(NA_character_, length(test_breed)),
                           params = list(), inner_acc = NA_real_,
                           inner_sd = NA_real_, searched = NA, grid_size = NA_integer_,
                           error = conditionMessage(e))
)
t1 <- Sys.time()

acc <- mean(fit$pred == as.character(test_breed), na.rm = TRUE)
err <- if (!is.null(fit$error)) fit$error else NA_character_

res <- data.frame(fold = k, metric = tk$metric, n_markers = as.integer(tk$n_markers),
                  classifier = method, acc = acc,
                  time_sec = as.numeric(difftime(t1, t0, units = "secs")),
                  params = params_to_string(fit$params),
                  inner_acc = fit$inner_acc, inner_sd = fit$inner_sd,
                  inner_grid_size = fit$grid_size,
                  status = ifelse(is.na(err), "ok", "error"),
                  error = err, stringsAsFactors = FALSE)
tmp <- paste0(out_file, ".tmp")
fwrite(res, tmp)
file.rename(tmp, out_file)
cat(sprintf("task %d done: fold=%d %s p=%d %s acc=%.4f time=%.1fs\n",
            idx, k, tk$metric, as.integer(tk$n_markers), method, acc, res$time_sec))
