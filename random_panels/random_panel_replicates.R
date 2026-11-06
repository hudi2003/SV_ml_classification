# ---------------------------------------------------------------------------
# Random-panel control for the ranked SV panels (reviewer Major 5)
#
# Repeats the nested procedure with independently drawn, unordered panels of the
# same size:
#   outer : 5-fold stratified CV (seed 42) on the 339 samples
#   inner : 5-fold stratified CV inside each outer training fold, grid search
#   model : SVM linear, cost in {0.01, 0.1, 1, 10}
#   panel : RAND_N SVs
#
# Each replicate draws its own panel for the outer training fold and for every
# inner training split, so the only difference from a ranked panel is how the
# markers are chosen.
#
# Environment variables
#   RAND_N         panel size                                  (default 400)
#   RAND_REP       number of independent draws                 (default 1000)
#   NESTED_CORES   workers                                     (default 32)
#   NESTED_INNER_K inner folds                                 (default 5)
#
# Output, prefix random<RAND_N>:
#   <prefix>_detail.csv   one row per draw x outer fold
#   <prefix>_summary.txt  draw-level mean, SD and percentiles
# Significance testing is done by stratified_monte_carlo.py.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(e1071); library(foreach)
  library(doParallel)
})
source("00_functions.R")

N_REP    <- as.integer(Sys.getenv("RAND_REP", "1000"))
N_CORES  <- as.integer(Sys.getenv("NESTED_CORES", "32"))
INNER_K  <- as.integer(Sys.getenv("NESTED_INNER_K", "5"))
N_MARK   <- as.integer(Sys.getenv("RAND_N", "400"))
SEED_OUT <- 42L
META     <- "SVM_linear"
PREFIX   <- sprintf("random%d", N_MARK)

geno_mat  <- readRDS("geno_mat.rds")
breed_vec <- as.factor(readRDS("breed_vec.rds"))
snp_names <- readRDS("snp_names.rds")

bnd_idx <- grepl("\\[|\\]", snp_names)
if (any(bnd_idx)) {
  geno_mat  <- geno_mat[, !bnd_idx, drop = FALSE]
  snp_names <- snp_names[!bnd_idx]
}
cat(sprintf("Samples: %d | sites: %d | breeds: %d | markers per panel: %d | reps: %d\n",
            nrow(geno_mat), ncol(geno_mat), nlevels(breed_vec), N_MARK, N_REP))

folds <- stratified_folds(breed_vec, 5L, seed = SEED_OUT)
n_folds <- 5L

registerDoParallel(cores = N_CORES)

t0_all <- Sys.time()
res <- foreach(r = seq_len(N_REP), .combine = rbind,
               .packages = c("e1071", "data.table"),
               .export = c("fit_predict", "inner_cv_select", "nested_fit_predict",
                           "hyper_grid", "params_to_string", "%||%",
                           "geno_mat", "breed_vec", "snp_names", "folds",
                           "INNER_K", "N_MARK", "META", "n_folds",
                           "stratified_folds"),
               .options.multicore = list(preschedule = FALSE)) %dopar% {
  out <- vector("list", n_folds)
  for (k in seq_len(n_folds)) {
    train_idx <- which(folds != k)
    test_idx  <- which(folds == k)
    tr_geno  <- geno_mat[train_idx, , drop = FALSE]
    tr_breed <- droplevels(breed_vec[train_idx])
    te_geno  <- geno_mat[test_idx, , drop = FALSE]
    te_breed <- breed_vec[test_idx]

    # random panel for the outer training fold
    set.seed(300000L + 1000L * r + k)
    perm_out <- snp_names[sample.int(length(snp_names))]
    outer_rank <- list(Random = perm_out)

    # random panel for every inner training split (used only for tuning)
    inner_folds <- stratified_folds(tr_breed, INNER_K, seed = 2024L + k)
    ids <- lapply(seq_len(INNER_K), function(f) which(inner_folds == f))
    ranks <- lapply(seq_len(INNER_K), function(f) {
      set.seed(300000L + 1000L * r + 100L * k + f)
      list(Random = snp_names[sample.int(length(snp_names))])
    })
    setup <- list(folds = inner_folds, ids = ids, rankings = ranks, k = INNER_K)

    t0 <- Sys.time()
    fit <- tryCatch(
      nested_fit_predict(META, "Random", N_MARK, tr_geno, tr_breed, te_geno,
                         outer_rank, setup, seed = 42L + 10000L * k + r),
      error = function(e) list(pred = rep(NA_character_, length(te_breed)),
                               params = list(), inner_acc = NA_real_,
                               inner_sd = NA_real_, grid_size = NA_integer_,
                               error = conditionMessage(e)))
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    acc <- mean(fit$pred == as.character(te_breed), na.rm = TRUE)

    out[[k]] <- data.frame(
      rep = r, fold = k, metric = "Random", n_markers = N_MARK,
      classifier = META, acc = acc,
      params = params_to_string(fit$params),
      inner_acc = fit$inner_acc, inner_sd = fit$inner_sd,
      time_sec = dt,
      status = if (is.null(fit$error)) "ok" else "error",
      error = if (is.null(fit$error)) NA_character_ else fit$error,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}
stopImplicitCluster()

fwrite(res, sprintf("%s_detail.csv", PREFIX))

# ---- summarise -------------------------------------------------------------
ok <- res[res$status == "ok", ]
rep_mean <- ok %>% group_by(rep) %>% summarise(mean_acc = mean(acc),
                                               sd_acc = sd(acc), .groups = "drop")
q <- quantile(rep_mean$mean_acc, c(0.025, 0.05, 0.5, 0.95, 0.975))
summ <- sprintf(paste0(
  "Random %d-SV panels, %d independent draws, nested 5x5, SVM linear\n",
  "  draw-level mean accuracy  : mean %.2f%%  sd %.2f%%  min %.2f%%  max %.2f%%\n",
  "  percentiles of draw mean  : 2.5%% %.2f%% | 5%% %.2f%% | 50%% %.2f%% | 95%% %.2f%% | 97.5%% %.2f%%\n",
  "  fold-level accuracy       : mean %.2f%%  sd %.2f%%  (n=%d fold estimates)\n"),
  N_MARK, N_REP,
  100 * mean(rep_mean$mean_acc), 100 * sd(rep_mean$mean_acc),
  100 * min(rep_mean$mean_acc), 100 * max(rep_mean$mean_acc),
  100 * q[1], 100 * q[2], 100 * q[3], 100 * q[4], 100 * q[5],
  100 * mean(ok$acc), 100 * sd(ok$acc), nrow(ok))

writeLines(summ, sprintf("%s_summary.txt", PREFIX))
cat(paste(summ, collapse = "\n"), "\n")
cat(sprintf("Elapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t0_all, units = "mins"))))
