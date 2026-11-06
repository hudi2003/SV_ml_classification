# ---------------------------------------------------------------------------
# Repeated outer cross-validation: is the ranking gain stable across splits?
#
# The main analysis uses one fixed 5-fold split (seed 42), so every panel has
# only five outer-fold accuracies. Here the whole procedure is repeated on
# N_SPLIT independently drawn, breed-stratified 5-fold splits. For each split we
# evaluate the four ranked panels and N_RAND random panels at each size, using
# exactly the same inner tuning, so the paired difference
#     d_b = acc(ranked, split b) - acc(random, split b)
# is available for 1,000 independent splits instead of five folds.
#
# Environment
#   REP_SPLITS  number of outer splits                 (default 1000)
#   REP_RAND    random panels per split                (default 5)
#   REP_SIZES   panel sizes, comma separated           (default 400,600)
#   NESTED_CORES workers                               (default 48)
#   NESTED_INNER_K inner folds                         (default 5)
#
# Output: repeated_outer_cv_detail.csv and a printed summary.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(foreach); library(doParallel); library(e1071)
})
source("00_functions.R")

N_SPLIT <- as.integer(Sys.getenv("REP_SPLITS", "1000"))
N_RAND  <- as.integer(Sys.getenv("REP_RAND", "5"))
SIZES   <- as.integer(strsplit(Sys.getenv("REP_SIZES", "400,600"), ",")[[1]])
N_CORES <- as.integer(Sys.getenv("NESTED_CORES", "48"))
INNER_K <- as.integer(Sys.getenv("NESTED_INNER_K", "5"))
META    <- "SVM_linear"
METRICS <- c("FST", "Delta", "AED", "PIC")

geno_mat  <- readRDS("geno_mat.rds")
breed_vec <- as.factor(readRDS("breed_vec.rds"))
snp_names <- readRDS("snp_names.rds")
bnd <- grepl("\\[|\\]", snp_names)
if (any(bnd)) {
  geno_mat  <- geno_mat[, !bnd, drop = FALSE]
  snp_names <- snp_names[!bnd]
}
cat(sprintf("splits %d | random panels per split %d | sizes %s | sites %d\n",
            N_SPLIT, N_RAND, paste(SIZES, collapse = ","), length(snp_names)))

registerDoParallel(cores = N_CORES)
t0 <- Sys.time()

res <- foreach(b = seq_len(N_SPLIT), .combine = rbind,
               .packages = c("e1071", "data.table"),
               .export = c("fit_predict", "inner_cv_select", "nested_fit_predict",
                           "hyper_grid", "params_to_string", "%||%",
                           "geno_mat", "breed_vec", "snp_names",
                           "INNER_K", "META", "METRICS", "SIZES", "N_RAND",
                           "stratified_folds", "rank_this_fold",
                           "make_inner_setup", "compute_fst_wc"),
               .options.multicore = list(preschedule = FALSE)) %dopar% {
  folds <- stratified_folds(breed_vec, 5L, seed = 500000L + b)
  rows <- list(); i <- 0L
  for (k in seq_len(5)) {
    tr_idx <- which(folds != k); te_idx <- which(folds == k)
    tr <- geno_mat[tr_idx, , drop = FALSE]; tb <- droplevels(breed_vec[tr_idx])
    te <- geno_mat[te_idx, , drop = FALSE]; eb <- breed_vec[te_idx]

    rk    <- rank_this_fold(tr, tb, snp_names, rng_seed = 600000L + 1000L * b + k)
    inner <- make_inner_setup(tr, tb, k = INNER_K,
                              seed = 700000L + 1000L * b + k, snp_names = snp_names)

    for (n in SIZES) for (m in METRICS) {
      fit <- tryCatch(
        nested_fit_predict(META, m, n, tr, tb, te, rk, inner,
                           seed = 42L + 1000L * b + 10L * k + match(m, METRICS)),
        error = function(e) list(pred = rep(NA_character_, length(eb))))
      i <- i + 1L
      rows[[i]] <- data.frame(split = b, fold = k, n_markers = n, metric = m,
                              acc = mean(fit$pred == as.character(eb), na.rm = TRUE))
    }

    # random panels: same folds, same inner tuning, only the panel is redrawn
    for (r in seq_len(N_RAND)) {
      set.seed(800000L + 10000L * b + 100L * k + r)
      perm_out <- snp_names[sample.int(length(snp_names))]
      outer_rank <- list(Random = perm_out)
      inner_folds <- stratified_folds(tb, INNER_K, seed = 700000L + 1000L * b + k)
      ids <- lapply(seq_len(INNER_K), function(f) which(inner_folds == f))
      ranks <- lapply(seq_len(INNER_K), function(f) {
        set.seed(800000L + 10000L * b + 100L * k + 1000L * r + f)
        list(Random = snp_names[sample.int(length(snp_names))])
      })
      rsetup <- list(folds = inner_folds, ids = ids, rankings = ranks, k = INNER_K)
      for (n in SIZES) {
        fit <- tryCatch(
          nested_fit_predict(META, "Random", n, tr, tb, te, outer_rank, rsetup,
                             seed = 42L + 500000L * b + 10L * k + r),
          error = function(e) list(pred = rep(NA_character_, length(eb))))
        i <- i + 1L
        rows[[i]] <- data.frame(split = b, fold = k, n_markers = n, metric = "Random",
                                acc = mean(fit$pred == as.character(eb), na.rm = TRUE))
      }
    }
  }
  do.call(rbind, rows)
}
stopImplicitCluster()

setDT(res)
fwrite(res, "repeated_outer_cv_detail.csv")
cat(sprintf("done: %d rows, %.1f min\n", nrow(res),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- per-split panel scores and paired differences -------------------------
res <- fread("repeated_outer_cv_detail.csv")
sm <- res[, .(acc = mean(acc, na.rm = TRUE)), by = .(split, n_markers, metric)]
wide <- dcast(sm, split + n_markers ~ metric, value.var = "acc")
setnames(wide, "Random", "rand")

cat("\nMean accuracy across splits (%):\n")
print(wide[, .(n_splits = .N,
               FST = round(100 * mean(FST), 2), PIC = round(100 * mean(PIC), 2),
               AED = round(100 * mean(AED), 2), Delta = round(100 * mean(Delta), 2),
               random = round(100 * mean(rand), 2)), by = n_markers])

cat("\nRanked minus random, paired across splits:\n")
for (n in SIZES) {
  sub <- wide[n_markers == n]
  for (m in METRICS) {
    d <- 100 * (sub[[m]] - sub$rand)
    tt <- t.test(d)
    cat(sprintf("  %d %-5s gain %+.2f pp  SD %.2f  95%% CI [%+.2f, %+.2f]  t(%d) %.1f  P = %.2g  positive in %d/%d splits\n",
                n, m, mean(d), sd(d), tt$conf.int[1], tt$conf.int[2],
                length(d) - 1, unname(tt$statistic), tt$p.value,
                sum(d > 0), length(d)))
  }
}
