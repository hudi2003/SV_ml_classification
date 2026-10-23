# ---------------------------------------------------------------------------
# Leakage / sanity checks for the nested cross-validation benchmark.
# Run BEFORE the full benchmark.  All checks must print PASS.
# ---------------------------------------------------------------------------
library(data.table)
source("00_functions.R")

geno_mat  <- readRDS("geno_mat.rds")
breed_vec <- as.factor(readRDS("breed_vec.rds"))
snp_names <- readRDS("snp_names.rds")

folds <- stratified_folds(breed_vec, 5, seed = 42)
k <- 1
train_idx <- which(folds != k)
test_idx  <- which(folds == k)

train_geno  <- geno_mat[train_idx, , drop = FALSE]
train_breed <- droplevels(breed_vec[train_idx])
test_geno   <- geno_mat[test_idx, , drop = FALSE]
test_breed  <- breed_vec[test_idx]

ok <- TRUE
check <- function(label, cond) {
  cat(sprintf("[%s] %s\n", ifelse(isTRUE(cond), "PASS", "FAIL"), label))
  if (!isTRUE(cond)) ok <<- FALSE
}

cat("=== 1. Split disjointness ===\n")
check("outer training and test indices are disjoint", length(intersect(train_idx, test_idx)) == 0)
check("outer folds cover all samples", length(union(train_idx, test_idx)) == length(breed_vec))

ranking <- rank_this_fold(train_geno, train_breed, snp_names)
setup   <- make_inner_setup(train_geno, train_breed, k = 5, seed = 2024L + k,
                            snp_names = snp_names)

check("inner validation ids index only outer-training rows",
      all(unlist(setup$ids) <= nrow(train_geno)))
check("inner folds are disjoint and cover the outer training set", {
  all_ids <- unlist(setup$ids)
  length(all_ids) == nrow(train_geno) && !anyDuplicated(all_ids)
})
for (f in 1:5) {
  check(sprintf("inner fold %d excludes its own validation rows when ranking", f),
        length(intersect(which(setup$folds != f), setup$ids[[f]])) == 0)
}

cat("\n=== 2. Inner rankings are recomputed, not inherited ===\n")
same_as_outer <- vapply(1:5, function(f)
  identical(head(setup$rankings[[f]][["FST"]], 100), head(ranking[["FST"]], 100)),
  logical(1))
check("inner rankings differ from the outer-training ranking (top-100)",
      !all(same_as_outer))

cat("\n=== 3. Hyper-parameter selection is invariant to the outer test fold ===\n")
set.seed(99)
noise_geno <- matrix(sample(0:2, nrow(test_geno) * ncol(test_geno), replace = TRUE),
                     nrow = nrow(test_geno), dimnames = list(NULL, snp_names))
shuffled_y <- sample(test_breed)

for (cl in c("RF", "SVM_linear")) {
  for (n in c(20, 100)) {
    seed_i <- 42L + 10L * n + match(cl, c("RF", "SVM_linear"))
    a <- nested_fit_predict(cl, "FST", n, train_geno, train_breed, test_geno,
                            ranking, setup, seed = seed_i)
    b <- nested_fit_predict(cl, "FST", n, train_geno, train_breed, noise_geno,
                            ranking, setup, seed = seed_i)
    check(sprintf("%s n=%d: identical selected params with scrambled test data", cl, n),
          identical(a$params, b$params) && isTRUE(all.equal(a$inner_acc, b$inner_acc)))
  }
}

cat("\n=== 4. Negative control: permuted outer test labels ===\n")
for (cl in c("RF", "SVM_linear")) {
  accs <- sapply(1:5, function(kk) {
    ti <- which(folds != kk)
    tg <- geno_mat[ti, , drop = FALSE]
    tb <- droplevels(breed_vec[ti])
    te <- geno_mat[which(folds == kk), , drop = FALSE]
    rk <- rank_this_fold(tg, tb, snp_names)
    st <- make_inner_setup(tg, tb, k = 5, seed = 2024L + kk, snp_names = snp_names)
    f  <- nested_fit_predict(cl, "FST", 100, tg, tb, te, rk, st,
                             seed = 42L + 10L * 100 + match(cl, c("RF", "SVM_linear")))
    mean(f$pred == as.character(sample(breed_vec[which(folds == kk)])))
  })
  cat(sprintf("      %s permuted-label accuracy per fold: %s (mean %.3f)\n",
              cl, paste(sprintf("%.3f", accs), collapse = ", "), mean(accs)))
  check(sprintf("%s permuted-label accuracy <= 0.15 (chance = %.3f)", cl, 1 / nlevels(breed_vec)),
        mean(accs) <= 0.15)
}

cat("\n=== 5. Inner grid definitions ===\n")
g_rf  <- hyper_grid("RF", 10)
g_rf2 <- hyper_grid("RF", 3000)
g_svm <- hyper_grid("SVM_linear", 100)
check("SVM linear grid = {0.01, 0.1, 1, 10}", identical(g_svm$cost, c(0.01, 0.1, 1, 10)))
check("RF grid searches ntree and mtry", all(c("ntree", "mtry") %in% names(g_rf)) &&
        length(unique(g_rf$ntree)) == 3)
check("RF mtry never exceeds p at p=10", all(g_rf$mtry <= 10))
check("RF mtry never exceeds p at p=3000", all(g_rf2$mtry <= 3000))

cat(sprintf("\nOVERALL: %s\n", ifelse(ok, "ALL CHECKS PASSED", "SOME CHECKS FAILED")))
if (!ok) quit(status = 1)
