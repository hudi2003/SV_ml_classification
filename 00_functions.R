stratified_folds <- function(y, k, seed = 42) {
  set.seed(seed)
  folds <- integer(length(y))
  for (lv in levels(y)) {
    idx <- which(y == lv)
    idx <- sample(idx)
    folds[idx] <- rep(1:k, length.out = length(idx))
  }
  folds
}
compute_fst_wc <- function(geno_sub, breed_sub) {
  breeds <- levels(droplevels(breed_sub))
  r <- length(breeds)
  L <- ncol(geno_sub)

  n_i   <- numeric(r)
  p_mat <- matrix(NA_real_, nrow = r, ncol = L)
  h_mat <- matrix(NA_real_, nrow = r, ncol = L)

  for (i in seq_len(r)) {
    idx <- which(breed_sub == breeds[i])
    n_i[i] <- length(idx)
    g <- geno_sub[idx, , drop = FALSE]
    p_mat[i, ] <- colMeans(g) / 2
    h_mat[i, ] <- colMeans(g == 1)
  }

  sum_n <- sum(n_i)
  nbar  <- mean(n_i)
  nc    <- (sum_n - sum(n_i^2) / sum_n) / (r - 1)

  w    <- n_i / sum_n
  pbar <- colSums(w * p_mat)
  hbar <- colSums(w * h_mat)
  s2   <- colSums((n_i * (p_mat - matrix(pbar, r, L, byrow = TRUE))^2)) / ((r - 1) * nbar)

  a <- (nbar / nc) * (s2 - (1 / (nbar - 1)) *
         (pbar * (1 - pbar) - ((r - 1) / r) * s2 - hbar / 4))
  b <- (nbar / (nbar - 1)) *
         (pbar * (1 - pbar) - ((r - 1) / r) * s2 - ((2 * nbar - 1) / (4 * nbar)) * hbar)
  c <- hbar / 2

  fst <- a / (a + b + c)
  fst[is.na(fst) | is.infinite(fst)] <- 0
  pmax(fst, 0)
}
compute_delta_aed_pic <- function(geno_sub, breed_sub, chunk_size = 50000) {
  breeds <- levels(droplevels(breed_sub))
  B <- length(breeds)
  L <- ncol(geno_sub)

  freq <- matrix(NA_real_, nrow = B, ncol = L)
  for (i in seq_len(B)) {
    idx <- which(breed_sub == breeds[i])
    freq[i, ] <- colMeans(geno_sub[idx, , drop = FALSE]) / 2
  }

  pop_var <- apply(freq, 2, function(x) mean((x - mean(x))^2))
  aed <- sqrt(2 * B / (B - 1) * pop_var)

  pair_idx <- which(lower.tri(matrix(0, B, B)), arr.ind = TRUE)
  delta <- numeric(L)
  starts <- seq(1, L, by = chunk_size)
  for (st in starts) {
    en <- min(st + chunk_size - 1, L)
    block <- freq[, st:en, drop = FALSE]
    diffs <- abs(block[pair_idx[, 1], , drop = FALSE] - block[pair_idx[, 2], , drop = FALSE])
    delta[st:en] <- colMeans(diffs)
  }

  overall_freq <- colMeans(geno_sub) / 2
  pic <- 2 * overall_freq * (1 - overall_freq) * (1 - overall_freq * (1 - overall_freq))

  list(Delta = delta, AED = aed, PIC = pic)
}
"%||%" <- function(a, b) if (is.null(a)) b else a

# Hyper-parameter grid used by the inner cross-validation.
# RF   : ntree x mtry            (mtry spans sqrt(p), p/5 and p/2)
# SVM  : cost in {0.01, 0.1, 1, 10}
# Other classifiers keep their pre-specified settings (grid = NULL).
hyper_grid <- function(method, p) {
  switch(method,
    "RF" = {
      ntrees <- c(300L, 500L, 1000L)
      mtries <- unique(as.integer(pmax(1L, pmin(p, c(floor(sqrt(p)),
                                                      floor(p / 5),
                                                      floor(p / 2))))))
      expand.grid(ntree = ntrees, mtry = mtries,
                  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    },
    "SVM_linear" = data.frame(cost = c(0.01, 0.1, 1, 10),
                              stringsAsFactors = FALSE),
    NULL
  )
}

fit_predict <- function(method, tr_x, tr_y, te_x, params = list()) {
  switch(method,
    "RF" = {
      ntree <- params$ntree %||% 500L
      mtry  <- params$mtry  %||% max(1L, floor(sqrt(ncol(tr_x))))
      fit <- randomForest::randomForest(x = tr_x, y = tr_y,
                                        ntree = ntree, mtry = mtry)
      as.character(predict(fit, te_x))
    },
    "SVM_linear" = {
      cost <- params$cost %||% 1
      fit <- e1071::svm(tr_x, tr_y, kernel = "linear", cost = cost,
                        scale = FALSE)
      as.character(predict(fit, te_x))
    },
    "SVM_RBF" = {
      fit <- e1071::svm(tr_x, tr_y, kernel = "radial", cost = 1, scale = FALSE)
      as.character(predict(fit, te_x))
    },
    "NB" = {
      fit <- e1071::naiveBayes(tr_x, tr_y)
      as.character(predict(fit, te_x))
    },
    "LDA" = {
      fit <- MASS::lda(tr_x, grouping = tr_y)
      as.character(predict(fit, te_x)$class)
    },
    "KNN" = {
      k <- params$k %||% 5L
      as.character(class::knn(tr_x, te_x, tr_y, k = k))
    },
    "DecisionTree" = {
      tr_df <- data.frame(y = tr_y, tr_x, check.names = FALSE)
      te_df <- data.frame(te_x, check.names = FALSE)
      fit <- rpart::rpart(y ~ ., data = tr_df, method = "class")
      as.character(predict(fit, te_df, type = "class"))
    },
    "XGBoost" = {
      y_tr_num <- as.integer(tr_y) - 1
      dtr <- xgboost::xgb.DMatrix(data = tr_x, label = y_tr_num)
      dte <- xgboost::xgb.DMatrix(data = te_x)
      fit <- xgboost::xgb.train(params = list(objective = "multi:softmax",
                                               num_class = nlevels(tr_y),
                                               max_depth = 6, eta = 0.3,
                                               nthread = 1),
                                 data = dtr, nrounds = 100, verbose = 0)
      pred_num <- predict(fit, dte)
      levels(tr_y)[pred_num + 1]
    },
    "MLP" = {
      y_ind <- nnet::class.ind(tr_y)
      fit <- nnet::nnet(x = tr_x, y = y_ind, size = 10, softmax = TRUE,
                         decay = 5e-4, maxit = 300, trace = FALSE, MaxNWts = 200000)
      pred_mat <- predict(fit, te_x)
      levels(tr_y)[apply(pred_mat, 1, which.max)]
    },
    "Lasso" = {
      # Stratified internal folds keep every breed represented in each fold;
      # the number of folds is capped by the rarest breed in the training set.
      k_lasso <- max(3L, min(5L, min(table(tr_y))))
      foldid  <- stratified_folds(tr_y, k_lasso, seed = 1234L)
      cv_fit <- glmnet::cv.glmnet(x = tr_x, y = tr_y, family = "multinomial",
                                  type.measure = "class", nfolds = k_lasso,
                                  foldid = foldid, alpha = 1)
      pred <- predict(cv_fit, newx = te_x, s = "lambda.min", type = "class")
      as.character(pred)
    },
    "LightGBM" = {
      y_tr_num <- as.integer(tr_y) - 1
      dtr <- lightgbm::lgb.Dataset(data = tr_x, label = y_tr_num)
      fit <- lightgbm::lgb.train(params = list(objective = "multiclass",
                                                num_class = nlevels(tr_y),
                                                learning_rate = 0.1,
                                                num_leaves = 31, verbose = -1),
                                  data = dtr, nrounds = 100, verbose = -1)
      pred_mat <- predict(fit, te_x, reshape = TRUE)
      pred_num <- apply(pred_mat, 1, which.max) - 1
      levels(tr_y)[pred_num + 1]
    },
    stop("Unknown classifier: ", method)
  )
}

# Build the inner (tuning) splits of one outer training fold.
# Rankings are recomputed inside every inner training split, so neither the
# outer test fold nor the inner validation samples contribute to the ranking
# used to fit/tune the classifier.
make_inner_setup <- function(train_geno, train_breed, k = 3, seed = 42,
                             snp_names = colnames(train_geno)) {
  inner_folds <- stratified_folds(train_breed, k, seed = seed)
  ids <- lapply(seq_len(k), function(f) which(inner_folds == f))
  rankings <- lapply(seq_len(k), function(f) {
    itr <- which(inner_folds != f)
    rank_this_fold(train_geno[itr, , drop = FALSE],
                   droplevels(train_breed[itr]), snp_names,
                   rng_seed = seed + f)
  })
  list(folds = inner_folds, ids = ids, rankings = rankings, k = k)
}

# Grid search on the inner folds only.  The outer test fold is never passed to
# this function.
inner_cv_select <- function(method, grid, setup, metric, n,
                            train_geno, train_breed, seed = 1L) {
  K <- setup$k
  acc_mat <- matrix(NA_real_, nrow = nrow(grid), ncol = K)
  for (g in seq_len(nrow(grid))) {
    params <- as.list(grid[g, , drop = FALSE])
    for (f in seq_len(K)) {
      va <- setup$ids[[f]]
      tr <- setdiff(seq_along(train_breed), va)
      snps <- head(setup$rankings[[f]][[metric]], n)
      tr_x <- train_geno[tr, snps, drop = FALSE]
      va_x <- train_geno[va, snps, drop = FALSE]
      set.seed(seed + 1000L * g + f)
      pred <- tryCatch(
        fit_predict(method, tr_x, train_breed[tr], va_x, params),
        error = function(e) rep(NA_character_, length(va))
      )
      acc_mat[g, f] <- mean(pred == as.character(train_breed[va]), na.rm = TRUE)
    }
  }
  mean_acc <- rowMeans(acc_mat, na.rm = TRUE)
  if (all(is.na(mean_acc))) {
    return(list(params = as.list(grid[1, , drop = FALSE]),
                inner_acc = NA_real_, inner_sd = NA_real_,
                inner_acc_mat = acc_mat))
  }
  best <- which.max(mean_acc)          # ties -> first (smallest) configuration
  list(params = as.list(grid[best, , drop = FALSE]),
       inner_acc = mean_acc[best],
       inner_sd = if (K > 1) stats::sd(acc_mat[best, ], na.rm = TRUE) else NA_real_,
       inner_acc_mat = acc_mat)
}

# Full nested prediction for one (metric, panel size, classifier) cell:
# tune on the inner folds of the outer training set, refit on the whole outer
# training set, then predict the untouched outer test set.
nested_fit_predict <- function(method, metric, n, train_geno, train_breed,
                               test_geno, outer_ranking, setup, seed = 1L) {
  snps <- head(outer_ranking[[metric]], n)
  grid <- hyper_grid(method, length(snps))

  if (is.null(grid)) {
    sel <- list(params = list(), inner_acc = NA_real_, inner_sd = NA_real_)
  } else {
    sel <- inner_cv_select(method, grid, setup, metric, n,
                           train_geno, train_breed, seed = seed)
  }

  set.seed(seed + 7L)
  pred <- fit_predict(method,
                      train_geno[, snps, drop = FALSE], train_breed,
                      test_geno[, snps, drop = FALSE], sel$params)

  list(pred = pred, params = sel$params,
       inner_acc = sel$inner_acc, inner_sd = sel$inner_sd,
       searched = !is.null(grid),
       grid_size = if (is.null(grid)) 0L else nrow(grid))
}

params_to_string <- function(params) {
  if (length(params) == 0) return("fixed")
  paste(names(params), unlist(params), sep = "=", collapse = ";")
}

rank_this_fold <- function(train_geno, train_breed, snp_names, add_rf_mda = FALSE,
                           rng_seed = NULL) {
  fst_train <- compute_fst_wc(train_geno, train_breed)
  dap       <- compute_delta_aed_pic(train_geno, train_breed)

  if (!is.null(rng_seed)) set.seed(rng_seed)
  out <- list(
    FST    = snp_names[order(-fst_train)],
    Delta  = snp_names[order(-dap$Delta)],
    AED    = snp_names[order(-dap$AED)],
    PIC    = snp_names[order(-dap$PIC)],
    Random = snp_names[sample.int(length(snp_names))]
  )

  if (add_rf_mda) {
    rf_imp_model <- randomForest::randomForest(x = train_geno, y = train_breed,
                                                ntree = 500, importance = TRUE)
    mda <- randomForest::importance(rf_imp_model, type = 1)[, "MeanDecreaseAccuracy"]
    out$RF_MDA <- colnames(train_geno)[order(-mda)]
  }
  out
}
