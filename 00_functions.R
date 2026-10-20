#!/usr/bin/env Rscript
# ================================================================
# 脚本：00_functions.R
# 功能：公共函数库，被 02_rigorous_nested_benchmark.R 和
#       03_visualization.R 共同 source()。
#       新增：概率级 Stacking 集成方法（候选池扩展 + glmnet正则化元学习器）。
# 不要单独运行这个脚本，它只提供函数定义。
# ================================================================

# ---- 分层分折 ----
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

# ---- FST：Weir & Cockerham (1984) ----
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

# ---- Delta / AED / PIC ----
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

# ---- 单个分类器训练/预测封装（硬标签输出） ----
fit_predict <- function(method, tr_x, tr_y, te_x) {
  switch(method,
    "RF" = {
      fit <- randomForest::randomForest(x = tr_x, y = tr_y, ntree = 500)
      as.character(predict(fit, te_x))
    },
    "SVM_linear" = {
      fit <- e1071::svm(tr_x, tr_y, kernel = "linear", scale = FALSE)
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
      as.character(class::knn(tr_x, te_x, tr_y, k = 5))
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
      fit <- xgboost::xgb.train(data = dtr, nrounds = 100, max_depth = 6, eta = 0.3,
                                 objective = "multi:softmax",
                                 num_class = nlevels(tr_y), verbose = 0)
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
      cv_fit <- glmnet::cv.glmnet(x = tr_x, y = tr_y, family = "multinomial",
                                  type.measure = "class", nfolds = 5, alpha = 1)
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
    stop("未知分类器: ", method)
  )
}

# ================================================================
# ---- 概率级 Stacking 相关函数 ----
# ================================================================

# ---- 输出概率矩阵版本的基学习器封装：返回 n_test x n_class 矩阵，列名=品种 ----

# ---- 概率级 Stacking：内层K折生成OOF概率 -> 拼接为元特征 -> glmnet正则化元学习器自动筛选权重 ----

# ---- 给定一折训练集，重新算四个指标的排序，返回 named list ----
rank_this_fold <- function(train_geno, train_breed, snp_names, add_rf_mda = FALSE) {
  fst_train <- compute_fst_wc(train_geno, train_breed)
  dap       <- compute_delta_aed_pic(train_geno, train_breed)

  out <- list(
    FST   = snp_names[order(-fst_train)],
    Delta = snp_names[order(-dap$Delta)],
    AED   = snp_names[order(-dap$AED)],
    PIC   = snp_names[order(-dap$PIC)]
  )

  if (add_rf_mda) {
    rf_imp_model <- randomForest::randomForest(x = train_geno, y = train_breed,
                                                ntree = 500, importance = TRUE)
    mda <- randomForest::importance(rf_imp_model, type = 1)[, "MeanDecreaseAccuracy"]
    out$RF_MDA <- colnames(train_geno)[order(-mda)]
  }
  out
}