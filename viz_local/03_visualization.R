library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)

source("00_functions.R")

classifier_colors <- c(
  RF             = "#3B8BD4",
  SVM_linear     = "#D85A30",
  SVM_RBF        = "#E8A33D",
  NB             = "#1D9E75",
  LDA            = "#7F77DD",
  KNN            = "#C25CA8",
  DecisionTree   = "#8C6D3F",
  XGBoost        = "#2C2C54",
  MLP            = "#54A0FF",
  Lasso          = "#10AC84",
  LightGBM       = "#EE5253"
)

res <- fread("nested_cv_full_summary.csv")
x_breaks <- sort(unique(res$n_markers))

ACC_THRESHOLDS  <- c(0.90, 0.95, 0.99)
PRIMARY_THRESHOLD <- 0.95

economical_list <- lapply(ACC_THRESHOLDS, function(thr) {
  res %>%
    filter(mean_acc >= thr) %>%
    group_by(metric, classifier) %>%
    slice_min(n_markers, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(threshold = thr)
})
economical_tab <- bind_rows(economical_list)
fwrite(economical_tab, "economical_panel_summary.csv")
cat("Economical panel summary saved: economical_panel_summary.csv\n")

primary_tab <- economical_tab %>% filter(threshold == PRIMARY_THRESHOLD)

if (nrow(primary_tab) == 0) {
  cat(sprintf("WARNING: no combination reached %.0f%% accuracy within tested panel sizes.\n",
              PRIMARY_THRESHOLD * 100))
  cat("Falling back to the globally highest-accuracy combination instead.\n")
  best_economical <- res %>% slice_max(mean_acc, n = 1, with_ties = FALSE)
} else {
  best_economical <- primary_tab %>% slice_min(n_markers, n = 1, with_ties = FALSE)
}

ECO_METRIC     <- best_economical$metric[1]
ECO_CLASSIFIER <- best_economical$classifier[1]
ECO_N          <- best_economical$n_markers[1]
cat(sprintf("\n*** Most economical combination (>= %.0f%% accuracy): metric=%s, classifier=%s, n_markers=%d, acc=%.4f ***\n",
            PRIMARY_THRESHOLD * 100, ECO_METRIC, ECO_CLASSIFIER, ECO_N, best_economical$mean_acc[1]))

if (nrow(economical_tab) > 0) {
  economical_tab$combo <- paste(economical_tab$metric, economical_tab$classifier, sep = " + ")
  economical_tab$threshold_label <- paste0(economical_tab$threshold * 100, "%")

  p_eco <- ggplot(economical_tab, aes(x = reorder(combo, n_markers), y = n_markers,
                                       fill = metric)) +
    geom_col() +
    geom_text(aes(label = n_markers), hjust = -0.1, size = 2.8) +
    coord_flip() +
    facet_wrap(~ threshold_label, ncol = 1, scales = "free_y") +
    scale_y_log10() +
    labs(title = "Minimum number of markers needed to reach each accuracy threshold",
         x = "Ranking method + Classifier", y = "Minimum number of markers (log scale)",
         fill = "Ranking method") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          strip.text = element_text(face = "bold"))

  ggsave("fig_economical_selection.pdf", p_eco,
         width = 9, height = 3 + 1.5 * length(ACC_THRESHOLDS))
  ggsave("fig_economical_selection.png", p_eco,
         width = 9, height = 3 + 1.5 * length(ACC_THRESHOLDS), dpi = 300)
  cat("Economical selection plot saved: fig_economical_selection.png/pdf\n")
}
plot_accuracy_curve <- function(metric_name) {
  df <- res %>% filter(metric == metric_name)

  p <- ggplot(df, aes(x = n_markers, y = mean_acc,
                       color = classifier, fill = classifier)) +
    geom_ribbon(aes(ymin = mean_acc - sd_acc, ymax = mean_acc + sd_acc),
                alpha = 0.08, color = NA) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 1.8) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_hline(yintercept = 0.99, linetype = "dotted", color = "gray50", linewidth = 0.5) +
    annotate("text", x = min(x_breaks), y = 0.958, label = "95%",
             size = 3, color = "gray50", hjust = 0) +
    annotate("text", x = min(x_breaks), y = 0.997, label = "99%",
             size = 3, color = "gray50", hjust = 0) +
    scale_x_log10(breaks = x_breaks, labels = x_breaks) +
    scale_y_continuous(limits = c(0, 1.05), labels = scales::percent_format(accuracy = 1)) +
    scale_color_manual(values = classifier_colors) +
    scale_fill_manual(values = classifier_colors) +
    labs(title = paste0("Classifier accuracy under ", metric_name, " ranking"),
         x = "Number of markers (log scale)", y = "Nested CV accuracy",
         color = "Classifier", fill = "Classifier") +
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 8),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

  ggsave(paste0("fig_accuracy_", metric_name, ".pdf"), p, width = 9, height = 5.5)
  ggsave(paste0("fig_accuracy_", metric_name, ".png"), p, width = 9, height = 5.5, dpi = 300)
}

for (m in unique(res$metric)) plot_accuracy_curve(m)
cat("Accuracy curve figures saved (one per metric: fig_accuracy_*.png/pdf)\n")

p_time <- ggplot(res, aes(x = classifier, y = mean_time, fill = classifier)) +
  geom_boxplot(outlier.size = 0.6) +
  scale_fill_manual(values = classifier_colors) +
  scale_y_log10() +
  labs(title = "Average train+predict runtime per classifier (all panels & metrics pooled)",
       x = NULL, y = "Runtime (seconds, log scale)") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 11))

ggsave("fig_runtime_comparison.pdf", p_time, width = 8, height = 5)
ggsave("fig_runtime_comparison.png", p_time, width = 8, height = 5, dpi = 300)
cat("Runtime comparison figure saved: fig_runtime_comparison.png/pdf\n")

best_tab <- res %>%
  group_by(metric, n_markers) %>%
  slice_max(mean_acc, n = 1, with_ties = FALSE) %>%
  ungroup()

p_heat <- ggplot(best_tab, aes(x = factor(n_markers), y = metric, fill = mean_acc)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(classifier, "\n", sprintf("%.3f", mean_acc))),
            size = 2.6, color = "white") +
  scale_fill_gradient(low = "#9ecae1", high = "#08306b", limits = c(0, 1)) +
  labs(title = "Best-performing classifier at each marker count (rigorous nested CV)",
       x = "Number of markers", y = "Ranking method", fill = "Accuracy") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 11))

ggsave("fig_best_classifier_heatmap.pdf", p_heat, width = 10, height = 4)
ggsave("fig_best_classifier_heatmap.png", p_heat, width = 10, height = 4, dpi = 300)
cat("Best-classifier heatmap saved: fig_best_classifier_heatmap.png/pdf\n")

geno_mat  <- readRDS("geno_mat.rds")
breed_vec <- as.factor(readRDS("breed_vec.rds"))
snp_names <- readRDS("snp_names.rds")
breeds <- sort(levels(breed_vec))

avail_sizes <- sort(unique(res$n_markers[res$metric == ECO_METRIC & res$classifier == ECO_CLASSIFIER]))
idx_eco <- which(avail_sizes == ECO_N)
CONF_PANEL_SIZES <- avail_sizes[idx_eco:min(length(avail_sizes), idx_eco + 2)]
CONF_PANEL_SIZES <- c(100, 200, 300, 400, 500, 600, 1000)
CONF_METRIC     <- ECO_METRIC
CONF_CLASSIFIER <- ECO_CLASSIFIER

cat(sprintf("\nDrawing confusion matrices for metric=%s, classifier=%s, panel sizes=%s\n",
            CONF_METRIC, CONF_CLASSIFIER, paste(CONF_PANEL_SIZES, collapse = ",")))

n_folds_conf <- 5
folds <- stratified_folds(breed_vec, n_folds_conf)

plot_list <- list()

# Pre-compute the outer ranking and the inner tuning setup once per outer fold.
conf_cache <- lapply(1:n_folds_conf, function(k) {
  train_idx <- which(folds != k)
  train_geno  <- geno_mat[train_idx, , drop = FALSE]
  train_breed <- droplevels(breed_vec[train_idx])
  list(train_geno  = train_geno,
       train_breed = train_breed,
       test_geno   = geno_mat[which(folds == k), , drop = FALSE],
       test_breed  = breed_vec[which(folds == k)],
       ranking     = rank_this_fold(train_geno, train_breed, snp_names),
       inner       = make_inner_setup(train_geno, train_breed, k = 5,
                                        seed = 2024L + k, snp_names = snp_names))
})

for (n in CONF_PANEL_SIZES) {
  cat(sprintf("\n=== Panel size: %d ===\n", n))

  all_true <- c(); all_pred <- c()
  for (k in 1:n_folds_conf) {
    cc <- conf_cache[[k]]
    train_geno  <- cc$train_geno
    train_breed <- cc$train_breed
    test_geno   <- cc$test_geno
    test_breed  <- cc$test_breed

    fit <- nested_fit_predict(CONF_CLASSIFIER, CONF_METRIC, n,
                              train_geno, train_breed, test_geno,
                              cc$ranking, cc$inner, seed = 42L + 10000L * k + 10L * n)
    pred <- fit$pred
    cat(sprintf("  fold %d: %s (inner acc %.3f)\n",
                k, params_to_string(fit$params), fit$inner_acc))

    all_true <- c(all_true, as.character(test_breed))
    all_pred <- c(all_pred, as.character(pred))
  }

  overall_acc <- mean(all_true == all_pred)
  cat(sprintf("Overall accuracy: %.4f\n", overall_acc))

  conf_mat  <- table(True = factor(all_true, levels = breeds),
                      Pred = factor(all_pred, levels = breeds))
  conf_prop <- sweep(conf_mat, 1, rowSums(conf_mat), "/")
  conf_df   <- as.data.frame(conf_prop)
  colnames(conf_df) <- c("True", "Pred", "Proportion")

  conf_df$label     <- sprintf("%.2f", conf_df$Proportion)
  conf_df$txt_color <- ifelse(conf_df$Proportion > 0.5, "white", "gray20")

  p <- ggplot(conf_df, aes(x = Pred, y = True, fill = Proportion)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = label, color = txt_color), size = 2.6, show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_gradient2(low = "white", mid = "#85B7EB", high = "#0C447C",
                          midpoint = 0.5, limits = c(0, 1),
                          name = "Proportion") +
    labs(x = "Predicted breed", y = "True breed",
         title = sprintf("%d markers (%s + %s)", n, CONF_METRIC, CONF_CLASSIFIER),
         subtitle = sprintf("Accuracy: %.1f%%", overall_acc * 100)) +
    theme_classic(base_size = 20) +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y     = element_text(size = 14),
          plot.title      = element_text(size = 18, face = "bold"),
          plot.subtitle   = element_text(size = 14),
          legend.key.size = unit(1, "cm"),
          legend.title    = element_text(size = 16),
          legend.text     = element_text(size = 10))

  plot_list[[as.character(n)]] <- p
}

for (n in CONF_PANEL_SIZES) {
  p <- plot_list[[as.character(n)]]
  ggsave(sprintf("fig_confusion_%dSNPs.pdf", n), p, width = 9, height = 8)
  ggsave(sprintf("fig_confusion_%dSNPs.png", n), p, width = 9, height = 8, dpi = 300)
  cat(sprintf("Saved: fig_confusion_%dSNPs.pdf / .png\n", n))
}

cat(format(Sys.time(), "%H:%M:%S"), "Checking ranking stability across folds...\n")

stability_check_n <- max(CONF_PANEL_SIZES)
fold_rankings <- list()
for (k in 1:n_folds_conf) {
  train_idx <- which(folds != k)
  train_geno  <- geno_mat[train_idx, , drop = FALSE]
  train_breed <- droplevels(breed_vec[train_idx])
  ranking_this_fold <- rank_this_fold(train_geno, train_breed, snp_names)
  fold_rankings[[k]] <- head(ranking_this_fold[[CONF_METRIC]], stability_check_n)
}

overlap_mat <- matrix(NA, n_folds_conf, n_folds_conf)
for (i in 1:n_folds_conf) {
  for (j in 1:n_folds_conf) {
    overlap_mat[i, j] <- length(intersect(fold_rankings[[i]], fold_rankings[[j]])) / stability_check_n
  }
}
overlap_df <- as.data.frame(overlap_mat)
colnames(overlap_df) <- paste0("Fold", 1:n_folds_conf)
overlap_df$Fold <- paste0("Fold", 1:n_folds_conf)
overlap_long <- pivot_longer(overlap_df, -Fold, names_to = "Fold2", values_to = "Overlap")

p_stab <- ggplot(overlap_long, aes(x = Fold2, y = Fold, fill = Overlap)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Overlap)), size = 3.2, color = "white") +
  scale_fill_gradient(low = "#fee5d9", high = "#a50f15", limits = c(0, 1)) +
  labs(title = sprintf("Overlap of top-%d %s markers across folds", stability_check_n, CONF_METRIC),
       x = NULL, y = NULL, fill = "Overlap proportion") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 11))

ggsave("fig_ranking_stability.pdf", p_stab, width = 6, height = 5)
ggsave("fig_ranking_stability.png", p_stab, width = 6, height = 5, dpi = 300)
cat("Ranking stability figure saved: fig_ranking_stability.png/pdf\n")

cat(format(Sys.time(), "%H:%M:%S"), "All visualizations complete.\n")
