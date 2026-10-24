# -*- coding: utf-8 -*-
library(data.table); library(ggplot2); library(dplyr)
source("00_functions.R")

geno_mat  <- readRDS("geno_mat.rds")
breed_vec <- as.factor(readRDS("breed_vec.rds"))
snp_names <- readRDS("snp_names.rds")
breeds <- sort(levels(breed_vec))

CONF_METRIC <- "FST"; CONF_CLASSIFIER <- "SVM_linear"; N <- 600
folds <- stratified_folds(breed_vec, 5, seed = 42)

all_true <- c(); all_pred <- c()
for (k in 1:5) {
  tr <- which(folds != k)
  tg <- geno_mat[tr, , drop = FALSE]; tb <- droplevels(breed_vec[tr])
  eg <- geno_mat[which(folds == k), , drop = FALSE]; eb <- breed_vec[which(folds == k)]
  rk <- rank_this_fold(tg, tb, snp_names)
  st <- make_inner_setup(tg, tb, k = 5, seed = 2024L + k, snp_names = snp_names)
  fit <- nested_fit_predict(CONF_CLASSIFIER, CONF_METRIC, N, tg, tb, eg, rk, st,
                            seed = 42L + 10000L * k + 10L * N)
  cat(sprintf("fold %d: %s (inner %.3f)\n", k, params_to_string(fit$params), fit$inner_acc))
  all_true <- c(all_true, as.character(eb)); all_pred <- c(all_pred, as.character(fit$pred))
}
acc <- mean(all_true == all_pred)
cat(sprintf("Overall accuracy: %.4f\n", acc))

conf_mat  <- table(True = factor(all_true, levels = breeds), Pred = factor(all_pred, levels = breeds))
conf_prop <- sweep(conf_mat, 1, rowSums(conf_mat), "/")
conf_df <- as.data.frame(conf_prop); colnames(conf_df) <- c("True", "Pred", "Proportion")
conf_df$label <- sprintf("%.2f", conf_df$Proportion)
conf_df$txt_color <- ifelse(conf_df$Proportion > 0.5, "white", "gray20")

p <- ggplot(conf_df, aes(x = Pred, y = True, fill = Proportion)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = label, color = txt_color), size = 3.0, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_gradient2(low = "white", mid = "#85B7EB", high = "#0C447C",
                       midpoint = 0.5, limits = c(0, 1), name = "Proportion") +
  labs(x = "Predicted breed", y = "True breed",
       title = sprintf("%d markers (%s + %s)", N, CONF_METRIC, CONF_CLASSIFIER),
       subtitle = sprintf("Accuracy: %.1f%%", acc * 100)) +
  theme_classic(base_size = 23) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 21, face = "bold"),
        plot.subtitle = element_text(size = 17),
        legend.key.size = unit(1, "cm"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12))

ggsave("Figure4b_confusion_600SV_FST_SVM.pdf", p, width = 9, height = 8, device = cairo_pdf)
ggsave("Figure4b_confusion_600SV_FST_SVM.png", p, width = 9, height = 8, dpi = 300)
cat("saved Figure4b (FST 600, SVM_linear)\n")
