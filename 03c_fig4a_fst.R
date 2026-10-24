library(data.table); library(ggplot2); library(scales)

res <- fread("nested_cv_full_summary.csv")
x_breaks <- sort(unique(res$n_markers))
classifier_colors <- c(
  RF = "#3B8BD4", SVM_linear = "#D85A30", SVM_RBF = "#E8A33D", NB = "#1D9E75",
  LDA = "#7F77DD", KNN = "#C25CA8", DecisionTree = "#8C6D3F", XGBoost = "#2C2C54",
  MLP = "#54A0FF", Lasso = "#10AC84", LightGBM = "#EE5253")

df <- res[metric == "FST"]
p <- ggplot(df, aes(x = n_markers, y = mean_acc, color = classifier, fill = classifier)) +
  geom_ribbon(aes(ymin = mean_acc - sd_acc, ymax = mean_acc + sd_acc), alpha = 0.08, color = NA) +
  geom_line(linewidth = 1.0) + geom_point(size = 1.8) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_hline(yintercept = 0.99, linetype = "dotted", color = "gray50", linewidth = 0.5) +
  annotate("text", x = min(x_breaks), y = 0.958, label = "95%", size = 3, color = "gray50", hjust = 0) +
  annotate("text", x = min(x_breaks), y = 0.997, label = "99%", size = 3, color = "gray50", hjust = 0) +
  scale_x_log10(breaks = x_breaks, labels = x_breaks) +
  scale_y_continuous(limits = c(0, 1.05), labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = classifier_colors) + scale_fill_manual(values = classifier_colors) +
  labs(title = "Classifier accuracy under FST ranking", x = "Number of markers (log scale)",
       y = "Nested CV accuracy", color = "Classifier", fill = "Classifier") +
  theme_classic(base_size = 19) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 11),
        legend.position = "right", plot.title = element_text(hjust = 0.5, size = 19, face = "bold"))
ggsave("Figure4a_accuracy_FST_nestedCV.pdf", p, width = 9, height = 5.5, device = cairo_pdf)
ggsave("Figure4a_accuracy_FST_nestedCV.png", p, width = 9, height = 5.5, dpi = 300)
cat("4a done\n")
