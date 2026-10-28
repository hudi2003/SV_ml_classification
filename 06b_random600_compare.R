# ---------------------------------------------------------------------------
# 06b – statistics and figure for the 100 random 600-SV panels vs the four
# ranking methods at the same panel size (SVM linear, nested 5x5 CV).
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2); library(tidyr); library(patchwork)
})
source("00_functions.R")

res <- fread("random600_detail.csv")
det <- fread("nested_cv_full_detail.csv")
N_MARK <- 600L; META <- "SVM_linear"; N_REP <- max(res$rep)
stopifnot(all(res$status == "ok"))

ref <- det[det$n_markers == N_MARK & det$classifier == META &
             det$metric %in% c("FST", "Delta", "AED", "PIC") & det$status == "ok"]
ref <- ref[order(match(ref$metric, c("FST", "AED", "Delta", "PIC")), ref$fold), ]

rep_mean <- res %>% group_by(rep) %>%
  summarise(mean_acc = mean(acc), sd_acc = sd(acc), .groups = "drop")
fold_random <- res %>% group_by(fold) %>%
  summarise(r_mean = mean(acc), r_sd = sd(acc), r_lo = quantile(acc, .025),
            r_hi = quantile(acc, .975), .groups = "drop")

# ---- paired comparison on the shared outer folds ---------------------------
pairs <- ref %>% left_join(fold_random, by = "fold") %>%
  mutate(diff_pp = 100 * (acc - r_mean))
stat <- pairs %>% group_by(metric) %>%
  summarise(ranked = 100 * mean(acc), ranked_sd = 100 * sd(acc),
            random = 100 * mean(r_mean), random_sd = 100 * sd(r_sd),
            mean_diff_pp = mean(diff_pp),
            p_paired_t = tryCatch(t.test(diff_pp)$p.value, error = function(e) NA),
            p_wilcoxon = tryCatch(suppressWarnings(wilcox.test(diff_pp)$p.value),
                                  error = function(e) NA),
            .groups = "drop") %>%
  mutate(p_empirical = sapply(ranked, function(m) mean(100 * rep_mean$mean_acc >= m)),
         draws_ge = round(p_empirical * N_REP),
         pctile = 100 * sapply(ranked, function(m) mean(100 * rep_mean$mean_acc < m))) %>%
  rename(diff_pp = mean_diff_pp) %>%
  select(metric, ranked, ranked_sd, random, diff_pp, draws_ge, p_empirical,
         pctile, p_paired_t, p_wilcoxon)
fwrite(stat, "random600_stats.csv")

# ---- marker overlap: random panels vs FST top-600 (same outer fold) --------
geno_mat  <- readRDS("geno_mat.rds"); breed_vec <- as.factor(readRDS("breed_vec.rds"))
snp_names <- readRDS("snp_names.rds")
bnd <- grepl("\\[|\\]", snp_names)
geno_mat <- geno_mat[, !bnd, drop = FALSE]; snp_names <- snp_names[!bnd]
folds <- stratified_folds(breed_vec, 5L, seed = 42L)
ov <- data.frame()
for (k in 1:5) {
  tr <- which(folds != k)
  fst_k <- compute_fst_wc(geno_mat[tr, , drop = FALSE],
                          droplevels(breed_vec[tr]))
  top_fst <- snp_names[order(-fst_k)][1:N_MARK]
  for (r in 1:N_REP) {
    set.seed(300000L + 1000L * r + k)
    perm <- snp_names[sample.int(length(snp_names))][1:N_MARK]
    ov <- rbind(ov, data.frame(fold = k, rep = r,
                               overlap = length(intersect(perm, top_fst))))
  }
}
ov_sum <- sprintf(paste0("Overlap of a random 600-SV panel with the FST top-600 panel of ",
                         "the same outer training fold: %.1f +/- %.1f SVs (%.1f%%), ",
                         "chance expectation %.1f (%.1f%%); max %d"),
                  mean(ov$overlap), sd(ov$overlap),
                  100 * mean(ov$overlap) / N_MARK,
                  N_MARK^2 / length(snp_names),
                  100 * N_MARK / length(snp_names), max(ov$overlap))

# ---- figure ---------------------------------------------------------------
lvl <- c("FST", "PIC", "AED", "Delta")
cols <- c(FST = "#B2182B", PIC = "#EF8A62", AED = "#2166AC", Delta = "#67A9CF")
grid_cols <- c(FST = "#B2182B", AED = "#2166AC", Delta = "#67A9CF", PIC = "#EF8A62")

pA <- ggplot(rep_mean, aes(x = 100 * mean_acc)) +
  geom_histogram(aes(y = after_stat(density)), bins = 18,
                 fill = "grey80", colour = "white") +
  geom_density(linewidth = 0.7, colour = "grey20") +
  geom_vline(data = stat, aes(xintercept = ranked, colour = metric),
             linetype = "dashed", linewidth = 0.8) +
  geom_text(data = stat, aes(x = ranked, y = Inf, label = metric, colour = metric),
            angle = 90, vjust = -0.3, hjust = 1.05, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = grid_cols, breaks = lvl) +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     limits = c(floor(min(100 * rep_mean$mean_acc)) - 0.2, 97)) +
  labs(x = "Accuracy of a random 600-SV panel (mean of 5 outer folds)",
       y = "Density", colour = "Ranking",
       title = sprintf("A  %d random 600-SV panels vs. ranking methods", N_REP)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

pB <- ggplot() +
  geom_boxplot(data = res, aes(x = factor(fold), y = 100 * acc),
               width = 0.62, outlier.size = 0.7, fill = "grey90",
               colour = "grey35") +
  geom_point(data = ref, aes(x = factor(fold), y = 100 * acc, colour = metric),
             position = position_dodge(width = 0.55), size = 2.1) +
  scale_colour_manual(values = cols, breaks = lvl) +
  labs(x = "Outer fold", y = "Accuracy", colour = "Ranking",
       title = "B  Paired by outer fold") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p <- pA + pB + plot_layout(widths = c(1.35, 1))
ggsave("FigureS_random600_panels.png", p, width = 10, height = 4.2, dpi = 300)
ggsave("FigureS_random600_panels.pdf", p, width = 10, height = 4.2)

# ---- console report -------------------------------------------------------
txt <- c(
  sprintf("Random 600-SV panels: %d draws, nested 5x5, SVM linear", N_REP),
  sprintf("  draw mean: %.2f%% +/- %.2f%%  (range %.2f-%.2f%%, 2.5-97.5%%: %.2f-%.2f%%)",
          100 * mean(rep_mean$mean_acc), 100 * sd(rep_mean$mean_acc),
          100 * min(rep_mean$mean_acc), 100 * max(rep_mean$mean_acc),
          100 * quantile(rep_mean$mean_acc, .025), 100 * quantile(rep_mean$mean_acc, .975)),
  "",
  "Ranking method vs random null:",
  sprintf("  %-5s ranked %.2f%% (SD %.2f) | random %.2f%% | diff %+.2f pp | draws >= ranked %d/%d (p=%.2f) | paired t p=%.3f",
          stat$metric, stat$ranked, stat$ranked_sd, stat$random, stat$diff_pp,
          stat$draws_ge, N_REP, stat$p_empirical, stat$p_paired_t),
  "", ov_sum)
writeLines(txt, "random600_compare_report.txt")
cat(paste(txt, collapse = "\n"), "\n")
