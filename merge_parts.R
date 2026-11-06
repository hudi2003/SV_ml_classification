suppressPackageStartupMessages({ library(data.table); library(dplyr) })
files <- list.files("parts", pattern = "^row_[0-9]+\\.csv$", full.names = TRUE)
cat("row files:", length(files), "\n")
detailed <- rbindlist(lapply(files, fread), fill = TRUE)
setorder(detailed, fold, metric, n_markers, classifier)
fwrite(detailed, "nested_cv_detail_full.csv")

detailed$params[is.na(detailed$params) | detailed$params == ""] <- "fixed"
results <- detailed %>%
  group_by(metric, n_markers, classifier) %>%
  summarise(mean_acc = mean(acc, na.rm = TRUE),
            sd_acc = sd(acc, na.rm = TRUE),
            mean_time = mean(time_sec, na.rm = TRUE),
            n_ok = sum(status == "ok"), .groups = "drop")
sel_tab <- detailed %>% filter(status == "ok") %>%
  select(fold, metric, n_markers, classifier, params, inner_acc, inner_sd, inner_grid_size)

fwrite(results, "nested_cv_summary_ALL.csv")
fwrite(detailed, "nested_cv_detail_ALL.csv")
fwrite(sel_tab, "nested_cv_params_ALL.csv")
cat("combos:", nrow(detailed), "| errors:", sum(detailed$status != "ok"), "\n")
cat("summary rows:", nrow(results), "\n")
