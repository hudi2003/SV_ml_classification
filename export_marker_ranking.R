#!/usr/bin/env Rscript
# ================================================================
# 脚本：export_marker_ranking.R
# 功能：基于全量样本（不分折），一次性计算 FST / Delta / AED / PIC
#       四种标记信息量指标，并按PIC降序排序，导出完整排序表，
#       用于专利交底书"表2：结构变异标记信息量排序"。
#
# 依赖：00_functions.R 中的 compute_fst_wc() 与 compute_delta_aed_pic()
# 用法：在 ml_classification/stacking_analysis/ 目录下运行：
#       Rscript export_marker_ranking.R
# ================================================================

library(data.table)        # 用于 fwrite
source("00_functions.R")   # 复用 compute_fst_wc / compute_delta_aed_pic

geno_mat  <- readRDS("geno_mat.rds")
breed_vec <- as.factor(readRDS("breed_vec.rds"))
snp_names <- readRDS("snp_names.rds")

cat(format(Sys.time(), "%H:%M:%S"), "基于全量数据计算FST/Delta/AED/PIC...\n")

fst_all <- compute_fst_wc(geno_mat, breed_vec)
dap_all <- compute_delta_aed_pic(geno_mat, breed_vec)

ranking_tab <- data.frame(
  SV_ID = snp_names,
  FST   = fst_all,
  Delta = dap_all$Delta,
  AED   = dap_all$AED,
  PIC   = dap_all$PIC,
  stringsAsFactors = FALSE
)

# ---- 从SV_ID中解析染色体、起止坐标（若ID形如 Chr01_181638_181849_... ）----
id_parts <- strsplit(ranking_tab$SV_ID, "_")
ranking_tab$Chr   <- sapply(id_parts, function(x) x[1])
ranking_tab$Start <- suppressWarnings(sapply(id_parts, function(x) as.integer(x[2])))
ranking_tab$End   <- suppressWarnings(sapply(id_parts, function(x) as.integer(x[3])))
ranking_tab$SVLEN <- abs(ranking_tab$End - ranking_tab$Start)

# ---- 按PIC降序排序（如需按其他指标排序，把 order(-ranking_tab$PIC) 换成对应列即可）----
ranking_tab <- ranking_tab[order(-ranking_tab$PIC), ]
ranking_tab$PIC_rank <- seq_len(nrow(ranking_tab))

# 同时给出FST/Delta/AED各自的排名，供需要时切换排序依据
ranking_tab$FST_rank   <- rank(-ranking_tab$FST,   ties.method = "min")
ranking_tab$Delta_rank <- rank(-ranking_tab$Delta, ties.method = "min")
ranking_tab$AED_rank   <- rank(-ranking_tab$AED,   ties.method = "min")

out_cols <- c("PIC_rank", "SV_ID", "Chr", "Start", "End", "SVLEN",
              "PIC", "FST", "Delta", "AED",
              "FST_rank", "Delta_rank", "AED_rank")
fwrite(ranking_tab[, out_cols], "marker_ranking_full.csv")

cat(format(Sys.time(), "%H:%M:%S"), "完整排序表已保存: marker_ranking_full.csv\n")
cat(sprintf("共 %d 个标记，已按PIC降序排列。\n", nrow(ranking_tab)))

# ---- 额外导出前100个标记，可直接用于专利交底书表2 ----
fwrite(ranking_tab[1:min(100, nrow(ranking_tab)), out_cols], "marker_ranking_top100_PIC.csv")
cat("前100个标记（按PIC排序）已保存: marker_ranking_top100_PIC.csv\n")
