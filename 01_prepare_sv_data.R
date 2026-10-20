#!/usr/bin/env Rscript
# ================================================================
# 脚本：01_prepare_sv_data.R (v2 - 适配品种化ID版本)
# 功能：把 PLINK --recode A 输出的 .raw 文件转换成
#       geno_mat.rds / breed_vec.rds / snp_names.rds
#
# 前提: SV_filtered_final.fam 的 FID 列已经是品种代码
#       (由 run_filter_popgen.sh 里的 --update-ids 步骤生成)
#       所以品种信息直接从 .raw 文件的 FID 列取，不再需要额外的品种映射文件。
# ================================================================

library(data.table)

RAW_FILE <- "SV_geno_additive.raw"

cat(format(Sys.time(), "%H:%M:%S"), "读取PLINK raw文件...\n")
raw <- fread(RAW_FILE, header = TRUE)
cat(sprintf("原始维度: %d 样本 x %d 列(含前6列meta信息)\n", nrow(raw), ncol(raw)))

meta_cols <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
snp_cols  <- setdiff(colnames(raw), meta_cols)
cat(sprintf("SV位点数: %d\n", length(snp_cols)))

# ---- 品种信息直接取自 FID 列 ----
raw[, FID := as.character(FID)]
raw[, IID := as.character(IID)]

cat(sprintf("样本数: %d\n", nrow(raw)))
cat("品种分布(来自FID列):\n")
print(table(raw$FID))

# ---- 构建 geno_mat ----
geno_mat <- as.matrix(raw[, ..snp_cols])
rownames(geno_mat) <- paste(raw$FID, raw$IID, sep = "_")
mode(geno_mat) <- "numeric"

# ---- 缺失值填补 ----
n_na <- sum(is.na(geno_mat))
cat(sprintf("缺失基因型数: %d (%.4f%%)\n", n_na, 100 * n_na / length(geno_mat)))
if (n_na > 0) {
  cat("用各位点均值填补缺失值...\n")
  col_means <- colMeans(geno_mat, na.rm = TRUE)
  for (j in seq_len(ncol(geno_mat))) {
    na_idx <- which(is.na(geno_mat[, j]))
    if (length(na_idx) > 0) geno_mat[na_idx, j] <- col_means[j]
  }
}

# ---- 剔除零方差位点 ----
snp_var <- apply(geno_mat, 2, var, na.rm = TRUE)
zero_var_idx <- which(snp_var == 0 | is.na(snp_var))
if (length(zero_var_idx) > 0) {
  cat(sprintf("剔除 %d 个零方差位点(所有样本基因型相同)\n", length(zero_var_idx)))
  geno_mat <- geno_mat[, -zero_var_idx, drop = FALSE]
}

breed_vec <- as.factor(raw$FID)
snp_names <- colnames(geno_mat)

cat(sprintf("\n最终数据: %d 样本 x %d SV位点 x %d 品种\n",
            nrow(geno_mat), ncol(geno_mat), nlevels(breed_vec)))
cat("各品种样本数:\n")
print(table(breed_vec))

saveRDS(geno_mat,  "geno_mat.rds")
saveRDS(breed_vec, "breed_vec.rds")
saveRDS(snp_names, "snp_names.rds")

cat(format(Sys.time(), "%H:%M:%S"), "数据准备完成: geno_mat.rds / breed_vec.rds / snp_names.rds\n")
