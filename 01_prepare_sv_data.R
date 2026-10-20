library(data.table)

RAW_FILE <- "SV_geno_additive.raw"

cat(format(Sys.time(), "%H:%M:%S"), "Reading PLINK raw file...\n")
raw <- fread(RAW_FILE, header = TRUE)
cat(sprintf("Raw dimensions: %d samples x %d columns\n", nrow(raw), ncol(raw)))

meta_cols <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
snp_cols  <- setdiff(colnames(raw), meta_cols)
cat(sprintf("SV sites: %d\n", length(snp_cols)))

raw[, FID := as.character(FID)]
raw[, IID := as.character(IID)]

cat(sprintf("Samples: %d\n", nrow(raw)))
cat("Breed distribution (FID column):\n")
print(table(raw$FID))

geno_mat <- as.matrix(raw[, ..snp_cols])
rownames(geno_mat) <- paste(raw$FID, raw$IID, sep = "_")
mode(geno_mat) <- "numeric"

n_na <- sum(is.na(geno_mat))
cat(sprintf("Missing genotypes: %d (%.4f%%)\n", n_na, 100 * n_na / length(geno_mat)))
if (n_na > 0) {
  cat("Imputing missing genotypes by site mean...\n")
  col_means <- colMeans(geno_mat, na.rm = TRUE)
  for (j in seq_len(ncol(geno_mat))) {
    na_idx <- which(is.na(geno_mat[, j]))
    if (length(na_idx) > 0) geno_mat[na_idx, j] <- col_means[j]
  }
}

snp_var <- apply(geno_mat, 2, var, na.rm = TRUE)
zero_var_idx <- which(snp_var == 0 | is.na(snp_var))
if (length(zero_var_idx) > 0) {
  cat(sprintf("Removed %d zero-variance sites\n", length(zero_var_idx)))
  geno_mat <- geno_mat[, -zero_var_idx, drop = FALSE]
}

breed_vec <- as.factor(raw$FID)
snp_names <- colnames(geno_mat)

cat(sprintf("\nFinal data: %d samples x %d SV sites x %d breeds\n",
            nrow(geno_mat), ncol(geno_mat), nlevels(breed_vec)))
cat("Samples per breed:\n")
print(table(breed_vec))

saveRDS(geno_mat,  "geno_mat.rds")
saveRDS(breed_vec, "breed_vec.rds")
saveRDS(snp_names, "snp_names.rds")

cat(format(Sys.time(), "%H:%M:%S"), "Done: geno_mat.rds / breed_vec.rds / snp_names.rds\n")
