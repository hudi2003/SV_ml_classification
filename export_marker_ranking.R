library(data.table)
source("00_functions.R")

geno_mat  <- readRDS("geno_mat.rds")
breed_vec <- as.factor(readRDS("breed_vec.rds"))
snp_names <- readRDS("snp_names.rds")

cat(format(Sys.time(), "%H:%M:%S"), "Computing FST/Delta/AED/PIC on the full dataset...\n")

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

id_parts <- strsplit(ranking_tab$SV_ID, "_")
ranking_tab$Chr   <- sapply(id_parts, function(x) x[1])
ranking_tab$Start <- suppressWarnings(sapply(id_parts, function(x) as.integer(x[2])))
ranking_tab$End   <- suppressWarnings(sapply(id_parts, function(x) as.integer(x[3])))
ranking_tab$SVLEN <- abs(ranking_tab$End - ranking_tab$Start)

ranking_tab <- ranking_tab[order(-ranking_tab$PIC), ]
ranking_tab$PIC_rank <- seq_len(nrow(ranking_tab))

ranking_tab$FST_rank   <- rank(-ranking_tab$FST,   ties.method = "min")
ranking_tab$Delta_rank <- rank(-ranking_tab$Delta, ties.method = "min")
ranking_tab$AED_rank   <- rank(-ranking_tab$AED,   ties.method = "min")

out_cols <- c("PIC_rank", "SV_ID", "Chr", "Start", "End", "SVLEN",
              "PIC", "FST", "Delta", "AED",
              "FST_rank", "Delta_rank", "AED_rank")
fwrite(ranking_tab[, out_cols], "marker_ranking_full.csv")

cat(format(Sys.time(), "%H:%M:%S"), "Full ranking saved: marker_ranking_full.csv\n")
cat(sprintf("%d markers sorted by PIC.\n", nrow(ranking_tab)))

fwrite(ranking_tab[1:min(100, nrow(ranking_tab)), out_cols], "marker_ranking_top100_PIC.csv")
cat("Top 100 markers by PIC saved.\n")
