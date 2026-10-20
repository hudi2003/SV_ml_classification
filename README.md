# SV_ml_classification

**Breed assignment in domestic ducks using population-scale structural variants**

Code accompanying the manuscript *"Population-scale structural variants as genetic markers for duck breed tracing"* (PSJ-D-26-02294).

---

## Overview

This repository contains the R code used to evaluate structural variants (SVs) as breed-informative markers in 339 domestic ducks from 22 populations. Four ranking statistics are compared — FST, Delta, the average Euclidean distance (AED) and the polymorphism information content (PIC) — across ten classifiers and 20 panel sizes under a nested cross-validation design.

**Headline result.** With a linear-kernel SVM, a panel of 600 FST-ranked SVs assigned breeds with 98.8% accuracy; 400 FST-ranked markers already exceeded 95%, as did 200 AED- or Delta-ranked markers. Panels of randomly selected SVs required 600 markers to reach 95.3%.

---

## Repository contents

| File | Description |
|---|---|
| `00_functions.R` | Core functions: stratified folds, Weir–Cockerham FST, Delta/AED/PIC, classifier wrappers |
| `01_prepare_sv_data.R` | Converts a PLINK `--recode A` `.raw` file into `geno_mat.rds`, `breed_vec.rds` and `snp_names.rds` |
| `02_rigorous_nested_benchmark.R` | Main benchmark: four ranking statistics x 20 panel sizes x ten classifiers |
| `03_visualization.R` | Accuracy curves, confusion matrices and best-classifier heatmap |
| `export_marker_ranking.R` | Exports per-SV values of the four ranking statistics |
| `fst_confusion_only.R` | Confusion matrix for the FST-ranked panel |

---

## Input data

The scripts expect `SV_geno_additive.raw`, produced by `plink --recode A` from the final SV genotype set:

* **339 individuals x 9,384 SVs** (8,588 deletions, 153 duplications, 643 break ends)
* FID column holds the breed code (22 populations)
* Missing genotypes are imputed with the per-site mean; zero-variance sites are removed

The input files are not distributed here because they form part of the manuscript's data.

---

## Method

### Marker ranking

Four statistics are computed for every SV and used to rank the markers:

1. **FST** — multi-population Weir & Cockerham (1984) estimator
2. **Delta** — mean pairwise allele frequency difference across all breed pairs
3. **AED** — mean Euclidean distance of allele frequencies across all breed pairs
4. **PIC** — polymorphism information content, `2pq(1 - pq)`

### Nested cross-validation

* **Outer loop** — 5-fold stratified cross-validation (`stratified_folds`, seed 42)
* **Inner loop** — 3-fold cross-validation with grid search over the SVM cost (0.01, 0.1, 1, 10) and the number of random-forest features (`sqrt`)
* **Marker ranking is recomputed inside every training fold**, so held-out samples never influence marker selection
* **Panel sizes** — 10, 20, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1,000, 1,200, 1,400, 1,600, 1,800, 2,000, 2,500, 3,000

### Classifiers

Random forest (`ntree` = 500), linear and RBF SVM, naive Bayes, LDA, k-nearest neighbours (k = 5), decision tree, XGBoost (`nrounds` = 100, `max_depth` = 6, `eta` = 0.3), multilayer perceptron (`size` = 10, `decay` = 5e-4, `maxit` = 300), and lasso-penalised multinomial logistic regression.

---

## Usage

```bash
# 1. Build the genotype matrix from the PLINK raw file
Rscript 01_prepare_sv_data.R

# 2. Run the nested benchmark (this is the long step)
Rscript 02_rigorous_nested_benchmark.R

# 3. Produce the figures
Rscript 03_visualization.R
```

Tested with R 4.3.0. Required packages: `data.table`, `dplyr`, `randomForest`, `e1071`, `MASS`, `class`, `rpart`, `xgboost`, `nnet`, `glmnet`, `foreach`, `doParallel`.

---

## Citation

If you use this code, please cite:

> Hu D, Zhang Y, Fan Q, Hu Z, Hou S, Wang X. Population-scale structural variants as genetic markers for duck breed tracing. *Poultry Science* (under revision).

---

## License

Released for academic use. Please contact the corresponding author (Xia Wang, xiawang@nwafu.edu.cn) regarding reuse.
