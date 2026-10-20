# SV-based breed assignment in domestic ducks

Code for the marker-ranking and machine-learning analyses in
**"Population-scale structural variants as genetic markers for duck breed tracing"** (PSJ-D-26-02294).

## Contents

| File | Purpose |
|---|---|
| `00_functions.R` | Stratified fold construction, Weir–Cockerham FST, Delta/AED/PIC, classifier wrappers |
| `01_prepare_sv_data.R` | Converts a PLINK `--recode A` `.raw` file into `geno_mat.rds`, `breed_vec.rds`, `snp_names.rds` |
| `02_rigorous_nested_benchmark.R` | Main benchmark: panels ranked by FST / Delta / AED / PIC, ten classifiers, 20 panel sizes |
| `03_visualization.R` | Accuracy curves, confusion matrices, best-classifier heatmap |
| `export_marker_ranking.R` | Exports per-SV ranking statistics |
| `fst_confusion_only.R` | Confusion matrix for the FST-ranked panel |

## Input data

`SV_geno_additive.raw` is the PLINK additive-encoding output of `plink --recode A` applied to
`SV_filtered_final.bed/.bim/.fam`:

* 339 individuals x 9,384 autosomal SVs (8,588 DEL, 153 DUP, 643 BND)
* FID column contains the breed code (22 populations)
* missing genotypes are imputed with the per-site mean; zero-variance sites are dropped

The input files are not included in this repository.

## Resampling design

* Outer loop: 5-fold stratified cross-validation (`stratified_folds`, seed 42)
* Inner loop: 3-fold cross-validation with grid search over the SVM cost
  (C in 0.01, 0.1, 1, 10) and the number of random-forest features (`sqrt`)
* Marker ranking is **recomputed inside every training fold**, so held-out samples
  contribute no information to marker selection
* Panel sizes: 10, 20, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
  1200, 1400, 1600, 1800, 2000, 2500, 3000

## Classifiers

Random forest (`ntree` = 500), linear and RBF SVM, naive Bayes, LDA, k-NN (k = 5),
decision tree, XGBoost (`nrounds` = 100, `max_depth` = 6, `eta` = 0.3),
multilayer perceptron (`size` = 10, `decay` = 5e-4, `maxit` = 300),
and lasso-penalised multinomial logistic regression.

## Usage

```bash
Rscript 01_prepare_sv_data.R
Rscript 02_rigorous_nested_benchmark.R
Rscript 03_visualization.R
```

R 4.3.0 with data.table, dplyr, randomForest, e1071, MASS, class, rpart, xgboost, nnet,
glmnet, foreach and doParallel.
