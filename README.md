# SV_ml_classification

**Breed assignment in domestic ducks using population-scale structural variants**

Code accompanying the manuscript *"Population-scale structural variants as genetic markers for duck breed tracing"* (PSJ-D-26-02294).

---

## Overview

This repository contains the R code used to evaluate structural variants (SVs) as breed-informative markers in 339 domestic ducks from 22 populations. Four ranking statistics are compared — FST, Delta, the average Euclidean distance (AED) and the polymorphism information content (PIC) — across ten classifiers and 20 panel sizes under a nested cross-validation design.

**Headline result.** With a linear-kernel SVM, a panel of 600 FST-ranked SVs assigned breeds with 95.3% accuracy (94.4% for AED-, 93.8% for Delta- and 95.3% for PIC-ranked SVs); the 95% level was exceeded with 600 FST-, 800 AED-, 800 Delta- and 600 PIC-ranked markers. Randomly selected panels of the same size reached 93.7 +/- 1.0% over 1,000 independent draws, so ranking buys a modest (0.2-1.6 percentage point) but reproducible gain rather than being the source of the accuracy; random panels required 1,200 markers to reach 95%.

---

## Repository contents

| File | Description |
|---|---|
| `00_functions.R` | Core functions: stratified folds, Weir–Cockerham FST, Delta/AED/PIC, classifier wrappers |
| `01_prepare_sv_data.R` | Converts a PLINK `--recode A` `.raw` file into `geno_mat.rds`, `breed_vec.rds` and `snp_names.rds` |
| `02_rigorous_nested_benchmark.R` | Main benchmark: four ranking statistics x 20 panel sizes x ten classifiers |
| `03_visualization.R` | Accuracy curves, confusion matrices, economical-panel table and best-classifier heatmap |
| `04_verify_no_leakage.R` | Leakage checks: split disjointness, invariance of the inner search to the test fold, permuted-label negative control |
| `05_summarize_results.py` | Extracts the numbers cited in the manuscript from the benchmark tables |
| `random_panels/random_panel_replicates.R` | Random-panel control for any panel size and draw count (`RAND_N`, `RAND_REP`) |
| `random_panels/stratified_monte_carlo.py` | Stratified Monte Carlo test, effect size, 95% CI and t-test sensitivities |
| `random_panels/repeated_outer_cv.R` | Repeats the whole nested procedure on many independently drawn 5-fold splits |

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
* **Inner loop** — 5-fold stratified cross-validation inside every outer training fold, with a grid search over the SVM cost (0.01, 0.1, 1, 10) and over the random forest (300, 500 or 1,000 trees x mtry = sqrt(p), p/5 or p/2)
* **Marker ranking is recomputed inside every training fold and inside every inner training split**, so neither held-out nor inner-validation samples influence marker selection or the hyper-parameter choice
* **Panel sizes** — 10, 20, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1,000, 1,200, 1,400, 1,600, 1,800, 2,000, 2,500, 3,000

### Random-panel control

The control for the ranking step uses exactly the same nested design: a panel of N SVs (400 or 600) is re-drawn at random inside every outer training fold and again inside every inner training split, and the whole procedure is repeated for every draw (1,000 draws, and 10,000 for the 600-SV panel). Accuracy across draws is summarised as mean +/- SD and compared with the four ranked panels of the same size.

Significance is assessed with a stratified Monte Carlo test rather than a t-test. The ranked panel yields only five outer-fold accuracies, which are not independent (the folds share about 80% of their training data) and are discrete proportions bounded at 100%; a t distribution with four degrees of freedom is therefore both thick-tailed and built on an assumption the data do not meet. Instead, one random panel is drawn from each fold, the five accuracies are averaged, and the draw is repeated 100,000 times; this is the distribution of the five-fold mean when the panel carries no ranking information. The p value is the fraction of these draws at least as large as the observed panel. The one-sample and paired t-tests are reported alongside as sensitivity analyses.

### Classifiers

Random forest (ntree and mtry tuned in the inner loop), linear SVM (cost tuned in the inner loop), RBF SVM, naive Bayes, LDA, k-nearest neighbours (k = 5), decision tree, XGBoost (`nrounds` = 100, `max_depth` = 6, `eta` = 0.3), multilayer perceptron (`size` = 10, `decay` = 5e-4, `maxit` = 300), and lasso-penalised multinomial logistic regression.

---

## Usage

```bash
# 1. Build the genotype matrix from the PLINK raw file
Rscript 01_prepare_sv_data.R

# 2. Run the nested benchmark (this is the long step)
Rscript 02_rigorous_nested_benchmark.R

# 3. Produce the figures
Rscript 03_visualization.R

# 4. Random-panel control and its significance test
cd random_panels
RAND_N=600 RAND_REP=1000 NESTED_CORES=48 Rscript random_panel_replicates.R
cd ..
python random_panels/stratified_monte_carlo.py \
  --random random600_detail.csv \
  --detail ../nested_cv_detail_ALL.csv --n 600

# 5. Optional: repeated outer cross-validation (1,000 splits, 32 workers)
cd random_panels
REP_SPLITS=1000 REP_RAND=2 NESTED_CORES=32 Rscript repeated_outer_cv.R
```

Tested with R 4.3.3. Required packages: `data.table`, `dplyr`, `randomForest`, `e1071`, `MASS`, `class`, `rpart`, `xgboost`, `nnet`, `glmnet`, `foreach`, `doParallel`.

---

## Citation

If you use this code, please cite:

> Hu D, Zhang Y, Fan Q, Hu Z, Hou S, Wang X. Population-scale structural variants as genetic markers for duck breed tracing. *Poultry Science* (under revision).

---

## License

Released for academic use. Please contact the corresponding author (Xia Wang, xiawang@nwafu.edu.cn) regarding reuse.

## All-classifier tuning (2026-11)

Hyper-parameters are now tuned by the inner cross-validation for all ten classifiers (grids in
`00_functions.R`): linear SVM C; RBF SVM C x gamma; RF trees x mtry; k-NN k; naive Bayes Laplace
smoothing; decision tree cp (minsplit = 5, no surrogate splits); XGBoost rounds x depth; MLP hidden
units x weight decay. LDA has no hyper-parameters and the lasso penalty is chosen by `cv.glmnet`
inside the training fold; the lasso branch also sweeps the elastic-net alpha from 0 to 1.

At 600 FST-ranked markers the linear SVM remains the most accurate (95.3%), followed by the RBF SVM
(94.7%) and the MLP (93.0%); LDA reaches 90.3%, k-NN 90.0%, the lasso 89.5%, the random forest 85.0%,
XGBoost 76.8% and the decision tree 50.5%. SVM_linear, RF, LDA and NB reproduce the earlier
fixed-setting run exactly, which also verifies reproducibility.

Each of the 5,000 (fold x ranking statistic x panel size x classifier) cells is run as an independent
single-core task:

| File | Description |
|---|---|
| `02_single.R` | Runs one cell; writes `parts/row_XXXXX.csv` |
| `make_tasks.py` | Builds `tasks.csv` (5,000 cells, most expensive first) |
| `merge_parts.R` | Merges `parts/` into the `nested_cv_*_ALL.csv` tables |
| `panel_selected_by_inner_cv.py` | Panel size chosen inside each training fold by the inner CV |
The merged tables (`nested_cv_summary_ALL.csv`, `nested_cv_detail_ALL.csv`, `nested_cv_params_ALL.csv`) and `panel_selected_by_inner_cv.csv` are written by these scripts. They are not tracked in the repository.


## Random-panel control (2026-11)

600-SV panel, 10,000 draws, nested 5x5, linear SVM. Random panels reached 93.75 +/- 1.03%
(2.5-97.5%: 91.74-95.66%; 1,114 of 10,000 draws reached 95%). One-sided stratified Monte Carlo
P values (100,000 resamples): FST 0.065, PIC 0.056, AED 0.232, Delta 0.455; gains of +1.53,
+1.55, +0.69 and +0.09 percentage points over the random mean.

At 400 markers (1,000 draws) the random mean is 91.48 +/- 1.19% and the gain is larger:
FST 93.81% (+2.34 pp, P = 0.024, positive in all five folds) and Delta 93.56% (+2.08 pp,
P = 0.037); AED and PIC are not significant. Ranking therefore matters more when the panel is
small, and the gain is modest in absolute terms.

The scripts in `random_panels/` accept any panel size through `RAND_N`, so the 400- and
600-SV controls are the same code with a different environment variable. The 10,000-draw run
used the same seeds as the 1,000-draw run for draws 1-1,000, which reproduces the earlier numbers
exactly.
