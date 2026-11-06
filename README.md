# SV_ml_classification

Code for *Population-scale structural variants as genetic markers for duck breed tracing*
(PSJ-D-26-02294).

Breed assignment from structural variants in 339 ducks (22 populations). Four ranking statistics
(FST, Delta, AED, PIC) are compared across ten classifiers and 20 panel sizes under nested
cross-validation, with a random-panel control and a stratified Monte Carlo test for the ranking
step.

## Files

| File | Purpose |
|---|---|
| `00_functions.R` | Folds, ranking statistics, classifier wrappers, nested CV |
| `01_prepare_sv_data.R` | PLINK `--recode A` output to `geno_mat.rds` / `breed_vec.rds` / `snp_names.rds` |
| `02_rigorous_nested_benchmark.R` | Main benchmark (4 statistics x 20 panel sizes x 10 classifiers) |
| `02_single.R`, `make_tasks.py`, `merge_parts.R` | Run one cell, build the task list, merge the parts |
| `03_visualization.R`, `03b`, `03c` | Accuracy curves, confusion matrices, Figure 4 |
| `04_verify_no_leakage.R` | Split disjointness, invariance of the inner search, permutation control |
| `05_summarize_results.py` | Numbers reported in the manuscript |
| `panel_selected_by_inner_cv.py` | Panel size chosen by the inner CV |
| `random_panels/*` | Random-panel control, Monte Carlo test, repeated outer CV |

## Run

```bash
Rscript 01_prepare_sv_data.R
Rscript 02_rigorous_nested_benchmark.R
Rscript 03_visualization.R

cd random_panels
RAND_N=600 RAND_REP=1000 NESTED_CORES=48 Rscript random_panel_replicates.R
cd ..
python random_panels/stratified_monte_carlo.py --random random600_detail.csv \
  --detail nested_cv_detail_ALL.csv --n 600
```

## Method in brief

* Outer loop: 5-fold stratified CV (seed 42). Inner loop: 5-fold CV inside each training fold,
  selecting hyper-parameters. Ranking statistics are recomputed inside every training fold and
  every inner split, so no test sample influences marker choice or tuning.
* Random-panel control: panels of the same size are redrawn inside the same folds and put through
  the same tuning, 1,000 draws (10,000 for the 600-SV panel).
* The ranked panel has only five outer-fold accuracies, so significance comes from a stratified
  Monte Carlo test rather than a t-test: one random panel per fold, averaged, 100,000 resamples.
  `repeated_outer_cv.R` re-runs the whole procedure on many independent splits.

## Requirements

R 4.3.3 with `data.table`, `dplyr`, `randomForest`, `e1071`, `MASS`, `class`, `rpart`,
`xgboost`, `nnet`, `glmnet`, `foreach`, `doParallel`; Python 3 with `pandas`, `numpy`, `scipy`.

## Data

Input is `SV_geno_additive.raw` from `plink --recode A` (339 individuals x 9,384 SVs, breed code
in FID). Genotypes are mean-imputed per site. Neither the input data nor the result tables are
included here.

## Citation

> Hu D, Zhang Y, Fan Q, Hu Z, Hou S, Wang X. Population-scale structural variants as genetic
> markers for duck breed tracing. *Poultry Science* (under revision).

## License

For academic use. Contact the corresponding author (Xia Wang, xiawang@nwafu.edu.cn) about reuse.
