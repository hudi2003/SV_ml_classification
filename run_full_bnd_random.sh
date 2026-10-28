#!/bin/bash
cd "$HOME/hudi/SV_ml_nested" || exit 1
export NESTED_QUICK=0
export NESTED_CORES=64
export NESTED_INNER_K=5
export OMP_NUM_THREADS=1
export R_LIBS_USER="$HOME/hudi/Rlib"
rm -f nested_cv_full_summary.csv nested_cv_full_detail.csv nested_cv_selected_params.csv nested_cv_run_info.txt
Rscript 02_rigorous_nested_benchmark.R > full.log 2>&1
echo "EXIT=$?" >> full.log
