#!/bin/bash
cd "$HOME/hudi/SV_ml_nested" || exit 1
export R_LIBS_USER="$HOME/hudi/Rlib"
export OMP_NUM_THREADS=1
export RAND_REP=100
export RAND_N=600
export NESTED_CORES=48
export NESTED_INNER_K=5
Rscript 06_random600_replicates.R > random600.log 2>&1
echo "EXIT=$?" >> random600.log
