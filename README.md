# SV_ml_classification

Code accompanying *Population-scale structural variants as genetic markers for duck breed tracing*
(PSJ-D-26-02294). It assigns 339 ducks from 22 populations to breeds from structural-variant
genotypes, comparing four ranking statistics (FST, Delta, AED, PIC) across ten classifiers and 20
panel sizes under nested cross-validation. Marker ranking and hyper-parameter tuning are done
inside every training fold, the ranked panels are compared with randomly drawn panels of the same
size, and significance is assessed with a stratified Monte Carlo test. The input is a PLINK
`--recode A` genotype file with the breed code in FID; the data and the result tables are not
included.
