# Random-panel control at 600 markers (1,000 draws)

Numbers written into the revised manuscript (Methods 2.8, Results 3.4, Figure S5) and
into the Major 5 response. Generated on mserver with R 4.3.3, 48 workers, 18.4 min.

Design: 8,741 autosomal DEL+DUP (643 BND excluded), 339 samples, 22 breeds; outer
5-fold stratified CV (seed 42), inner 5-fold stratified CV with the linear-SVM cost
grid {0.01, 0.1, 1, 10}; the 600-SV panel is re-drawn inside every outer training
fold and every inner training split; 1,000 independent draws.

| Panel | Accuracy | SD (5 folds) | Gain vs random | Random draws >= it | Paired t (5 folds) |
|---|---|---|---|---|---|
| FST-ranked | 95.29% | 1.86 | +1.59 pp | 61 / 1,000 | P = 0.053 |
| PIC-ranked | 95.31% | 1.82 | +1.61 pp | 55 / 1,000 | P = 0.024 |
| AED-ranked | 94.44% | 2.05 | +0.75 pp | 201 / 1,000 | P = 0.196 |
| Delta-ranked | 93.85% | 1.73 | +0.15 pp | 437 / 1,000 | P = 0.873 |
| Random (1,000 draws) | 93.69% | 1.04 (between draws) | - | - | - |

Random draws: range 89.39-97.09%, 2.5-97.5th percentile 91.71-95.81%; the single
random ordering used in the main benchmark (93.51%) lies inside the distribution.
A random panel shares 41.3 +/- 6.0 of its 600 SVs with the FST top-600 panel of the
same fold (6.9%), against 6.9% expected by chance.

Outputs: `random600_detail.csv` (rep x fold), `random600_stats.csv`,
`random600_vs_ranked.csv`, `random600_summary.txt`,
`FigureS5_random600_panels.png/.pdf`.
