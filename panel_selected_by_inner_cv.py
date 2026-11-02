"""Panel size chosen by the inner cross-validation (reviewer Major 5, point 2).

For every (metric, classifier) and every outer fold, the panel size is selected
inside the training fold as the smallest n whose inner-CV accuracy reaches 95%
(falling back to the n with the highest inner accuracy when 95% is not reached).
The selected panel is then evaluated on the held-out outer test fold, so the
reported accuracy is an honest estimate of the whole 'choose the panel, then
classify' procedure.
"""
import pandas as pd, numpy as np, os

P = r"C:/Users/hudi/rev_work/SV_ml_classification/"
d = pd.read_csv(P + "nested_cv_detail_ALL.csv")
d = d[(d.status == "ok") & d.inner_acc.notna()].copy()
print("rows usable:", len(d), "classifiers:", d.classifier.nunique())

rows = []
for (m, c), g in d.groupby(["metric", "classifier"]):
    picks = []
    for f, gf in g.groupby("fold"):
        gf = gf.sort_values("n_markers")
        hit = gf[gf.inner_acc >= 0.95]
        if len(hit):
            row, src = hit.iloc[0], "95%"
        else:
            row, src = gf.loc[gf.inner_acc.idxmax()], "best"
        picks.append((int(row.n_markers), float(row.acc), src, float(row.inner_acc)))
    ns = [p[0] for p in picks]
    rows.append(dict(metric=m, classifier=c,
                     n_selected_median=int(np.median(ns)),
                     n_selected_range="%d-%d" % (min(ns), max(ns)),
                     outer_acc=np.mean([p[1] for p in picks]) * 100,
                     outer_sd=np.std([p[1] for p in picks], ddof=1) * 100,
                     reached95=sum(1 for p in picks if p[2] == "95%"),
                     mean_inner_acc=np.mean([p[3] for p in picks]) * 100))
t = pd.DataFrame(rows).sort_values(["classifier", "metric"])
t.to_csv(P + "panel_selected_by_inner_cv.csv", index=False)
print(t.to_string(index=False, float_format=lambda x: f"{x:.1f}"))

# fixed-600 reference for the same cells
fixed = (d[d.n_markers == 600].groupby(["metric", "classifier"]).acc.agg(["mean", "std"]) * 100).round(1)
fixed.columns = ["acc600_mean", "acc600_sd"]
print("\n== fixed 600 markers ==")
print(fixed.reset_index().pivot(index="classifier", columns="metric", values="acc600_mean").to_string())
