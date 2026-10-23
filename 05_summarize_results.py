import sys, pandas as pd
pd.set_option("display.width", 200)
summ = pd.read_csv("nested_cv_full_summary.csv")
det  = pd.read_csv("nested_cv_full_detail.csv")

print("rows summary/detail:", len(summ), len(det))
print("failures:", (det.status != "ok").sum() if "status" in det else "n/a")
print("\n== SVM_linear at 600 ==")
print(summ[(summ.classifier=="SVM_linear") & (summ.n_markers==600)][["metric","mean_acc","sd_acc"]].to_string(index=False))

print("\n== best classifier at 600 per metric ==")
b = summ[summ.n_markers==600].sort_values("mean_acc", ascending=False).groupby("metric").head(2)
print(b[["metric","classifier","mean_acc","sd_acc"]].to_string(index=False))

print("\n== 95% threshold (minimum n with mean_acc >= 0.95) ==")
for cl in ["SVM_linear","RF"]:
    t = (summ[(summ.classifier==cl) & (summ.mean_acc>=0.95)]
         .groupby("metric").n_markers.min())
    print(cl, t.to_dict())

print("\n== overall max accuracy ==")
print(summ.sort_values("mean_acc", ascending=False).head(5)[["metric","n_markers","classifier","mean_acc"]].to_string(index=False))

print("\n== selected hyper-parameters (share) ==")
sel = pd.read_csv("nested_cv_selected_params.csv")
for cl in ["SVM_linear","RF"]:
    sub = sel[sel.classifier==cl]
    print(cl, sub.params.value_counts().head(6).to_dict())
