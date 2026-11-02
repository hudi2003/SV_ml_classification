#!/usr/bin/env python3
import csv
panels = [10,20,50,100,200,300,400,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2500,3000]
metrics = ["FST","Delta","AED","PIC","Random"]
cls = ["RF","SVM_linear","SVM_RBF","NB","LDA","KNN","DecisionTree","XGBoost","MLP","Lasso"]
w = dict(RF=6, SVM_linear=1, SVM_RBF=1.5, NB=0.5, LDA=0.2, KNN=0.2,
         DecisionTree=1, XGBoost=2, MLP=30, Lasso=1)
rows = []
for k in range(1,6):
    for m in metrics:
        for n in panels:
            for c in cls:
                rows.append((k, m, n, c, n*w[c]))
rows.sort(key=lambda r: -r[4])
with open("tasks.csv","w",newline="") as f:
    wr = csv.writer(f); wr.writerow(["idx","fold","metric","n_markers","classifier","cost"])
    for i,r in enumerate(rows,1):
        wr.writerow([i,r[0],r[1],r[2],r[3],round(r[4],1)])
print("tasks:", len(rows))
