import pandas as pd, numpy as np
d = pd.read_csv(r"C:/Users/hudi/rev_work/SV_ml_classification/random600/random600_detail.csv")
g = d.groupby("rep")["acc"].mean()*100
st = pd.read_csv(r"C:/Users/hudi/rev_work/SV_ml_classification/random600/random600_stats.csv")
for _, r in st.iterrows():
    diff = r["ranked"] - g
    print(f'{r["metric"]:5s} mean diff {diff.mean():+.2f} pp  sd {diff.std(ddof=1):.2f}  '
          f'95%CI [{np.percentile(diff,2.5):+.2f},{np.percentile(diff,97.5):+.2f}]  '
          f'share of draws beaten {(diff>0).mean()*100:.0f}%')
print("draw means q:", np.round(np.percentile(g,[0,2.5,25,50,75,97.5,100]),2))
print("draws >=95%:", int((g>=95).sum()), ">=94%:", int((g>=94).sum()))
