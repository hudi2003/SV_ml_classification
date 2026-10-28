# -*- coding: utf-8 -*-
"""Figure S5: 100 random 600-SV panels vs. the four ranked 600-SV panels.
Style matched to Figure S3 / Figure 4 (same rcParams, palette, legend style)."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["font.size"] = 13
matplotlib.rcParams["axes.labelsize"] = 14
matplotlib.rcParams["xtick.labelsize"] = 12
matplotlib.rcParams["ytick.labelsize"] = 12
matplotlib.rcParams["legend.fontsize"] = 11
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FormatStrFormatter

BASE = os.environ.get("RAND_DIR", ".")
det = pd.read_csv(BASE + r"\random600_detail.csv")
rep = det.groupby("rep")["acc"].mean().values * 100
NREP = len(rep)
st = pd.read_csv(BASE + r"\random600_stats.csv").set_index("metric")
detall = pd.read_csv(os.environ.get("NESTED_DETAIL", "nested_cv_full_detail.csv"))
ref = detall[(detall.n_markers == 600) & (detall.classifier == "SVM_linear") &
             (detall.status == "ok") &
             detall.metric.isin(["FST", "AED", "Delta", "PIC"])]

col = {"FST": "#1f77b4", "Delta": "#d62728", "AED": "#2ca02c", "PIC": "#9467bd"}
mark = {"FST": "o", "Delta": "s", "AED": "^", "PIC": "D"}
ls = {"FST": (0, (4, 2)), "PIC": (0, (1.2, 1.4)), "AED": (0, (4, 2, 1, 2)),
      "Delta": (0, (6, 2, 1, 2))}
order = ["FST", "PIC", "AED", "Delta"]
lab = {"FST": r"$\mathit{F}_{\mathrm{ST}}$-ranked", "Delta": "Delta-ranked",
       "AED": "AED-ranked", "PIC": "PIC-ranked"}
plain = {"FST": "FST-ranked", "Delta": "Delta-ranked", "AED": "AED-ranked",
         "PIC": "PIC-ranked"}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.4),
                               gridspec_kw={"width_ratios": [1.18, 1]})

# ---- (a) distribution of the 100 random panels -----------------------------
ax1.hist(rep, bins=np.arange(90.3, 96.75, 0.25), color="#9e9e9e", alpha=0.45,
         edgecolor="#555555", linewidth=0.7)
ax1.axvline(rep.mean(), color="#555555", lw=1.0, ls=":")
ax1.text(rep.mean() - 0.02, 1.0, "mean of the %d draws = %.1f%%" % (NREP, rep.mean()),
         fontsize=9, color="#333333", rotation=90, va="bottom", ha="right")

handles = [Patch(facecolor="#9e9e9e", edgecolor="#555555", alpha=0.45,
                 label="random panels (n = %d)" % NREP)]
for m in order:
    v = float(st.loc[m, "ranked"])
    ax1.axvline(v, color=col[m], lw=1.5, ls=ls[m], alpha=0.95)
    handles.append(Line2D([0], [0], color=col[m], lw=1.5, ls=ls[m],
                          label="%s  %.1f%%" % (plain[m], v)))
ax1.legend(handles=handles, loc="upper left", frameon=False, fontsize=9.8,
           handlelength=2.2, labelspacing=0.35, borderaxespad=0.3)
ax1.set_xlabel("Accuracy of a random 600-SV panel")
ax1.set_ylabel("Number of random panels")
ax1.set_xlim(90.3, 96.5)
ax1.set_ylim(0, ax1.get_ylim()[1] * 1.02)
ax1.set_xticks([91, 92, 93, 94, 95, 96])
ax1.xaxis.set_major_formatter(FormatStrFormatter("0.%02d"))
ax1.grid(alpha=0.25, lw=0.6)
ax1.text(0.0, 1.02, "(a)", transform=ax1.transAxes, fontsize=15,
         va="bottom", ha="left", fontweight="bold")

# ---- (b) paired by outer fold ---------------------------------------------
folds = [1, 2, 3, 4, 5]
data = [100 * det.loc[det.fold == f, "acc"].values for f in folds]
ax2.boxplot(data, widths=0.55, patch_artist=True, showfliers=True,
            medianprops=dict(color="#333333", lw=1.2),
            boxprops=dict(facecolor="#d9d9d9", edgecolor="#555555", lw=0.9),
            whiskerprops=dict(color="#555555", lw=0.9),
            capprops=dict(color="#555555", lw=0.9),
            flierprops=dict(marker="o", ms=2.6, mfc="#555555",
                            mec="#555555", alpha=0.6))
off = {"FST": -0.19, "PIC": -0.065, "AED": 0.065, "Delta": 0.19}
for m in order:
    rv = [100 * float(ref[(ref.metric == m) & (ref.fold == f)].acc.iloc[0]) for f in folds]
    ax2.scatter([f + off[m] for f in folds], rv, marker=mark[m], s=34,
                color=col[m], zorder=5, edgecolor="white", linewidth=0.6,
                label=lab[m])
ax2.set_xlabel("Outer cross-validation fold")
ax2.set_ylabel("Breed-assignment accuracy (nested CV)")
ax2.set_ylim(81.3, 100.5)
ax2.set_yticks([85, 90, 95, 100])
ax2.set_yticklabels(["0.85", "0.90", "0.95", "1.00"])
ax2.set_xticks(folds)
ax2.grid(alpha=0.25, lw=0.6, axis="y")
ax2.text(0.0, 1.02, "(b)", transform=ax2.transAxes, fontsize=15,
         va="bottom", ha="left", fontweight="bold")
ax2.legend(loc="lower right", frameon=False, fontsize=9.8, ncol=2,
           handletextpad=0.3, borderaxespad=0.3, columnspacing=0.8)

plt.tight_layout()
out = os.environ.get("FIG_OUT", "FigureS5_random600_panels")
plt.savefig(out + ".png", dpi=300, bbox_inches="tight")
plt.savefig(out + ".pdf", bbox_inches="tight")
print("saved", out)
