import argparse
import glob
import sys

import numpy as np
import pandas as pd
from scipy import stats

METRICS = ["FST", "Delta", "AED", "PIC"]


def read_random(patterns):
    files = []
    for p in patterns:
        hits = sorted(glob.glob(p))
        files.extend(hits if hits else [p])
    if not files:
        sys.exit("no random-panel files matched")
    frames = [pd.read_csv(f) for f in files]
    d = pd.concat(frames, ignore_index=True)
    if "status" in d.columns:
        bad = int((d["status"] != "ok").sum())
        if bad:
            print("WARNING: dropping %d failed fold evaluations" % bad, flush=True)
        d = d[d["status"] == "ok"]
    print("random panels: %d file(s), %d rows, %d draws"
          % (len(files), len(d), d["rep"].nunique()), flush=True)
    return d


def stratified_mc(fold_values, observed, B, seed):
    """One random panel per fold per resample."""
    rng = np.random.default_rng(seed)
    k = len(fold_values)
    n = min(len(v) for v in fold_values)
    mat = np.column_stack([v[:n] for v in fold_values])       # (n_draws, k)
    idx = rng.integers(0, n, size=(B, k))
    null = mat[idx, np.arange(k)].mean(axis=1)
    p_one = (1.0 + np.sum(null >= observed)) / (B + 1.0)
    lo, hi = np.percentile(null, [2.5, 97.5])
    return dict(p_one=p_one, p_two=min(1.0, 2 * p_one),
                null_mean=float(null.mean()),
                ci_lo=float(observed - hi), ci_hi=float(observed - lo))


def main():
    ap = argparse.ArgumentParser(description="Stratified Monte Carlo test of ranked panels against random panels",
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--random", nargs="+", required=True,
                    help="random-panel detail CSV(s); globs allowed")
    ap.add_argument("--detail", required=True,
                    help="nested_cv_full_detail.csv from the main benchmark")
    ap.add_argument("--n", type=int, default=None, help="panel size")
    ap.add_argument("--classifier", default="SVM_linear")
    ap.add_argument("--metrics", nargs="+", default=METRICS)
    ap.add_argument("--B", type=int, default=100000, help="resamples (default 100000)")
    ap.add_argument("--seed", type=int, default=20261104)
    ap.add_argument("--out", default=None, help="output CSV")
    args = ap.parse_args()

    rnd = read_random(args.random)
    if args.n is None:
        args.n = int(rnd["n_markers"].mode().iloc[0])
    det = pd.read_csv(args.detail)
    if "status" in det.columns:
        det = det[det["status"] == "ok"]
    ref = det[(det["n_markers"] == args.n) &
              (det["classifier"] == args.classifier) &
              (det["metric"].isin(args.metrics))]

    # panel-level score: mean of the five outer folds
    draw = rnd.groupby("rep")["acc"].mean().values * 100
    folds = sorted(rnd["fold"].unique())
    fold_values = [rnd.loc[rnd["fold"] == k, "acc"].values * 100 for k in folds]
    fold_means = np.array([v.mean() for v in fold_values])

    print("panel size %d | %d draws | draw mean %.2f +/- %.2f | "
          "2.5-97.5%% %.2f-%.2f"
          % (args.n, len(draw), draw.mean(), draw.std(ddof=1),
             np.percentile(draw, 2.5), np.percentile(draw, 97.5)), flush=True)
    print("draws reaching 95%%: %d / %d (%.1f%%)"
          % ((draw >= 95).sum(), len(draw), 100 * (draw >= 95).mean()), flush=True)

    rows = []
    for m in args.metrics:
        sub = ref[ref["metric"] == m].sort_values("fold")
        if len(sub) != len(folds):
            print("skip %s: %d fold rows" % (m, len(sub)), flush=True)
            continue
        a = sub["acc"].values * 100
        T = float(a.mean())
        mc = stratified_mc(fold_values, T, args.B, args.seed)

        # sensitivity analyses, not the primary test
        t1, p1 = stats.ttest_1samp(a, draw.mean())
        diff = a - fold_means[:len(a)]
        t2, p2 = stats.ttest_rel(a, fold_means[:len(a)])

        rows.append(dict(
            n_markers=args.n, metric=m,
            ranked=round(T, 2), ranked_sd=round(float(a.std(ddof=1)), 2),
            random=round(float(draw.mean()), 2),
            random_sd=round(float(draw.std(ddof=1)), 2),
            diff_pp=round(T - float(draw.mean()), 2),
            p_mc_one_sided=round(mc["p_one"], 4),
            p_mc_two_sided=round(mc["p_two"], 4),
            ci_lo=round(mc["ci_lo"], 2), ci_hi=round(mc["ci_hi"], 2),
            draws_ge=int((draw >= T).sum()), n_draws=int(len(draw)),
            pctile=round(100 * float((draw < T).mean()), 1),
            p_one_sample_t=round(float(p1), 4),
            p_paired_t=round(float(p2), 4)))

    out = pd.DataFrame(rows)
    cols = ["n_markers", "metric", "ranked", "random", "diff_pp",
            "p_mc_one_sided", "p_mc_two_sided", "ci_lo", "ci_hi",
            "draws_ge", "n_draws", "pctile", "p_one_sample_t", "p_paired_t"]
    out = out[cols]
    print(out.to_string(index=False), flush=True)
    if args.out:
        out.to_csv(args.out, index=False)
        print("saved " + args.out, flush=True)


if __name__ == "__main__":
    main()
