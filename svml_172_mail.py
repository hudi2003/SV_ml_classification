#!/usr/bin/env python3
"""Progress mail for the 5000-task nested-CV run on 172 (runs on the local PC).

Remote side: metrics.sh prints STATE / DONE / WORK / ERRS / NONRUN / SUBMITTED / IDS.
Progress is reported both as finished combinations and as finished workload
(cost units = panel size x classifier weight), because the expensive MLP tasks
dominate the remaining wall time.
"""
import subprocess, smtplib, datetime, os, re
from email.mime.text import MIMEText

SMTP_HOST, SMTP_PORT = "smtp.qq.com", 465
USER = "874600336@qq.com"
PASS = "lsprampnpuqfbdbd"
TO = "874600336@qq.com"
HOST = "2022110014@172.16.115.107"
ROOT = "/storage/public/home/2022110014/hudi/SV_PK_OTHERS/08_SURVIVOR/2_caller/downstream/ml_classification/nested_cv_allcls_20261028"
TOTAL = 5000

PANELS = [10, 20, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900,
          1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000]
METRICS = ["FST", "Delta", "AED", "PIC", "Random"]
CLS = ["RF", "SVM_linear", "SVM_RBF", "NB", "LDA", "KNN", "DecisionTree",
       "XGBoost", "MLP", "Lasso"]
W = dict(RF=6, SVM_linear=1, SVM_RBF=1.5, NB=0.5, LDA=0.2, KNN=0.2,
         DecisionTree=1, XGBoost=2, MLP=30, Lasso=1)

tasks = []
for k in range(1, 6):
    for m in METRICS:
        for n in PANELS:
            for c in CLS:
                tasks.append((k, m, n, c, n * W[c]))
tasks.sort(key=lambda r: -r[4])
TASK_COST = {i + 1: t[4] for i, t in enumerate(tasks)}
TASK_CLS = {i + 1: t[3] for i, t in enumerate(tasks)}
TOTAL_COST = sum(TASK_COST.values())


def ssh(cmd, timeout=300):
    p = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=25", HOST, cmd],
                       capture_output=True, text=True, timeout=timeout)
    return (p.stdout or "") + (p.stderr or "")


def main():
    out = ssh("bash -lc 'bash %s/metrics.sh'" % ROOT)
    states, done_ids, work, errs, nonrun, submitted, maxt = {}, [], 0.0, 0, "", 0, 0.0
    for ln in out.splitlines():
        s = ln.strip()
        if s.startswith("STATE "):
            for cnt, st in re.findall(r"(\d+)\s+(\w+)", s[6:]):
                states[st] = states.get(st, 0) + int(cnt)
        elif s.startswith("DONE "):
            pass
        elif s.startswith("WORK "):
            work = float(s[5:] or 0)
        elif s.startswith("MAXT "):
            maxt = float(s[5:] or 0)
        elif s.startswith("ERRS "):
            errs = int(s[5:] or 0)
        elif s.startswith("NONRUN "):
            nonrun = s[7:].strip()
        elif s.startswith("SUBMITTED "):
            submitted = int(s[10:] or 0)
        elif s.startswith("IDS "):
            done_ids = [int(x) for x in s[4:].split() if x.isdigit()]

    done = len(done_ids)
    done_cost = sum(TASK_COST.get(i, 0) for i in done_ids)
    pct_n = 100.0 * done / TOTAL
    pct_c = 100.0 * done_cost / TOTAL_COST if TOTAL_COST else 0.0
    running, pend = states.get("RUN", 0), states.get("PEND", 0)

    eta_txt = "待估（尚无已完成任务的耗时样本）"
    if done_cost > 0 and work > 0 and running > 0:
        sec_per_cost = work / done_cost
        rem_core = (TOTAL_COST - done_cost) * sec_per_cost
        eta_s = rem_core / running
        # a single MLP@3000 task is a serial chain of 31 nnet fits, so the wall
        # time cannot drop below the cost of the heaviest unfinished task
        heaviest = max(TASK_COST.get(i, 0) for i in TASK_COST if i not in set(done_ids))
        eta_s = max(eta_s, 0.9 * heaviest * sec_per_cost)
        eta_txt = "%.1f h（剩余 %.0f core-h ÷ %d 并发；下限由最慢单任务 %.0f min 决定），预计 %s 完成" % (
            eta_s / 3600.0, rem_core / 3600.0, running, maxt / 60.0,
            (datetime.datetime.now() + datetime.timedelta(seconds=eta_s)).strftime("%m-%d %H:%M"))

    status = ("全部完成" if done >= TOTAL and running == 0 else
              "运行中" if running > 0 else
              "排队中" if pend > 0 else "无运行作业（需检查）")

    by_cls = {}
    for i in done_ids:
        c = TASK_CLS.get(i)
        by_cls[c] = by_cls.get(c, 0) + 1
    cls_line = "  ".join("%s %d/%d" % (c, by_cls.get(c, 0), 500) for c in CLS)

    lines = [
        "SV_ml_nested（172 集群 / 10 分类器全调参）进度报告",
        "生成时间: %s" % datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "",
        "状态        : %s" % status,
        "组合进度    : %d/%d (%.1f%%)" % (done, TOTAL, pct_n),
        "工作量进度  : %.1f%%（按预估耗时加权；剩余主要是 MLP）" % pct_c,
        "作业        : RUN=%d  PEND=%d  已提交=%d/5000" % (running, pend, submitted),
        "剩余 ETA    : %s" % eta_txt,
        "已耗 CPU    : %.0f core-h" % (work / 3600.0),
        "失败组合    : %d" % errs,
        "",
        "各分类器完成: %s" % cls_line,
    ]
    if nonrun:
        lines.append("非运行作业  : %s" % nonrun)
    lines += ["", "目录: %s" % ROOT]
    body = "\n".join(lines)
    subject = "[SV_ml_nested@172] %d/5000 组合 / %.0f%% 工作量 | %s | RUN=%d" % (done, pct_c, status, running)

    if os.environ.get("DRY"):
        print(subject)
        print(body)
        return
    msg = MIMEText(body, "plain", "utf-8")
    msg["Subject"] = subject
    msg["From"] = USER
    msg["To"] = TO
    with smtplib.SMTP_SSL(SMTP_HOST, SMTP_PORT, timeout=30) as s:
        s.login(USER, PASS)
        s.sendmail(USER, [TO], msg.as_string())
    print("sent:", subject)


if __name__ == "__main__":
    main()
