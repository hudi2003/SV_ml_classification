#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Email a progress report for the all-classifier nested-CV run on mserver.

Reads full_allcls.log + nested_cv_full_detail.csv, estimates the remaining
core-seconds from the combinations already finished, and sends the report to
the address in ~/.mail/config.json (QQ SMTP, SSL 465).
"""
import csv, json, os, re, smtplib, subprocess, statistics, sys
from datetime import datetime
from email.mime.text import MIMEText
from email.header import Header

BASE = "/home/wangxia/hudi/SV_ml_nested"
N_FOLDS, COMBOS_PER_FOLD = 5, 1000
TOTAL = N_FOLDS * COMBOS_PER_FOLD
LOG = os.path.join(BASE, "full_allcls.log")
DETAIL = os.path.join(BASE, "nested_cv_full_detail.csv")


def sh(cmd, timeout=20):
    try:
        return subprocess.run(cmd, shell=True, capture_output=True, text=True,
                              timeout=timeout).stdout.strip()
    except Exception:
        return ""


def read_text(path):
    if not os.path.exists(path):
        return ""
    with open(path, errors="replace") as fh:
        return fh.read()


def n_r_procs():
    return int(sh("ps -u $(id -u) -o comm= | grep -c '^R$'") or 0)


def load_detail():
    if not os.path.exists(DETAIL):
        return []
    try:
        with open(DETAIL, newline="") as fh:
            return list(csv.DictReader(fh))
    except Exception:
        return []


def fmt_dur(sec):
    if sec is None or sec != sec or sec <= 0:
        return "未知"
    sec = int(sec)
    d, rem = divmod(sec, 86400)
    h, rem = divmod(rem, 3600)
    m = rem // 60
    if d:
        return f"{d}d{h:02d}h{m:02d}m"
    if h:
        return f"{h}h{m:02d}m"
    return f"{m}m"


def predict_remaining_core_s(rows, done_keys):
    """Sum estimated core-seconds for the combinations not yet finished.

    Uses the per-classifier median wall time of finished combinations and
    rescales linearly in the panel size (cost ~ O(n_markers)).
    """
    by_cls = {}
    for r in rows:
        try:
            t = float(r["time_sec"])
            n = int(r["n_markers"])
        except (KeyError, ValueError):
            continue
        if t > 0:
            by_cls.setdefault(r["classifier"], []).append((n, t))
    if not by_cls:
        return None, 0
    med = {c: (statistics.median([n for n, _ in v]),
               statistics.median([t for _, t in v])) for c, v in by_cls.items()}
    total, n_missing = 0.0, 0
    for k in range(1, N_FOLDS + 1):
        for metric in ("FST", "Delta", "AED", "PIC", "Random"):
            for n in (10, 20, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900,
                      1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000):
                for cls in ("RF", "SVM_linear", "SVM_RBF", "NB", "LDA", "KNN",
                            "DecisionTree", "XGBoost", "MLP", "Lasso"):
                    key = f"{k}|{metric}|{n}|{cls}"
                    if key in done_keys:
                        continue
                    if cls in med:
                        n0, t0 = med[cls]
                        total += t0 * (n / n0 if n0 else 1.0)
                    else:
                        n_missing += 1
                        total += 120.0 * (n / 500.0)
    return total, n_missing


def main():
    log = read_text(LOG)
    rows = load_detail()
    procs = n_r_procs()
    exit_m = re.findall(r"EXIT=(\d+)", log)
    finished = len(rows) >= TOTAL or (exit_m and exit_m[-1] == "0")

    starts = re.findall(r"(\d{2}:\d{2}:\d{2}) Recomputing marker ranking on the outer training set", log)
    folds_done = sorted({int(r["fold"]) for r in rows}) if rows else []
    cur_fold = max(folds_done) if folds_done else (len(starts) if starts else 1)
    in_fold = sum(1 for r in rows if int(r["fold"]) == cur_fold) if rows else 0

    done_keys = {f"{r['fold']}|{r['metric']}|{r['n_markers']}|{r['classifier']}" for r in rows}
    rem_core, n_missing = predict_remaining_core_s(rows, done_keys)

    now = datetime.now()
    # elapsed of the whole run: earliest log timestamp of today's session
    elapsed = None
    if starts:
        try:
            t0 = datetime.strptime(starts[0], "%H:%M:%S").replace(
                year=now.year, month=now.month, day=now.day)
            if t0 > now:
                t0 = t0.replace(day=now.day - 1)
            elapsed = (now - t0).total_seconds()
        except Exception:
            elapsed = None

    cores = 32
    used_core = sum(float(r["time_sec"]) for r in rows if r.get("time_sec")) if rows else 0.0
    eff = None
    if elapsed and used_core:
        eff = max(0.05, min(1.0, used_core / (cores * elapsed)))

    eta_s = None
    if rem_core is not None and eff:
        eta_s = rem_core / (cores * eff)
    eta_txt = "未知（首个分块尚未落盘）" if eta_s is None else (
        fmt_dur(eta_s) + " → 预计 " + (now + __import__("datetime").timedelta(seconds=eta_s)).strftime("%m-%d %H:%M"))

    if finished:
        status = "已完成"
    elif procs > 0:
        status = f"运行中（{procs} 个 R 进程，{cores} 线程）"
    elif exit_m:
        status = f"已停止（EXIT={exit_m[-1]}）"
    else:
        status = "未检测到 R 进程"

    tail = "\n".join([ln for ln in log.splitlines() if ln.strip()][-6:])
    body = f"""SV_ml_nested 全分类器嵌套 CV 进度报告
生成时间: {now:%Y-%m-%d %H:%M:%S}

状态     : {status}
进度     : {len(rows)}/{TOTAL} 组合（第 {cur_fold}/{N_FOLDS} 折，本折 {in_fold}/{COMBOS_PER_FOLD}）
已用时长 : {fmt_dur(elapsed)}
剩余 ETA : {eta_txt}
并行效率 : {('%.0f%%' % (eff * 100)) if eff else '未知'}
已完成折 : {', '.join(map(str, folds_done)) or '无'}

失败组合 : {sum(1 for r in rows if r.get('status') not in (None, '', 'ok'))}
错误信息 : {'; '.join(sorted({r.get('error','') for r in rows if r.get('status') not in (None,'','ok') and r.get('error')}))[:400] or '无'}

日志末尾:
{tail}

目录: {BASE}
日志: {LOG}
"""
    cfg = json.load(open(os.path.expanduser("~/.mail/config.json")))
    msg = MIMEText(body, "plain", "utf-8")
    subj = f"[SV_ml_nested] {len(rows)}/{TOTAL} | 第{cur_fold}折 | {status} | ETA {fmt_dur(eta_s) if eta_s else '待定'}"
    msg["Subject"] = Header(subj, "utf-8")
    msg["From"] = cfg["from"]
    msg["To"] = cfg["user"]
    with smtplib.SMTP_SSL(cfg["host"], int(cfg["port"]), timeout=30) as s:
        s.login(cfg["user"], cfg["pass"])
        s.sendmail(cfg["from"], [cfg["user"]], msg.as_string())
    print(f"{now:%H:%M:%S} mail sent: {subj}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # keep cron quiet but traceable
        print(f"{datetime.now():%H:%M:%S} ERROR {exc!r}", file=sys.stderr)
        sys.exit(1)
