#!/usr/bin/env python3
"""Offline evaluation of the relative-stagnation stop rule for PALNS.

Reads [palns-shared-trace] and [palns] lines from a reference run made with the
most patient setting (gamma_max, Nmin_max). For any gamma <= gamma_max and
Nmin <= Nmin_max the rule stops no later than the reference run, so its stop
time and cost can be read from the trace.

Stop rule (Phase II): stop when t - t_last >= gamma * t and the thread has done
Nmin iterations since the last shared improvement, or (optional cap) when the
thread has done Nmax iterations since the last shared improvement. Nmin is converted to time
with the Phase II iteration rate per thread r: Nmin / r seconds.

usage: stag_sim.py REF_DIR BASELINE_DIR OUT_DIR
"""
import sys, os, re, glob, math, csv, statistics

def read(p):
    with open(p, errors="ignore") as f:
        return f.read()

def baseline(bdir):
    out = {}
    for f in glob.glob(os.path.join(bdir, "RecDep_day_*_p0_log.txt")):
        s = read(f); nm = os.path.basename(f)[11:-11]
        c = re.search(r"\[initial-incumbent\] original_cost=([0-9.]+)", s)
        t = re.search(r"\[initial-incumbent-attempt\] index=0 .*? runtime_seconds=([0-9.]+)", s)
        if c and t:
            out[nm] = (float(c.group(1)), float(t.group(1)))
    return out

def traces(rdir):
    out = {}
    for f in glob.glob(os.path.join(rdir, "RecDep_day_*_p0_log.txt")):
        nm = os.path.basename(f)[11:-11]; s = read(f)
        m = re.search(r"^\[palns\] k=(\d+) .*runtime_seconds=([0-9.e+-]+).*time_to_feasible=([0-9.e+-]+).*phase2_iterations=(\d+) threads", s, re.M)
        tr = re.search(r"^\[palns-shared-trace\] .*points=(.*)$", s, re.M)
        if not m or not tr:
            continue
        pts = [tuple(map(float, p.split(":"))) for p in tr.group(1).split(";") if p]
        rt = float(m.group(2)); ttf = float(m.group(3)); it2 = int(m.group(4))
        rate = it2 / 8.0 / max(rt - ttf, 1e-6) if it2 > 0 else 1e9
        out[nm] = dict(K=int(m.group(1)), pts=pts, rt=rt, ttf=ttf, rate=rate)
    return out

def simulate(x, g, nmin, tl=1200.0, nmax=None):
    pts = x["pts"]
    for i, (t, c) in enumerate(pts):
        s = max(t + nmin / x["rate"], t / (1.0 - g))
        if nmax:
            s = min(s, t + nmax / x["rate"])  # absolute cap: stop after nmax stagnant iterations
        nxt = pts[i + 1][0] if i + 1 < len(pts) else math.inf
        if s < nxt:
            return min(s, tl, x["rt"]), c
    return x["rt"], pts[-1][1]

def first_hit(x, target):
    for t, c in x["pts"]:
        if c <= target + 1e-6:
            return t
    return math.inf

def key(nm):
    return (nm[0], int(nm[1:]))

def main():
    rdir, bdir, odir = sys.argv[1:4]
    os.makedirs(odir, exist_ok=True)
    B = baseline(bdir); T = traces(rdir)
    gs = [round(0.05 * k, 2) for k in range(1, 17)]
    ns = [500, 1000, 2000, 3000, 5000, 7500, 10000, 15000, 20000, 30000, 50000]
    nmaxs = [None, 50000, 75000, 100000, 150000, 200000]
    rows, summ = [], []
    for g in gs:
      for nmax in nmaxs:
        for n in ns:
            if nmax and nmax <= n: continue
            fails = []; tsum = 0.0; tall = 0.0; ratios = []; passt = 0; worst_ratio = 0.0
            worst_excess = -math.inf; margins = []
            for nm in sorted(T, key=key):
                x = T[nm]; t, c = simulate(x, g, n, nmax=nmax); tall += t
                rec = dict(gamma=g, nmin=n, nmax=nmax or "", inst=nm, K=x["K"], stop_s=round(t, 3), cost=c)
                if nm in B:
                    bc, bt = B[nm]; tsum += t
                    hit = first_hit(x, bc)
                    margin = t - hit if c <= bc + 1e-6 else -1.0
                    rec.update(base_cost=bc, base_s=bt, gap_pct=round((c - bc) / bc * 100, 4),
                               time_ratio=round(t / bt, 4), time_ok=int(t <= max(1.2 * bt, 1.0)),
                               hit_base_s=round(hit, 3) if math.isfinite(hit) else "", margin_s=round(margin, 3))
                    if c > bc + 1e-6: fails.append(nm)
                    ratios.append(max(t, 1e-3) / bt); passt += t <= max(1.2 * bt, 1.0)
                    worst_ratio = max(worst_ratio, t / bt); worst_excess = max(worst_excess, t - bt)
                    if c <= bc + 1e-6: margins.append((margin, nm))
                rows.append(rec)
            mm = min(margins) if margins else (math.nan, "")
            summ.append(dict(gamma=g, nmin=n, nmax=nmax or "", fails=len(fails), fail_list=" ".join(fails),
                             total_s_62=round(tsum, 1), total_s_65=round(tall, 1),
                             median_ratio=round(statistics.median(ratios), 3),
                             geomean_ratio=round(math.exp(sum(math.log(r) for r in ratios) / len(ratios)), 3),
                             time_ok=passt, worst_ratio=round(worst_ratio, 2), worst_excess_s=round(worst_excess, 1),
                             min_margin_s=round(mm[0], 2), min_margin_inst=mm[1]))
    with open(os.path.join(odir, "per_instance.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=sorted({k for r in rows for k in r}, key=lambda k: ["gamma","nmin","nmax","inst","K","base_cost","base_s","cost","gap_pct","stop_s","time_ratio","time_ok","hit_base_s","margin_s"].index(k) if k in ["gamma","nmin","nmax","inst","K","base_cost","base_s","cost","gap_pct","stop_s","time_ratio","time_ok","hit_base_s","margin_s"] else 99))
        w.writeheader(); w.writerows(rows)
    with open(os.path.join(odir, "summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summ[0].keys())); w.writeheader(); w.writerows(summ)
    print("instances:", len(T), "combos:", len(summ), "feasible:", sum(s["fails"] == 0 for s in summ))

if __name__ == "__main__":
    main()
