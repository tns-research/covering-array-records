#!/usr/bin/env python3
"""CA(4,k,6) for k = A..B, one column at a time.
For each k: extend the best array for k-1 by one column (extend.c, then greedy rows), then shrink it with the C generator of
run 20261103-230206-ca-4-34-6 recompiled with K=k. Each array goes to out/CA_N_4_k_6.txt as soon as it exists.
usage: chain.py start.txt k_start k_end ext_seconds gen_seconds"""
import glob, json, os, re, subprocess, sys, time
H = os.path.dirname(os.path.abspath(__file__)); W = os.path.join(H, "work"); OUT = os.path.join(H, "out")
os.makedirs(W, exist_ok=True); os.makedirs(OUT, exist_ok=True)
start, a, b, es, gs = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
subprocess.run(["gcc", "-O3", "-o", os.path.join(W, "extend"), os.path.join(H, "extend.c")], check=True)
cur = os.path.join(H, start); src = open(os.path.join(H, "gen.c")).read()
rows = lambda p: sum(1 for l in open(p) if l.strip())
for k in range(a, b + 1):
    t0 = time.time(); ext = os.path.join(W, f"ext_{k}.txt")
    with open(ext, "w") as o:
        subprocess.run([os.path.join(W, "extend"), cur, str(k - 1), "6", str(es), str(k)], stdout=o, stderr=open(os.path.join(W, f"ext_{k}.err"), "w"))
    n_ext = rows(ext); gd = os.path.join(W, f"gen_{k}"); os.makedirs(gd, exist_ok=True)
    g = re.sub(r"#define K \d+", f"#define K {k}", src, count=1)
    g = re.sub(r"#define MAXN \d+", f"#define MAXN {2 * n_ext + 1000}", g, count=1)
    open(os.path.join(gd, "gen.c"), "w").write(g)
    subprocess.run(["gcc", "-O3", "-std=gnu11", "-o", "gen", "gen.c", "-lm"], cwd=gd, check=True)
    with open(os.path.join(gd, "cand.txt"), "w") as o:
        p = subprocess.Popen(["./gen"], cwd=gd, stdout=o, stderr=open(os.path.join(gd, "err.txt"), "w"), env={"PATH": "/usr/bin:/bin", "HOME": gd, "BEST_CANDIDATE": ext})
        try:
            p.wait(timeout=gs)
        except subprocess.TimeoutExpired:
            p.terminate()
            try:
                p.wait(timeout=60)
            except subprocess.TimeoutExpired:
                p.kill(); p.wait()
    cand = os.path.join(gd, "cand.txt")
    best = cand if os.path.getsize(cand) and rows(cand) <= n_ext and all(len(l.split()) == k for l in open(cand) if l.strip()) else ext
    n = rows(best)
    open(os.path.join(OUT, f"CA_{n}_4_{k}_6.txt"), "w").write(open(best).read())
    print(json.dumps({"k": k, "extended": n_ext, "N": n, "seconds": round(time.time() - t0)}), flush=True)
    cur = best
print(json.dumps({"done": [a, b]}))
