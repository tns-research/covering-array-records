#!/usr/bin/env python3
"""Round 13 C search driver: compile starter.c, run it on one (family, k), develop and check any witness exactly.
usage: run.py FAMILY K SECONDS [SEED]   (C3: f=0, C4: f=1, C5: f=2). A miss decides nothing (search only)."""
import json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common import develop, ca_missing, LIVE, LAB, save_rows  # noqa: E402

fam, k, secs = sys.argv[1], int(sys.argv[2]), sys.argv[3]
seed = sys.argv[4] if len(sys.argv) > 4 else "1"
f = {"C3": 0, "C4": 1, "C5": 2}[fam]
exe = os.path.join("/tmp", f"starter13-{os.getpid()}")
subprocess.run(["gcc", "-O3", "-march=native", "-o", exe, os.path.join(HERE, "starter.c"), "-lm"], check=True)
p = subprocess.run([exe, str(k), str(f), secs, seed], capture_output=True, text=True)
os.unlink(exe)
lines = [json.loads(x) for x in p.stdout.splitlines() if x.startswith("{")]
best = lines[-1] if lines else {}
bar = min(LIVE.get(k, 10 ** 9), LAB.get(k, 10 ** 9))
out = {"family": fam, "k": k, "f": f, "value": 270 * (k - f) + 6, "bar": bar, "best_holes": best.get("holes"),
       "search": p.stderr.strip()}
if best.get("holes") == 0:
    rows = develop(best["starters"], k, f)
    miss, ex = ca_missing(rows)
    out.update(N=int(rows.shape[0]), missing=miss, starters=best["starters"])
    if miss == 0:
        path = os.path.join(HERE, "out", f"{fam}_CA_{rows.shape[0]}_4_{k}_6.txt")
        save_rows(path, rows.tolist())
        json.dump({"k": k, "f": f, "starters": best["starters"], "N": int(rows.shape[0])},
                  open(path.replace(".txt", "_starters.json"), "w"))
        out["file"] = os.path.basename(path)
        out["below_bar"] = rows.shape[0] < bar
    else:
        out["error"] = f"orbit reduction mismatch {ex}"
print(json.dumps(out), flush=True)
