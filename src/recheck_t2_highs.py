"""Second-solver recheck of split_t2.py verdicts with HiGHS (via scipy.optimize.milp).

Model built from scratch here (sparse matrix, words by itertools.product); minimises the code
size subject to covering and the block counts, so 'infeasible' is reported directly by HiGHS.
Usage: python src/recheck_t2_highs.py B M [--tlim 1800]
"""
import argparse, itertools, json, os, sys, time
import numpy as np
from scipy.optimize import LinearConstraint, milp
from scipy.sparse import lil_matrix
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from split_t2 import CELLS, distributions
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser(); ap.add_argument("B", type=int); ap.add_argument("M", type=int)
ap.add_argument("--tlim", type=float, default=1800.0); a = ap.parse_args()
words = list(itertools.product(range(3), range(3), *[range(2)] * a.B)); N = len(words)
cover = lil_matrix((N, N))
for i, w in enumerate(words):
    for j, u in enumerate(words):
        if sum(p != q for p, q in zip(w, u)) <= 1: cover[i, j] = 1
block = lil_matrix((9, N))
for j, w in enumerate(words): block[CELLS.index(w[:2]), j] = 1
dists = distributions(a.B, a.M)
print(f"HiGHS recheck K(2,{a.B},1) M={a.M}: {len(dists)} cases", flush=True)
tally = {}
for k, d in enumerate(dists):
    cons = [LinearConstraint(cover.tocsr(), 1, np.inf), LinearConstraint(block.tocsr(), d, d)]
    t = time.time()
    r = milp(np.ones(N), constraints=cons, integrality=np.ones(N), bounds=(0, 1),
             options={"time_limit": a.tlim, "disp": False})
    st = {2: "infeasible", 0: "feasible", 1: "timeout"}.get(r.status, f"status{r.status}")
    if r.status == 0: st = "feasible"
    tally[st] = tally.get(st, 0) + 1
    print(f"  {k}: {d} -> {st} ({time.time() - t:.0f}s) {r.message[:60]}", flush=True)
print("done:", tally, flush=True)
with open(os.path.join(ROOT, "results", "recheck_t2_highs.jsonl"), "a") as fh:
    fh.write(json.dumps({"B": a.B, "M": a.M, "cases": len(dists), "tally": tally, "tlim": a.tlim}) + "\n")
