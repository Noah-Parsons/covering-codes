"""Independent CP-SAT recheck of split_t2.py's 'no code' verdicts (integer arithmetic).

Takes the reduced block distributions from split_t2.distributions (the mapping check showed
every raw distribution is equivalent to one of them) and re-solves each with OR-tools CP-SAT,
built here from scratch: words by itertools.product, balls by direct distance.
Usage: python src/recheck_t2_cpsat.py B M [--tlim 1800]
"""
import argparse, itertools, json, os, sys, time
from ortools.sat.python import cp_model
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from split_t2 import CELLS, distributions
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ap = argparse.ArgumentParser(); ap.add_argument("B", type=int); ap.add_argument("M", type=int)
ap.add_argument("--tlim", type=float, default=1800.0); a = ap.parse_args()
words = list(itertools.product(range(3), range(3), *[range(2)] * a.B))
idx = {w: i for i, w in enumerate(words)}
dist = lambda u, v: sum(p != q for p, q in zip(u, v))
balls = [[idx[u] for u in words if dist(w, u) <= 1] for w in words]
dists = distributions(a.B, a.M)
print(f"CP-SAT recheck K(2,{a.B},1) M={a.M}: {len(dists)} cases", flush=True)
tally = {}
for k, d in enumerate(dists):
    md = cp_model.CpModel(); x = [md.NewBoolVar("") for _ in words]
    for bl in balls: md.AddBoolOr([x[j] for j in bl])
    for ci, c in enumerate(CELLS): md.Add(sum(x[idx[w]] for w in words if w[:2] == c) == d[ci])
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = a.tlim; s.parameters.num_workers = 1
    t = time.time(); st = s.StatusName(s.Solve(md)); tally[st] = tally.get(st, 0) + 1
    print(f"  {k}: {d} -> {st} ({time.time() - t:.0f}s)", flush=True)
print("done:", tally, flush=True)
with open(os.path.join(ROOT, "results", "recheck_t2_cpsat.jsonl"), "a") as fh:
    fh.write(json.dumps({"B": a.B, "M": a.M, "cases": len(dists), "tally": tally, "tlim": a.tlim}) + "\n")
