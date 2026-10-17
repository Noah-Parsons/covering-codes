"""Third-solver recheck of split_t2.py verdicts with Gurobi (gurobipy, size-limited licence).

Model built here from scratch: one binary per word, covering constraint per word (balls by
direct Hamming distance), fixed block counts. Single thread. Reports Gurobi's status per case.
Usage: python src/recheck_t2_gurobi.py B M [--tlim 600]
"""
import argparse, itertools, json, os, sys, time
import gurobipy as gp
from gurobipy import GRB
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from split_t2 import CELLS, distributions
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ap = argparse.ArgumentParser(); ap.add_argument("B", type=int); ap.add_argument("M", type=int)
ap.add_argument("--tlim", type=float, default=600.0); a = ap.parse_args()
words = list(itertools.product(range(3), range(3), *[range(2)] * a.B))
balls = [[j for j, u in enumerate(words) if sum(p != q for p, q in zip(w, u)) <= 1] for w in words]
names = {GRB.INFEASIBLE: "infeasible", GRB.OPTIMAL: "feasible", GRB.TIME_LIMIT: "timeout"}
env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
dists = distributions(a.B, a.M)
print(f"Gurobi {'.'.join(map(str, gp.gurobi.version()))} recheck K(2,{a.B},1) M={a.M}: {len(dists)} cases",
      flush=True)
tally = {}
for k, d in enumerate(dists):
    m = gp.Model(env=env)
    m.Params.TimeLimit = a.tlim; m.Params.Threads = 1
    x = m.addVars(len(words), vtype=GRB.BINARY)
    for bl in balls: m.addConstr(gp.quicksum(x[j] for j in bl) >= 1)
    for ci, c in enumerate(CELLS):
        m.addConstr(gp.quicksum(x[j] for j, w in enumerate(words) if w[:2] == c) == d[ci])
    t = time.time(); m.optimize()
    st = names.get(m.Status, f"status{m.Status}")
    tally[st] = tally.get(st, 0) + 1
    print(f"  {k}: {d} -> {st} ({time.time() - t:.1f}s)", flush=True)
print("done:", tally, flush=True)
with open(os.path.join(ROOT, "results", "recheck_t2_gurobi.jsonl"), "a") as fh:
    fh.write(json.dumps({"B": a.B, "M": a.M, "cases": len(dists), "tally": tally, "tlim": a.tlim}) + "\n")
