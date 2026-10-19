"""CPLEX (Community Edition) recheck of split_t2.py verdicts for K(2,B,1).

Model built here from scratch: one binary per word, a covering row per word (balls by Hamming
distance), a count row per block. Single thread. Community Edition allows 1000 variables and
1000 constraints, enough for B <= 6 (576 variables, 585 constraints).
Usage: python src/recheck_t2_cplex.py B M [--tlim 600]
"""
import argparse, itertools, json, os, sys, time
import cplex
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from split_t2 import CELLS, distributions
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser(); ap.add_argument("B", type=int); ap.add_argument("M", type=int)
ap.add_argument("--tlim", type=float, default=600.0); a = ap.parse_args()
words = list(itertools.product(range(3), range(3), *[range(2)] * a.B)); N = len(words)
balls = [[j for j, u in enumerate(words) if sum(p != q for p, q in zip(w, u)) <= 1] for w in words]
dists = distributions(a.B, a.M)
print(f"CPLEX {cplex.Cplex().get_version()} recheck K(2,{a.B},1) M={a.M}: {len(dists)} cases", flush=True)
tally = {}
for k, d in enumerate(dists):
    c = cplex.Cplex()
    for s in (c.set_log_stream, c.set_results_stream, c.set_warning_stream, c.set_error_stream): s(None)
    c.parameters.timelimit.set(a.tlim); c.parameters.threads.set(1)
    c.variables.add(types="B" * N)
    rows = [[bl, [1.0] * len(bl)] for bl in balls]
    senses, rhs = "G" * N, [1.0] * N
    for ci, cell in enumerate(CELLS):
        idx = [j for j, w in enumerate(words) if w[:2] == cell]
        rows.append([idx, [1.0] * len(idx)]); senses += "E"; rhs.append(float(d[ci]))
    c.linear_constraints.add(lin_expr=rows, senses=senses, rhs=rhs)
    t = time.time(); c.solve(); S = c.solution.status; s = c.solution.get_status()
    st = ("infeasible" if s in (S.MIP_infeasible, S.MIP_infeasible_or_unbounded) else
          "feasible" if s in (S.MIP_optimal, S.optimal_tolerance, S.MIP_feasible) else "timeout")
    tally[st] = tally.get(st, 0) + 1
    print(f"  {k}: {d} -> {st} ({time.time() - t:.1f}s)", flush=True)
print("done:", tally, flush=True)
with open(os.path.join(ROOT, "results", "recheck_t2_cplex.jsonl"), "a") as fh:
    fh.write(json.dumps({"B": a.B, "M": a.M, "cases": len(dists), "tally": tally, "tlim": a.tlim}) + "\n")
