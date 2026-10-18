"""Solve the leaf distributions of src/exact_dist.py (one ternary coordinate) with Gurobi or SCIP.

exact_dist enumerates codeword distributions over cells (ternary value x first L binary
coordinates) and reduces them by symmetry; here each leaf becomes an exact covering problem with
fixed cell counts. A code of size M exists iff some leaf is feasible.
Usage: python src/leaves_t1.py B M --leaf L --solver gurobi|scip [--tlim 600]
"""
import argparse, itertools, json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exact_dist as E
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ap = argparse.ArgumentParser()
ap.add_argument("B", type=int); ap.add_argument("M", type=int)
ap.add_argument("--leaf", type=int, default=2); ap.add_argument("--solver", default="gurobi")
ap.add_argument("--tlim", type=float, default=600.0)
ap.add_argument("--only", type=int, nargs="*", default=None, help="solve only these leaf indices")
a = ap.parse_args()

dists = E.level0(a.B, a.M)
for m in range(a.leaf):
    G = E.group_perms(m + 1); seen = {}
    for d in dists:
        for c in E.children(d, m, a.B, a.M):
            seen.setdefault(E.canonical(c, G), c)
    dists = list(seen.values())
words = list(itertools.product(range(3), *[range(2)] * a.B))
balls = [[j for j, u in enumerate(words) if sum(p != q for p, q in zip(w, u)) <= 1] for w in words]
n = 1 << a.leaf
def cell(w):
    x = 0
    for bit in w[1:1 + a.leaf]: x = 2 * x + bit
    return w[0] * n + x
cells = [cell(w) for w in words]
print(f"K(1,{a.B},1) M={a.M}: {len(dists)} leaves at level {a.leaf}, solver {a.solver}", flush=True)

def solve(d):
    if a.solver == "gurobi":
        import gurobipy as gp
        from gurobipy import GRB
        env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
        m = gp.Model(env=env); m.Params.TimeLimit = a.tlim; m.Params.Threads = 1
        x = m.addVars(len(words), vtype=GRB.BINARY)
        for bl in balls: m.addConstr(gp.quicksum(x[j] for j in bl) >= 1)
        for c in range(3 * n): m.addConstr(gp.quicksum(x[j] for j in range(len(words)) if cells[j] == c) == int(d[c]))
        m.optimize()
        st = {GRB.INFEASIBLE: "infeasible", GRB.OPTIMAL: "feasible", GRB.TIME_LIMIT: "timeout"}.get(m.Status, str(m.Status))
        code = [words[j] for j in range(len(words)) if st == "feasible" and x[j].X > 0.5]
        return st, code
    if a.solver == "cplex":
        # CPLEX Community Edition (pip install cplex): limit 1000 variables / 1000 constraints.
        import cplex
        c = cplex.Cplex()
        for stream in (c.set_log_stream, c.set_results_stream, c.set_warning_stream, c.set_error_stream):
            stream(None)
        c.parameters.timelimit.set(a.tlim); c.parameters.threads.set(1)
        c.variables.add(types="B" * len(words))
        rows, senses, rhs = [], [], []
        for bl in balls:
            rows.append([bl, [1.0] * len(bl)]); senses.append("G"); rhs.append(1.0)
        for cc in range(3 * n):
            idx = [j for j in range(len(words)) if cells[j] == cc]
            rows.append([idx, [1.0] * len(idx)]); senses.append("E"); rhs.append(float(d[cc]))
        c.linear_constraints.add(lin_expr=rows, senses="".join(senses), rhs=rhs)
        c.solve()
        s = c.solution.get_status()
        S = c.solution.status
        if s in (S.MIP_infeasible, S.MIP_infeasible_or_unbounded):
            return "infeasible", []
        if s in (S.MIP_optimal, S.optimal_tolerance, S.MIP_feasible):
            vals = c.solution.get_values()
            return "feasible", [words[j] for j in range(len(words)) if vals[j] > 0.5]
        return "timeout", []
    from pyscipopt import Model, quicksum
    m = Model(); m.hideOutput(); m.setParam("limits/time", a.tlim)
    x = [m.addVar(vtype="B") for _ in words]
    for bl in balls: m.addCons(quicksum(x[j] for j in bl) >= 1)
    for c in range(3 * n): m.addCons(quicksum(x[j] for j in range(len(words)) if cells[j] == c) == int(d[c]))
    m.optimize(); st = m.getStatus()
    code = [words[j] for j in range(len(words)) if m.getNSols() and m.getVal(x[j]) > 0.5]
    return ("feasible" if code else ("infeasible" if st == "infeasible" else "timeout")), code

tally = {}; t0 = time.time()
for k, d in enumerate(dists):
    if a.only is not None and k not in a.only:
        continue
    t = time.time(); st, code = solve(d); tally[st] = tally.get(st, 0) + 1
    print(f"  {k}: {tuple(int(v) for v in d)} -> {st} ({time.time() - t:.1f}s)", flush=True)
    if code:
        fn = os.path.join(ROOT, "results", "codes", f"K_1_{a.B}_1_M{len(code)}_leaves.txt")
        open(fn, "w").write("\n".join("".join(map(str, w)) for w in code) + "\n")
        print(f"  CODE FOUND -> {fn}", flush=True); break
print("done:", tally, f"({time.time() - t0:.0f}s)", flush=True)
with open(os.path.join(ROOT, "results", "leaves_t1.jsonl"), "a") as fh:
    fh.write(json.dumps({"B": a.B, "M": a.M, "leaf": a.leaf, "solver": a.solver, "leaves": len(dists),
                         "tally": tally, "tlim": a.tlim}) + "\n")
