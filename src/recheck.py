"""Independent recheck of the symmetry exclusions in results/classify.jsonl.

For every class that src/classify.py reported as excluded, this script rebuilds the problem
from scratch and solves it again. It shares no code with search.py, orbit_ilp.py or
classify.py; the only thing taken from the first pass is the generator (perm, mult, shift)
of each class. Differences from the first pass, on purpose:

  * words are tuples from itertools.product, distances are computed by direct comparison;
  * the generator g is replaced by a random conjugate h g h^-1 (seeded), so a different
    but equivalent instance is solved;
  * the solver is selectable: OR-tools CP-SAT (integer arithmetic, capped feasibility),
    CBC (a second LP-based solver, run as a program on an LP file), or HiGHS again.

Output: results/recheck.jsonl, one row per class and solver. Decided rows are skipped on rerun.
"""
import argparse
import itertools
import json
import math
import os
import random
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def make_map(radix, perm, mult, shift):
    return lambda w: tuple((mult[i] * w[perm[i]] + shift[i]) % radix[i] for i in range(len(radix)))


def random_isometry(rng, t, b, radix):
    pt = list(range(t)); rng.shuffle(pt)
    pb = list(range(t, t + b)); rng.shuffle(pb)
    mult = [rng.choice((1, 2)) for _ in range(t)] + [1] * b
    shift = [rng.randrange(radix[i]) for i in range(t + b)]
    return make_map(radix, pt + pb, mult, shift)


def build(t, b, R, spec, seed):
    """Orbits of a random conjugate of the generator, and the cover relation between them."""
    radix = [3] * t + [2] * b
    words = list(itertools.product(*[range(q) for q in radix]))
    g = make_map(radix, *spec)
    h = random_isometry(random.Random(seed), t, b, radix)
    h_of = {w: h(w) for w in words}
    assert len(set(h_of.values())) == len(words)
    h_inv = {v: k for k, v in h_of.items()}
    conj = {w: h_of[g(h_inv[w])] for w in words}        # h g h^-1
    assert len(set(conj.values())) == len(words)
    orbit_of, orbs = {}, []
    for w in words:
        if w in orbit_of:
            continue
        o, x = [], w
        while x not in orbit_of:
            orbit_of[x] = len(orbs)
            o.append(x)
            x = conj[x]
        orbs.append(o)
    rows = []
    for o in orbs:
        rep = o[0]
        near = {orbit_of[w] for w in words if sum(1 for u, v in zip(rep, w) if u != v) <= R}
        rows.append(sorted(near))
    return [len(o) for o in orbs], rows


def solve_cpsat(sizes, rows, ub, tlim, workers):
    from ortools.sat.python import cp_model
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f"x{j}") for j in range(len(sizes))]
    for r in rows:
        m.AddBoolOr([x[j] for j in r])
    m.Add(sum(s * v for s, v in zip(sizes, x)) <= ub)
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = tlim
    s.parameters.num_workers = workers
    st = s.Solve(m)
    return {cp_model.INFEASIBLE: "infeasible", cp_model.OPTIMAL: "feasible",
            cp_model.FEASIBLE: "feasible"}.get(st, "timeout"), None


def solve_cbc(sizes, rows, ub, tlim, workers, cutoff=False):
    # CBC is run as a program on a plain LP file (CBC from the cbcbox package, which is what
    # `pip install pulp[cbc]` brings). PuLP 4 itself is not used: its API changed completely,
    # and calling CBC directly keeps one less layer between the model and the second solver.
    #
    # By default CBC minimises outright and reports the proven optimum, which can be compared
    # with HiGHS's. cutoff=True instead discards every branch that cannot reach a size <= ub
    # (cutoff ub + 0.5; sizes are integers). On two K(1,7,1) classes that was slower (80 s
    # against 28 s, and 300 s+ against 226 s), so it is not the default.
    import subprocess
    import tempfile
    import cbcbox
    with tempfile.TemporaryDirectory() as d:
        lp, sol = os.path.join(d, "cover.lp"), os.path.join(d, "cover.sol")
        with open(lp, "w") as f:
            terms = [f"{s} x{j}" for j, s in enumerate(sizes)]      # 10 terms a line: short lines
            f.write("Minimize\n obj: " + "\n + ".join(" + ".join(terms[i:i + 10])
                                                      for i in range(0, len(terms), 10)) + "\n")
            f.write("Subject To\n")
            for i, r in enumerate(rows):
                f.write(f" c{i}: " + " + ".join(f"x{j}" for j in r) + " >= 1\n")
            f.write("Binary\n " + " ".join(f"x{j}" for j in range(len(sizes))) + "\nEnd\n")
        cmd = [cbcbox.cbc_bin_path(), lp, "-sec", str(tlim), "-threads", str(workers)]
        if cutoff:
            cmd += ["-cutoff", str(ub + 0.5)]
        subprocess.run(cmd + ["-solve", "-solution", sol], capture_output=True, check=False)
        if not os.path.exists(sol):
            return "timeout", None
        with open(sol) as f:
            head = f.readline()
    if cutoff and (head.startswith("Infeasible") or head.startswith("Integer infeasible")):
        return "infeasible", "cutoff"
    # "Optimal - objective value 48.00000000" only when CBC proved optimality.
    if not head.startswith("Optimal"):
        return "timeout", head.strip()
    opt = round(float(head.rsplit(None, 1)[-1]))
    return ("infeasible" if opt > ub else "feasible"), opt


def solve_highs(sizes, rows, ub, tlim, workers):
    import numpy as np
    from scipy.optimize import LinearConstraint, milp
    from scipy.sparse import lil_matrix
    A = lil_matrix((len(rows), len(sizes)))
    for i, r in enumerate(rows):
        A[i, r] = 1
    res = milp(c=np.array(sizes, float), constraints=[LinearConstraint(A.tocsr(), 1, np.inf)],
               integrality=np.ones(len(sizes)), bounds=(0, 1), options={"time_limit": tlim})
    bound = getattr(res, "mip_dual_bound", None)
    if res.x is not None and round(float(res.fun)) <= ub:
        return "feasible", round(float(res.fun))
    if bound is not None and math.ceil(bound - 1e-6) > ub:
        return "infeasible", bound
    return "timeout", bound


def solve_scip(sizes, rows, ub, tlim, workers):
    # SCIP (pip install pyscipopt): minimise, as with HiGHS. If time runs out, a proven lower
    # bound above ub still settles the class. SCIP's search is single-threaded; workers unused.
    from pyscipopt import Model, quicksum
    m = Model()
    m.hideOutput()
    x = [m.addVar(vtype="B", name=f"x{j}") for j in range(len(sizes))]
    m.setObjective(quicksum(s * v for s, v in zip(sizes, x)), "minimize")
    for r in rows:
        m.addCons(quicksum(x[j] for j in r) >= 1)
    m.setParam("limits/time", tlim)
    m.optimize()
    status, bound = m.getStatus(), m.getDualbound()
    if status == "optimal":
        opt = round(m.getObjVal())
        return ("infeasible" if opt > ub else "feasible"), opt
    if m.getNSols() > 0 and round(m.getObjVal()) <= ub:
        return "feasible", round(m.getObjVal())
    if math.ceil(bound - 1e-6) > ub:
        return "infeasible", bound
    return "timeout", bound


SOLVERS = {"cpsat": solve_cpsat, "cbc": solve_cbc, "highs": solve_highs, "scip": solve_scip}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("t", type=int)
    ap.add_argument("b", type=int)
    ap.add_argument("R", type=int)
    ap.add_argument("--record", type=int, required=True)
    ap.add_argument("--solver", choices=sorted(SOLVERS), required=True)
    ap.add_argument("--tlim", type=float, default=300.0)
    ap.add_argument("--time", type=float, default=3600.0)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--only", default="", help="semicolon list of class names (the names contain commas)")
    ap.add_argument("--out", default="recheck.jsonl",
                    help="results file under results/ (one per parallel worker, merged afterwards)")
    a = ap.parse_args()

    name, ub = f"K({a.t},{a.b},{a.R})", a.record - 1
    first = {}
    with open(os.path.join(ROOT, "results", "classify.jsonl")) as f:
        for ln in f:
            r = json.loads(ln)
            if r["target"] == name and r["ub"] == ub and r["status"] == "infeasible":
                first[r["class"]] = r
    out_path = os.path.join(ROOT, "results", a.out)
    done = set()
    if os.path.exists(out_path):
        with open(out_path) as f:
            for ln in f:
                r = json.loads(ln)
                if r["target"] == name and r["ub"] == ub and r["solver"] == a.solver and r["status"] != "timeout":
                    done.add(r["class"])
    todo = sorted((c for c in first if c not in done and (not a.only or c in a.only.split(";"))),
                  key=lambda c: first[c]["orbits"])
    print(f"{name} ub={ub} solver={a.solver}: {len(first)} excluded in first pass, "
          f"{len(done)} already rechecked, {len(todo)} to do", flush=True)
    t0 = time.time()
    agree = 0
    for cname in todo:
        if time.time() - t0 > a.time:
            break
        sizes, rows = build(a.t, a.b, a.R, first[cname]["gen"], a.seed)
        assert len(sizes) == first[cname]["orbits"], (cname, len(sizes), first[cname]["orbits"])
        ts = time.time()
        status, info = SOLVERS[a.solver](sizes, rows, ub, a.tlim, a.workers)
        agree += status == "infeasible"
        row = {"target": name, "ub": ub, "class": cname, "solver": a.solver, "seed": a.seed,
               "orbits": len(sizes), "status": status, "info": info, "solve_s": round(time.time() - ts, 2)}
        with open(out_path, "a") as f:
            f.write(json.dumps(row) + "\n")
        flag = "" if status == "infeasible" else "   <-- does NOT confirm"
        print(f"  {cname:22s} orbits={len(sizes):4d} {status:10s} ({row['solve_s']}s){flag}", flush=True)
    print(f"done: {agree} confirmed this run", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
