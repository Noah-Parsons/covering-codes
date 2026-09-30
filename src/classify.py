"""Which symmetries can a record-beating covering code have?

Claim being tested, for a target K(t, b, R) with published upper bound `record`:

    no code of size <= record - 1 has a non-trivial automorphism.

A non-trivial automorphism group contains an element of prime order, and a code invariant
under g is invariant under every conjugate of g up to equivalence of codes. So it is enough
to take one representative g of every conjugacy class of prime-order elements of the full
isometry group  (S3 wr S_t) x (S2 wr S_b),  and to show that no union of <g>-orbits of total
size <= record - 1 covers the space. Each case is the weighted set cover of src/orbit_ilp.py
with the size cap added; "infeasible" from the solver settles that class.

Conjugacy classes of prime order p in a wreath product S_k wr S_n: the coordinate
permutation has cycles of length 1 and p; each p-cycle can be conjugated to act with no
symbol change; each fixed coordinate carries a symbol permutation of order 1 or p, and only
its conjugacy class in S_k matters. That gives the representatives built in `prime_classes`.
`tests/test_classify.py` checks the count against brute-force enumeration on small spaces.

Results go to results/classify.jsonl; decided classes are skipped on rerun.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from orbit_ilp import (closure, gen_perm, orbit_problem, orbits, solve,  # noqa: E402
                       solve_cpsat, solve_min, verify_and_save)
from search import all_digits, build_balls  # noqa: E402


def primes_upto(n):
    return [p for p in range(2, n + 1) if all(p % q for q in range(2, p))]


def prime_classes(t, b):
    """One (p, name, spec) per conjugacy class of prime-order isometries. spec = (perm, mult, shift)."""
    out = []
    n = t + b
    for p in primes_upto(max(t, b, 3)):
        for a in range(t // p + 1):                 # p-cycles on ternary coordinates
            jmax = t - p * a if p in (2, 3) else 0  # fixed ternary coords with a symbol map of order p
            for j in range(jmax + 1):
                for c in range(b // p + 1):         # p-cycles on binary coordinates
                    kmax = b - p * c if p == 2 else 0   # fixed binary coords that are flipped
                    for k in range(kmax + 1):
                        if a == j == c == k == 0:
                            continue
                        perm = list(range(n))
                        mult = [1] * n
                        shift = [0] * n
                        for i in range(a):
                            cyc = list(range(p * i, p * i + p))
                            for u, v in zip(cyc, cyc[1:] + cyc[:1]):
                                perm[u] = v
                        for i in range(p * a, p * a + j):
                            if p == 2:
                                mult[i] = 2      # x -> 2x: the transposition (1 2)
                            else:
                                shift[i] = 1     # x -> x+1: a 3-cycle
                        for i in range(c):
                            cyc = list(range(t + p * i, t + p * i + p))
                            for u, v in zip(cyc, cyc[1:] + cyc[:1]):
                                perm[u] = v
                        for i in range(t + p * c, t + p * c + k):
                            shift[i] = 1
                        out.append((p, f"p{p}:t{a}c{j}s,b{c}c{k}f", (perm, mult, shift)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("t", type=int)
    ap.add_argument("b", type=int)
    ap.add_argument("R", type=int)
    ap.add_argument("--record", type=int, required=True)
    ap.add_argument("--tlim", type=float, default=30.0, help="solver seconds per class")
    ap.add_argument("--time", type=float, default=170.0, help="total wall seconds")
    ap.add_argument("--max-orbits", type=int, default=10 ** 9, help="skip classes with more orbits")
    ap.add_argument("--solver", choices=["highs", "cpsat", "highs-min"], default="highs-min")
    ap.add_argument("--workers", type=int, default=4, help="CP-SAT threads")
    ap.add_argument("--recheck", action="store_true",
                    help="re-solve classes that only the OTHER solver has decided")
    ap.add_argument("--primes", default="", help="comma list: only classes of these orders")
    ap.add_argument("--list", action="store_true", help="only list the classes and orbit counts")
    a = ap.parse_args()

    t, b, R = a.t, a.b, a.R
    name = f"K({t},{b},{R})"
    ub = a.record - 1
    digs = all_digits(t, b)
    balls = build_balls(t, b, R)
    log_path = os.path.join(ROOT, "results", "classify.jsonl")
    decided = {}
    if os.path.exists(log_path):
        with open(log_path) as f:
            for ln in f:
                r = json.loads(ln)
                if r["target"] == name and r["ub"] == ub and r["status"] != "timeout":
                    if a.recheck and r.get("solver", "highs") != a.solver:
                        continue
                    decided[r["class"]] = r["status"]

    classes = prime_classes(t, b)
    if a.primes:
        classes = [c for c in classes if c[0] in {int(x) for x in a.primes.split(",")}]
    decided = {k: v for k, v in decided.items() if k in {c[1] for c in classes}}
    work = []
    for p, cname, spec in classes:
        elems = closure([gen_perm(digs, t, b, *spec)], 64)
        assert elems is not None and len(elems) == p, (cname, None if elems is None else len(elems))
        reps, inv, sizes, A = orbit_problem(balls, orbits(elems))
        work.append((len(reps), p, cname, spec, inv, sizes, A))
    work.sort(key=lambda w: w[0])
    print(f"{name}: record={a.record} ub={ub}; {len(classes)} prime-order classes, "
          f"{len(decided)} already decided", flush=True)
    if a.list:
        for norb, p, cname, *_ in work:
            print(f"  {cname:22s} orbits={norb:5d} {decided.get(cname, '')}")
        return 0

    t0 = time.time()
    for norb, p, cname, spec, inv, sizes, A in work:
        if cname in decided or norb > a.max_orbits:
            continue
        if time.time() - t0 > a.time:
            break
        ts = time.time()
        bound = None
        if a.solver == "highs-min":
            status, x, bound = solve_min(A, sizes, ub, a.tlim)
        elif a.solver == "cpsat":
            status, x = solve_cpsat(A, sizes, ub, a.tlim, a.workers)
        else:
            status, x = solve(A, sizes, ub, a.tlim)
        row = {"solver": a.solver, "target": name, "ub": ub, "class": cname, "p": p, "gen": spec, "orbits": norb,
               "status": status, "solve_s": round(time.time() - ts, 2), "tlim": a.tlim, "dual_bound": bound}
        if x is not None:
            words = np.flatnonzero(x[inv] == 1)
            ok, m, nbad, path = verify_and_save(words, t, b, R, f"invariant under {json.dumps(spec)}")
            row.update(M=m, checker="VALID" if ok else f"INVALID({nbad})",
                       file=os.path.relpath(path, ROOT) if path else "")
            print(f"  *** {cname}: code of size {m} <= {ub}, checker {row['checker']} {row['file']}", flush=True)
        if status != "timeout":
            decided[cname] = status
        with open(log_path, "a") as f:
            f.write(json.dumps(row) + "\n")
        print(f"  {cname:22s} orbits={norb:5d} {status:10s} bound={bound} ({row['solve_s']}s)", flush=True)
    left = [w[2] for w in work if w[2] not in decided]
    print(f"done: {len(decided)}/{len(classes)} classes decided, "
          f"{sum(v == 'infeasible' for v in decided.values())} infeasible; undecided: {len(left)}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
