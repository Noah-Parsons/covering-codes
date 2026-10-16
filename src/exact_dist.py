"""Exact search for covering codes of radius 1 in F_3 x F_2^b (the Östergård–Blass method).

Question: does a covering code of radius 1 with exactly M words exist? (A smaller code can be
padded to M words, so "no code of size M" proves K(1,b,1) >= M+1.)

Method (Östergård and Blass, IEEE Trans. Inform. Theory 47 (2001), 2556-2557), adapted to one
ternary coordinate:

* Level m (m = 0..b): fix the ternary coordinate and the first m binary coordinates. The space
  splits into 3 * 2^m cells (a, x), each holding the 2^f words with f = b - m free coordinates.
  M_{a,x} = number of codewords in the cell.
* Every word of cell (a,x) must be covered. A codeword inside the cell covers itself and its f
  neighbours inside the cell; a codeword in a neighbouring cell (other ternary value, or one
  fixed binary coordinate flipped) covers exactly one word of the cell. Hence
        (1 + f) M_{a,x} + sum_{neighbours y} M_y >= 2^f,      0 <= M_{a,x} <= 2^f.
* Level m+1 splits each cell in two: M_{a,x} = M_{a,2x} + M_{a,2x+1}.
* Distributions equivalent under S_3 (relabel ternary values) x (permute and flip the fixed
  binary coordinates) are kept once (canonical form = lexicographically least image).
* At level `leaf`, each surviving distribution is completed exactly by CP-SAT (integer
  arithmetic): 3 * 2^b binary variables, the per-cell counts, and every covering constraint.

All arithmetic is integer. If no distribution survives, or every leaf is infeasible, no code
of size M exists.
"""
import argparse
import itertools
import json
import os
import sys
import time

import numpy as np
from ortools.sat.python import cp_model

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def group_perms(m):
    """All cell permutations at level m: (pi in S3) x (sigma in S_m) x (flip f in F_2^m).
    Returns an int array G of shape (|group|, 3*2^m) with image cell index of each cell."""
    n = 1 << m
    xs = np.arange(n)
    bits = (xs[:, None] >> (m - 1 - np.arange(m))[None, :]) & 1 if m else np.zeros((n, 0), int)
    perms = []
    for sigma in itertools.permutations(range(m)):
        moved = bits[:, list(sigma)] if m else bits
        base = (moved * (1 << (m - 1 - np.arange(m)))[None, :]).sum(1) if m else np.zeros(n, int)
        for f in range(n):
            xmap = base ^ f
            for pi in itertools.permutations(range(3)):
                img = np.empty(3 * n, int)
                for a in range(3):
                    img[a * n + xs] = pi[a] * n + xmap
                perms.append(img)
    return np.array(perms)


def canonical(vec, G):
    """Lexicographically least image of a distribution under the group (G: image indices)."""
    imgs = np.empty((G.shape[0], len(vec)), dtype=vec.dtype)
    imgs[np.arange(G.shape[0])[:, None], G] = vec[None, :]   # image[g(c)] = vec[c]
    order = np.lexsort(imgs.T[::-1])
    return tuple(imgs[order[0]])


def neighbours(m):
    n = 1 << m
    nb = []
    for a in range(3):
        for x in range(n):
            ys = [b * n + x for b in range(3) if b != a]
            ys += [a * n + (x ^ (1 << (m - 1 - i))) for i in range(m)]
            nb.append(ys)
    return nb


def children(parent, m, b, M):
    """All level-(m+1) distributions refining `parent` (a level-m distribution), via CP-SAT."""
    n1 = 1 << (m + 1)
    f = b - (m + 1)
    nb = neighbours(m + 1)
    mdl = cp_model.CpModel()
    v = [mdl.NewIntVar(0, 1 << f, f"v{i}") for i in range(3 * n1)]
    for a in range(3):
        for x in range(1 << m):
            mdl.Add(v[a * n1 + 2 * x] + v[a * n1 + 2 * x + 1] == int(parent[a * (1 << m) + x]))
    for c in range(3 * n1):
        mdl.Add((1 + f) * v[c] + sum(v[y] for y in nb[c]) >= (1 << f))
    out = []

    class Collect(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            out.append(np.array([self.Value(x) for x in v], dtype=np.int16))

    s = cp_model.CpSolver()
    s.parameters.enumerate_all_solutions = True
    s.parameters.num_workers = 1
    s.Solve(mdl, Collect())
    return out


def level0(b, M):
    f = b
    out = []
    for m0 in range(M + 1):
        for m1 in range(m0, M + 1):
            m2 = M - m0 - m1
            if m2 < m1:
                continue
            t = (m0, m1, m2)
            if all((1 + f) * t[a] + sum(t) - t[a] >= (1 << f) and t[a] <= (1 << f) for a in range(3)):
                out.append(np.array(t, dtype=np.int16))
    return out


def complete(dist, leaf, b, M, tlim):
    """Exact completion: is there a code with these per-cell counts? CP-SAT, integer only."""
    words = list(itertools.product(range(3), *[range(2)] * b))
    idx = {w: i for i, w in enumerate(words)}
    mdl = cp_model.CpModel()
    x = [mdl.NewBoolVar(f"x{i}") for i in range(len(words))]
    for w in words:
        nbr = [w] + [w[:i] + (v,) + w[i + 1:] for i in range(b + 1)
                     for v in range(3 if i == 0 else 2) if v != w[i]]
        mdl.AddBoolOr([x[idx[u]] for u in nbr])
    n = 1 << leaf
    cells = {}
    for w in words:
        xb = 0
        for bit in w[1:1 + leaf]:
            xb = 2 * xb + bit
        cells.setdefault(w[0] * n + xb, []).append(x[idx[w]])
    for c, vs in cells.items():
        mdl.Add(sum(vs) == int(dist[c]))
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = tlim
    s.parameters.num_workers = 1
    st = s.Solve(mdl)
    if st == cp_model.INFEASIBLE:
        return "infeasible", None
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return "feasible", [w for w in words if s.Value(x[idx[w]])]
    return "timeout", None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("b", type=int, help="binary coordinates (one ternary coordinate)")
    ap.add_argument("M", type=int, help="code size to test")
    ap.add_argument("--leaf", type=int, default=3, help="level at which CP-SAT completes")
    ap.add_argument("--count-only", action="store_true", help="only count distributions per level")
    ap.add_argument("--tlim", type=float, default=600.0)
    a = ap.parse_args()

    t0 = time.time()
    dists = level0(a.b, a.M)
    print(f"K(1,{a.b},1), M={a.M}: level 0: {len(dists)} distributions", flush=True)
    for m in range(a.leaf):
        G = group_perms(m + 1)
        seen = {}
        for d in dists:
            for c in children(d, m, a.b, a.M):
                key = canonical(c, G)
                seen.setdefault(key, np.array(key, dtype=np.int16))
        dists = list(seen.values())
        print(f"  level {m + 1}: {len(dists)} inequivalent distributions ({time.time() - t0:.1f}s)", flush=True)
        if not dists:
            break
    if a.count_only or not dists:
        verdict = "no code" if not dists else "counted"
        print(f"done: {verdict} ({time.time() - t0:.1f}s)")
        return 0
    found = timeouts = 0
    for i, d in enumerate(dists):
        st, code = complete(d, a.leaf, a.b, a.M, a.tlim)
        if st == "feasible":
            found += 1
            path = os.path.join(ROOT, "results", "codes", f"K_1_{a.b}_1_M{a.M}_exact.txt")
            with open(path, "w") as fh:
                fh.write("\n".join("".join(map(str, w)) for w in code) + "\n")
            print(f"  leaf {i}: CODE FOUND, saved {path}", flush=True)
            break
        timeouts += st == "timeout"
    verdict = "CODE EXISTS" if found else ("no code" if not timeouts else f"undecided ({timeouts} timeouts)")
    print(f"done: {verdict}; {len(dists)} leaves at level {a.leaf} ({time.time() - t0:.1f}s)")
    with open(os.path.join(ROOT, "results", "exact_dist.jsonl"), "a") as fh:
        fh.write(json.dumps({"b": a.b, "M": a.M, "leaf": a.leaf, "leaves": len(dists),
                             "verdict": verdict, "seconds": round(time.time() - t0, 1)}) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
