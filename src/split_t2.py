"""Is there a covering code of radius 1 with exactly M words in F_3^2 x F_2^B?

(A smaller code can be padded to M words, so "no code of size M" proves K(2,B,1) >= M+1.)

Split the space into the 3x3 grid of blocks given by the two ternary coordinates; each block
is a copy of Q_B (2^B words). A codeword in block c covers itself and its B neighbours inside
c; a codeword in a block sharing a row or column with c covers exactly one word of c; blocks
differing in both ternary coordinates cover nothing in c. So the block counts M_c satisfy

    (1 + B) M_c + sum of M_d over the 4 blocks d sharing a row or column with c  >=  2^B,

with 0 <= M_c <= 2^B and sum M_c = M. Distributions equivalent under permuting the symbols of
either ternary coordinate and swapping the two coordinates (72 maps) are solved once. Each
distribution becomes an exact covering problem with fixed block counts, solved with SCIP.

Usage:  python src/split_t2.py B M [--tlim 120]
Writes a found code to results/codes/ and one summary line to results/split_t2.jsonl.
"""
import argparse
import itertools
import json
import os
import sys
import time

from pyscipopt import Model, quicksum

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CELLS = [(a, b) for a in range(3) for b in range(3)]
NB = {c: [d for d in CELLS if d != c and (d[0] == c[0] or d[1] == c[1])] for c in CELLS}


def symmetries():
    out = []
    for pr in itertools.permutations(range(3)):
        for pc in itertools.permutations(range(3)):
            for tr in (False, True):
                out.append({c: ((pc[c[1]], pr[c[0]]) if tr else (pr[c[0]], pc[c[1]])) for c in CELLS})
    return out


def canonical(d, G):
    """Least image of a {cell: count} distribution under the 72 symmetries."""
    return min(tuple(d[next(k for k in CELLS if g[k] == c)] for c in CELLS) for g in G)


def distributions(B, M):
    """Fast version. The block inequality (1+B) M_c + sum_{N(c)} M_d >= 2^B is the same as
        (B-1) M_c + R_a + C_b >= 2^B,
    with R_a, C_b the totals of c's row and column (the row and column each contain M_c once).
    Row and column permutations can sort both margins, so only sorted (R, C) are tried; each
    3x3 matrix with those margins is then built from its four top-left entries. Results are
    reduced by all 72 symmetries, exactly as distributions_slow does."""
    Q = 1 << B
    G = symmetries()
    found = set()

    def sorted_triples(total, cap):
        return [(x, y, total - x - y) for x in range(total + 1) for y in range(x, total + 1)
                if total - x - y >= y and total - x - y <= cap]

    for R in sorted_triples(M, 3 * Q):
        for C in sorted_triples(M, 3 * Q):
            lo = [[max(0, -(-(Q - R[a] - C[b]) // (B - 1))) for b in range(3)] for a in range(3)]
            hi = [[min(Q, R[a], C[b]) for b in range(3)] for a in range(3)]
            if sum(map(sum, lo)) > M:
                continue
            for m00 in range(lo[0][0], hi[0][0] + 1):
                for m01 in range(lo[0][1], hi[0][1] + 1):
                    m02 = R[0] - m00 - m01
                    if not lo[0][2] <= m02 <= hi[0][2]:
                        continue
                    for m10 in range(lo[1][0], hi[1][0] + 1):
                        m20 = C[0] - m00 - m10
                        if not lo[2][0] <= m20 <= hi[2][0]:
                            continue
                        for m11 in range(lo[1][1], hi[1][1] + 1):
                            m12 = R[1] - m10 - m11
                            m21 = C[1] - m01 - m11
                            m22 = R[2] - m20 - m21
                            if not (lo[1][2] <= m12 <= hi[1][2] and lo[2][1] <= m21 <= hi[2][1]
                                    and lo[2][2] <= m22 <= hi[2][2]):
                                continue
                            d = dict(zip(CELLS, (m00, m01, m02, m10, m11, m12, m20, m21, m22)))
                            found.add(canonical(d, G))
    return sorted(found)


def distributions_slow(B, M):
    """Original version: every split of M into nine cells, then the block inequality."""
    Q = 1 << B
    G = symmetries()
    found = set()

    def rec(i, d, left):
        if i == 8:
            d = dict(d)
            d[CELLS[8]] = left
            if left <= Q and all((1 + B) * d[c] + sum(d[e] for e in NB[c]) >= Q for c in CELLS):
                found.add(min(tuple(d[next(k for k in CELLS if g[k] == c)] for c in CELLS) for g in G))
            return
        for v in range(min(left, Q) + 1):
            d[CELLS[i]] = v
            rec(i + 1, d, left - v)

    rec(0, {}, M)
    return sorted(found)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("B", type=int)
    ap.add_argument("M", type=int)
    ap.add_argument("--tlim", type=float, default=120.0)
    a = ap.parse_args()

    words = [(x, y) + r for x in range(3) for y in range(3) for r in itertools.product(range(2), repeat=a.B)]
    idx = {w: i for i, w in enumerate(words)}

    def ball(w):
        return [w] + [w[:i] + (v,) + w[i + 1:] for i in range(2 + a.B)
                      for v in range(3 if i < 2 else 2) if v != w[i]]

    dists = distributions(a.B, a.M)
    print(f"K(2,{a.B},1), M={a.M}: {len(dists)} inequivalent block distributions", flush=True)
    t0 = time.time()
    tally, code_file = {}, None
    for k, dist in enumerate(dists):
        m = Model()
        m.hideOutput()
        x = [m.addVar(vtype="B") for _ in words]
        for w in words:
            m.addCons(quicksum(x[idx[u]] for u in ball(w)) >= 1)
        for ci, c in enumerate(CELLS):
            m.addCons(quicksum(x[idx[w]] for w in words if w[:2] == c) == dist[ci])
        m.setParam("limits/time", a.tlim)
        m.optimize()
        st = m.getStatus()
        tally[st] = tally.get(st, 0) + 1
        print(f"  {k}: {dist} -> {st} ({m.getSolvingTime():.1f}s)", flush=True)
        if m.getNSols() > 0:
            code = [w for w in words if m.getVal(x[idx[w]]) > 0.5]
            code_file = os.path.join(ROOT, "results", "codes", f"K_2_{a.B}_1_M{len(code)}_split.txt")
            with open(code_file, "w") as fh:
                fh.write("\n".join("".join(map(str, w)) for w in code) + "\n")
            print(f"  CODE FOUND ({len(code)} words): {code_file}", flush=True)
            break
    verdict = ("code found" if code_file else
               "no code" if set(tally) == {"infeasible"} else "undecided")
    print(f"done: {verdict} {tally} ({time.time() - t0:.0f}s)", flush=True)
    with open(os.path.join(ROOT, "results", "split_t2.jsonl"), "a") as fh:
        fh.write(json.dumps({"B": a.B, "M": a.M, "distributions": len(dists), "verdict": verdict,
                             "tally": tally, "tlim": a.tlim, "seconds": round(time.time() - t0, 1)}) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
