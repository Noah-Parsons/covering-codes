"""Definition-level check of the block-distribution step behind K(2,B,1) >= M+1.

1. Coefficients from the definition: for each block c and each word u, count the words w in
   block c whose radius-1 ball contains u (pure Hamming distance). Summing the covering
   constraints over the words of c gives  sum_u coef[c][u] x_u >= |c| , i.e. a valid inequality
   on block counts IF coef depends only on u's block. Report the coefficients found.
2. Enumerate every way to put M words into the 9 blocks (each <= 2^B) with numpy, apply the
   inequality using the coefficients from step 1 (not a hand formula), and compare the resulting
   set with split_t2's raw enumeration (re-derived here from its reduced list + all 72 symmetries).
"""
import itertools, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
B, M = int(sys.argv[1]), int(sys.argv[2])
Q = 1 << B
words = list(itertools.product(range(3), range(3), *[range(2)] * B))
blocks = sorted({w[:2] for w in words})
dist = lambda u, v: sum(p != q for p, q in zip(u, v))
# --- step 1 -------------------------------------------------------------------------
coef = np.zeros((9, 9), dtype=int)          # coef[c][d]: value for u in block d (must be constant)
for ci, c in enumerate(blocks):
    per_u = {}
    for u in words:
        per_u[u] = sum(1 for w in words if w[:2] == c and dist(w, u) <= 1)
    for di, d in enumerate(blocks):
        vals = {per_u[u] for u in words if u[:2] == d}
        assert len(vals) == 1, ("coefficient not constant on a block", c, d, vals)
        coef[ci, di] = vals.pop()
print("coefficient matrix (rows: block c, cols: block of u):")
print(coef)
# --- step 2 -------------------------------------------------------------------------
# all compositions of M into 9 parts, vectorised: choose 8 bar positions among M+8 slots
from itertools import combinations
count = 0; mine = set()
bars_iter = combinations(range(M + 8), 8)
CH = 200000
while True:
    chunk = list(itertools.islice(bars_iter, CH))
    if not chunk: break
    b = np.array(chunk)
    edges = np.concatenate([np.full((len(b), 1), -1), b, np.full((len(b), 1), M + 8)], axis=1)
    parts = np.diff(edges, axis=1) - 1                     # shape (n, 9), sums to M
    count += len(parts)
    ok = (parts.max(axis=1) <= Q) & np.all(parts @ coef.T >= Q, axis=1)
    for row in parts[ok]: mine.add(tuple(int(v) for v in row))
print(f"compositions scanned: {count}; passing the definition-derived inequality: {len(mine)}")
# split_t2's raw set = all images of its reduced list under the 72 symmetries
from split_t2 import distributions, symmetries, CELLS
assert CELLS == blocks
theirs = set()
for d in distributions(B, M):
    for g in symmetries():
        img = {g[c]: d[i] for i, c in enumerate(CELLS)}
        theirs.add(tuple(img[c] for c in CELLS))
print(f"split_t2 raw set: {len(theirs)};  identical sets: {mine == theirs};  "
      f"only here: {len(mine - theirs)};  only in split_t2: {len(theirs - mine)}")
