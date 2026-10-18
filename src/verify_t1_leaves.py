"""Definition-level check of exact_dist's leaf list for K(1,B,1) at level 2.

1. For level L = 0, 1, 2: cells are (ternary value, first L binary coordinates). For each cell c
   and word u, count the words of c whose radius-1 ball contains u (Hamming distance only);
   check it is constant on each cell, giving coefficient matrices A_L with A_L M >= |cell|.
2. Enumerate by brute force every level-2 distribution (12 cells, each <= 2^(B-2), sum M) whose
   level-2, level-1 and level-0 aggregates satisfy A_2, A_1, A_0 respectively.
3. Check every such distribution is equivalent, under the 48 maps (permute the 3 ternary values)
   x (permute/flip the 2 fixed binary coordinates), to one of exact_dist's leaves, and vice versa.
Usage: python src/verify_t1_leaves.py B M
"""
import itertools, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exact_dist as E
B, M = int(sys.argv[1]), int(sys.argv[2])
words = list(itertools.product(range(3), *[range(2)] * B))
dist = lambda u, v: sum(p != q for p, q in zip(u, v))
def cellof(w, L):
    x = 0
    for bit in w[1:1 + L]: x = 2 * x + bit
    return w[0] * (1 << L) + x
A = {}
for L in (0, 1, 2):
    nc = 3 << L
    mat = np.zeros((nc, nc), dtype=int)
    for c in range(nc):
        per = {}
        for u in words:
            per[u] = sum(1 for w in words if cellof(w, L) == c and dist(w, u) <= 1)
        for d in range(nc):
            vals = {per[u] for u in words if cellof(u, L) == d}
            assert len(vals) == 1, (L, c, d, vals)
            mat[c, d] = vals.pop()
    A[L] = mat
    print(f"level {L} coefficients (row 0): {mat[0].tolist()}")
size = {L: 1 << (B - L) for L in (0, 1, 2)}
# level-1 distributions by brute force: 6 cells; level-1 cell (a, b1) = index a*2 + b1
def compositions(n, k):                     # all k-tuples of naturals summing to n
    if k == 1:
        yield (n,); return
    for first in range(n + 1):
        for rest in compositions(n - first, k - 1): yield (first,) + rest
lvl1 = []
for d in compositions(M, 6):
    if max(d) > size[1]: continue
    v = np.array(d)
    l0 = v.reshape(3, 2).sum(1)
    if np.all(A[1] @ v >= size[1]) and np.all(A[0] @ l0 >= size[0]): lvl1.append(v)
raw = set()
for v in lvl1:                               # split each level-1 cell into its two level-2 cells
    for split in itertools.product(*[range(int(t) + 1) for t in v]):
        w = []
        for t, s in zip(v, split): w += [s, int(t) - s]
        w = np.array(w)
        if w.max() <= size[2] and np.all(A[2] @ w >= size[2]): raw.add(tuple(int(t) for t in w))
print(f"brute force: {len(lvl1)} level-1 distributions, {len(raw)} level-2 distributions")
# leaves from exact_dist
dists = E.level0(B, M)
for m in range(2):
    G = E.group_perms(m + 1); seen = {}
    for d in dists:
        for c in E.children(d, m, B, M): seen.setdefault(E.canonical(c, G), c)
    dists = list(seen.values())
leaves = {tuple(int(t) for t in d) for d in dists}
# the 48 symmetries written out directly: cell index a*4 + x, x = 2*b1 + b2
def images(d):
    for pi in itertools.permutations(range(3)):
        for swap in (False, True):
            for f1 in (0, 1):
                for f2 in (0, 1):
                    out = [0] * 12
                    for a in range(3):
                        for b1 in (0, 1):
                            for b2 in (0, 1):
                                c1, c2 = (b2, b1) if swap else (b1, b2)
                                out[pi[a] * 4 + 2 * (c1 ^ f1) + (c2 ^ f2)] = d[a * 4 + 2 * b1 + b2]
                    yield tuple(out)
unmatched = [d for d in raw if not any(im in leaves for im in images(d))]
hit = {im for d in raw for im in images(d) if im in leaves}
print(f"leaves from exact_dist: {len(leaves)}; raw unmatched: {len(unmatched)}; leaves never hit: {len(leaves - hit)}")
