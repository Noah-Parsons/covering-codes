"""Reconstruct record codes from their stated constructions.

  extend   C (t,b) -> C x Z2 (t, b+1)  or  Z3 x C (t+1, b)
           realises the table's generic bound K(t,b,R) <= 2K(t,b-1,R) / 3K(t-1,b,R).
  linear   random search for a ternary [n,k] linear code with covering radius R
           (every syndrome in F_3^(n-k) is a combination of <= R columns of H).

Outputs plain codeword files; verify them with src/checker.py before use.
"""
import argparse
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from checker import parse_code  # only for reading files


def extend(inp, out, t, b, kind):
    with open(inp) as f:
        code = ["".join(map(str, w)) for w in parse_code(f, t, b)]
    if kind == "binary":
        new = [w + s for w in code for s in "01"]
        nt, nb = t, b + 1
    else:
        new = [s + w for w in code for s in "012"]
        nt, nb = t + 1, b
    with open(out, "w") as f:
        f.write(f"# extension ({kind}) of {os.path.basename(inp)}: t={nt} b={nb}, M={len(new)}\n")
        f.write("\n".join(new) + "\n")
    print(f"wrote {out}: M={len(new)} t={nt} b={nb}")


def linear(n, k, R, out, seed, tries):
    rng = np.random.default_rng(seed)
    r = n - k
    nsyn = 3 ** r
    weights = [3 ** i for i in range(r)]
    for trial in range(tries):
        H = rng.integers(0, 3, size=(r, n))
        cols = [H[:, j] for j in range(n)]
        seen = {0}
        for w in range(1, R + 1):
            for pos in itertools.combinations(range(n), w):
                for coef in itertools.product((1, 2), repeat=w):
                    s = sum(c * cols[p] for c, p in zip(coef, pos)) % 3
                    seen.add(int(np.dot(s, weights)))
        if len(seen) == nsyn:
            words = [w for w in itertools.product(range(3), repeat=n)
                     if not (H.dot(np.array(w)) % 3).any()]
            with open(out, "w") as f:
                f.write(f"# ternary linear [{n},{k}] code, covering radius <= {R}, H rows: "
                        + " ".join("".join(map(str, row)) for row in H) + "\n")
                f.write("\n".join("".join(map(str, w)) for w in words) + "\n")
            print(f"wrote {out}: M={len(words)} after {trial + 1} random H")
            return 0
    print("no linear code found")
    return 1


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extend")
    e.add_argument("inp"); e.add_argument("out")
    e.add_argument("t", type=int); e.add_argument("b", type=int)
    e.add_argument("kind", choices=["binary", "ternary"])
    li = sub.add_parser("linear")
    li.add_argument("n", type=int); li.add_argument("k", type=int); li.add_argument("R", type=int)
    li.add_argument("out"); li.add_argument("--seed", type=int, default=1)
    li.add_argument("--tries", type=int, default=2000)
    a = ap.parse_args()
    if a.cmd == "extend":
        extend(a.inp, a.out, a.t, a.b, a.kind)
    else:
        return linear(a.n, a.k, a.R, a.out, a.seed, a.tries)


if __name__ == "__main__":
    sys.exit(main())
