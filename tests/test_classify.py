"""The prime-order class representatives in src/classify.py are checked against a brute-force
enumeration of the whole isometry group on small spaces."""
import itertools
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from classify import prime_classes  # noqa: E402
from orbit_ilp import gen_perm  # noqa: E402
from search import all_digits  # noqa: E402


def full_group(t, b):
    """Every isometry of Z3^t x Z2^b as a permutation of word indices."""
    digs = all_digits(t, b)
    out = []
    for pt in itertools.permutations(range(t)):
        for pb in itertools.permutations(range(t, t + b)):
            for mult in itertools.product((1, 2), repeat=t):
                for st in itertools.product(range(3), repeat=t):
                    for sb in itertools.product(range(2), repeat=b):
                        out.append(gen_perm(digs, t, b, list(pt) + list(pb),
                                            list(mult) + [1] * b, list(st) + list(sb)))
    return out


def order(p):
    q, k = p.copy(), 1
    ident = np.arange(len(p))
    while not np.array_equal(q, ident):
        q, k = p[q], k + 1
    return k


@pytest.mark.parametrize("t,b", [(2, 2), (1, 3), (3, 0), (0, 4), (2, 1)])
def test_prime_classes_match_brute_force(t, b):
    G = full_group(t, b)
    assert len({g.tobytes() for g in G}) == len(G)
    inv = [np.argsort(g) for g in G]
    prime = {g.tobytes(): g for g in G if order(g) in (2, 3, 5, 7)}
    # conjugacy classes of the prime-order elements
    cls = {}
    for key, g in prime.items():
        if key in cls:
            continue
        for h, hi in zip(G, inv):
            cls[h[g[hi]].tobytes()] = key
    n_classes = len(set(cls.values()))

    digs = all_digits(t, b)
    reps = [gen_perm(digs, t, b, *spec) for _, _, spec in prime_classes(t, b)]
    assert all(r.tobytes() in prime for r in reps)
    assert len({cls[r.tobytes()] for r in reps}) == len(reps) == n_classes
