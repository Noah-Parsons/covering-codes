"""Unit tests for the independent checker on known exact values.

Run:  python -m pytest tests -q     (or: python tests/test_checker.py)
"""
import itertools
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from checker import check, parse_code, hamming  # noqa: E402


def words(q, n):
    return list(itertools.product(range(q), repeat=n))


def test_hamming_distance():
    assert hamming((0, 1, 2), (0, 1, 2)) == 0
    assert hamming((0, 1, 2), (1, 1, 0)) == 2


def test_K2_3_1_is_2():
    code = parse_code(["000", "111"], 0, 3)
    assert check(code, 0, 3, 1)[0]
    # no single word covers F_2^3 with radius 1
    assert not any(check([w], 0, 3, 1)[0] for w in words(2, 3))


K5_1_CODE = ["00000", "00001", "00010", "01111", "10111", "11011", "11100"]


def test_K2_5_1_is_7():
    code = parse_code(K5_1_CODE, 0, 5)
    assert check(code, 0, 5, 1) == (True, 7, 0)
    # dropping any codeword breaks coverage
    for i in range(7):
        assert not check(code[:i] + code[i + 1:], 0, 5, 1)[0]
    # lower bound: exhaustively, no 6-subset covers (independent bitmask method)
    ws = words(2, 5)
    ball = [sum(1 << j for j, v in enumerate(ws) if hamming(w, v) <= 1) for w in ws]
    full = (1 << 32) - 1
    for S in itertools.combinations(range(32), 6):
        m = 0
        for i in S:
            m |= ball[i]
        assert m != full


def ternary_hamming_4():
    # tetracode [4,2,3]_3, generator rows (1,0,1,1), (0,1,1,2)
    g = [(1, 0, 1, 1), (0, 1, 1, 2)]
    return [tuple((a * g[0][k] + c * g[1][k]) % 3 for k in range(4))
            for a in range(3) for c in range(3)]


def test_K3_4_1_is_9():
    code = ternary_hamming_4()
    assert check(code, 4, 0, 1) == (True, 9, 0)
    # sphere-covering bound: 81 / 9 = 9, so 8 words can never suffice
    assert not check(code[:8], 4, 0, 1)[0]


def binary_hamming_7():
    # parity-check columns are 1..7 in binary; codeword iff XOR of positions of ones is 0
    out = []
    for w in words(2, 7):
        s = 0
        for i, bit in enumerate(w):
            if bit:
                s ^= i + 1
        if s == 0:
            out.append(w)
    return out


def test_K2_7_1_is_16():
    code = binary_hamming_7()
    assert len(code) == 16
    assert check(code, 0, 7, 1) == (True, 16, 0)
    assert not check(code[1:], 0, 7, 1)[0]


def test_mixed_and_parse_errors():
    # K(1,1,1)=2 (table): {00, 11} covers Z3 x Z2? 20 is at distance 1 from 00 -> yes
    code = parse_code(["00", "11"], 1, 1)
    ok, _, _ = check(code, 1, 1, 1)
    assert ok
    for bad in (["2"], ["02"], ["3 0"]):
        try:
            parse_code(bad, 1, 1)
        except ValueError:
            pass
        else:
            raise AssertionError(f"parse accepted {bad}")


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
