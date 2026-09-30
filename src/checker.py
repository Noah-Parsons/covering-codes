"""Independent brute-force coverage checker for mixed ternary/binary covering codes.

Deliberately simple and slow. It shares NO code with the search (src/search.py):
every word of the space is enumerated with itertools.product, and for each word
the Hamming distance to every codeword is computed directly until one within
radius R is found.

Word format: a string of t ternary digits (0/1/2) followed by b binary digits
(0/1). Binary K(n,R) is t=0, b=n. Ternary K_3(n,R) is t=n, b=0.

Code file format: one codeword per line; blank lines and lines starting with
'#' are ignored. Separators ' ', ',', '|' inside a codeword are ignored.

Usage:
    python src/checker.py FILE t b R
Exit code 0 iff the file is a valid covering code with covering radius <= R.
"""
import itertools
import sys


def parse_code(lines, t, b):
    code = []
    for raw in lines:
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        s = s.replace(" ", "").replace(",", "").replace("|", "")
        if len(s) != t + b:
            raise ValueError(f"codeword {s!r} has length {len(s)}, expected {t + b}")
        for i, ch in enumerate(s):
            allowed = "012" if i < t else "01"
            if ch not in allowed:
                raise ValueError(f"codeword {s!r}: symbol {ch!r} at position {i} not in {allowed}")
        code.append(tuple(int(ch) for ch in s))
    return code


def hamming(u, v):
    return sum(1 for a, c in zip(u, v) if a != c)


def uncovered_words(code, t, b, R):
    """Return the list of words at distance > R from every codeword."""
    alphabets = [range(3)] * t + [range(2)] * b
    bad = []
    for w in itertools.product(*alphabets):
        if not any(hamming(w, c) <= R for c in code):
            bad.append(w)
    return bad


def check(code, t, b, R):
    """Return (ok, n_distinct_codewords, n_uncovered)."""
    distinct = set(code)
    bad = uncovered_words(list(distinct), t, b, R)
    return (len(bad) == 0, len(distinct), len(bad))


def main(argv):
    if len(argv) != 5:
        print(__doc__)
        return 2
    path, t, b, R = argv[1], int(argv[2]), int(argv[3]), int(argv[4])
    with open(path) as f:
        code = parse_code(f, t, b)
    ok, m, nbad = check(code, t, b, R)
    dup = len(code) - m
    print(f"file={path} t={t} b={b} R={R} codewords={m} (duplicates dropped: {dup}) "
          f"space={3**t * 2**b} uncovered={nbad} -> {'VALID' if ok else 'INVALID'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
