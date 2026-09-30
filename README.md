# covering-codes

This is a small project that tries to beat published upper bounds on covering codes
K(t, b, R): t ternary and b binary coordinates, covering radius R. Binary K(n,R) is t=0 and
ternary K_3(n,R) is b=0. The bounds come from Gerzson Kéri's tables
(https://old.sztaki.hu/~keri/codes/index.htm).

## Layout

- `src/checker.py`: independent brute-force verifier. **Nothing counts unless this passes.**
- `src/search.py`: tabu local search (numba) with exact incremental move evaluation and a
  descending mode.
- `src/orbit_ilp.py`: codes invariant under a random symmetry group, solved exactly as a
  weighted set cover on orbits (scipy/HiGHS); logs to `results/orbit_runs.jsonl`.
- `FINDINGS.md`: what has been tried and what came of it.
- `src/constructions.py`: rebuilds record codes (C×Z2, Z3×C extensions; random ternary
  linear codes with covering radius R).
- `tests/test_checker.py`: checker unit tests on known exact values.
- `results/runs.csv`: one row per search run (target, start, M, seed, tenure, iterations,
  wall time, best uncovered count).
- `results/codes/`: coverings found by the search; `results/records/`: reconstructed
  record codes.
- `sources/`: downloaded copies of Kéri's table PDFs (accessed 2026-09-29).
- `TARGETS.md`: the target shortlist.

## Usage

Requirements: Python 3 with numpy and numba (the checker and tests need only the standard
library).

```
python -m pytest tests -q                         # checker unit tests
python src/checker.py FILE t b R                  # verify a code file (exit 0 = valid)

# descending search from a random start of size 56 down to target 47, 60 s limit
python src/search.py 4 5 2 --start 56 --target 47 --time 60 --seed 1

# take a record code, delete one codeword, try to repair coverage
python src/search.py 4 5 2 --init results/records/K_4_5_2_M48.txt --delete 1 --target 47 --time 60

# rebuild records
python src/constructions.py linear 8 4 2 out.txt            # [8,4]_3 code, R=2 (M=81)
python src/constructions.py extend in.txt out.txt 4 4 binary  # C x Z2
```

Code files contain one codeword per line: ternary digits first, then binary digits. Lines
starting with `#` are comments. `--time` is capped at 600 s per run.
