# covering-codes

Computer searches and lower-bound proofs for covering codes K(t, b, R): the minimum number of
codewords in F_3^t x F_2^b such that every word is within Hamming distance R of a codeword
(t ternary and b binary coordinates). The reference values are Gerzson Kéri's tables
(https://old.sztaki.hu/~keri/codes/index.htm, last updated 2011).

## Results

| Entry | Kéri's table | This repository | Write-up |
|---|---|---|---|
| K(2,5,1) | 32–36 | **K(2,5,1) >= 34** | `paper/k251.tex` |
| K(1,7,1) | 42–48 | **K(1,7,1) >= 45** | `paper/k171_45.tex` (earlier version, >= 44: `paper/k171.tex`) |
| K(1,7,1) | — | any code with at most 47 words has trivial automorphism group | `paper/note.tex` |
| rainbow domination | — | gamma_r3(Q_n) = K(1,n,1) and gamma_r2(Q_n) = gamma(Q_{n+1}), with a table of values | `paper/rainbow.tex` |

The two lower bounds were found with a codeword-distribution method in the style of Östergård
and Blass (IEEE Trans. Inform. Theory 47, 2001): split the space into cells, use a counting
inequality to list the possible numbers of codewords per cell, and show every remaining case
infeasible by integer programming. Each bound is supported by:

- two independent solvers (K(2,5,1): SCIP 10 and Gurobi 13, HiGHS as a third check;
  K(1,7,1): Gurobi 13 and CPLEX 22.2, SCIP as a third check up to 43 words);
- the counting inequality recomputed from Hamming distances alone;
- a separate brute-force enumeration of the cell distributions, all of which map onto the cases
  that were solved;
- a model check: a real covering code is accepted by the per-case model, the same code with one
  word removed is rejected;
- controls on known values (K(2,3,1) = 12, K(2,4,1) = 20, K(1,5,1) = 16, K(1,6,1) = 24).

The solvers work in floating-point arithmetic, so their verdicts are strong evidence, not exact
certificates. `FINDINGS.md` records every round of the project, including the attempts that
failed.

## Layout

- `src/checker.py`: independent brute-force verifier for code files. Nothing counts unless this
  passes.
- `src/split_t2.py`: distribution search for K(2,B,1) (3x3 grid of blocks); SCIP.
- `src/recheck_t2_gurobi.py`, `src/recheck_t2_highs.py`: second- and third-solver rechecks.
- `src/verify_t2_distributions.py`: definition-level check of the K(2,B,1) distribution step.
- `src/exact_dist.py`: distribution search with one ternary coordinate (levels of cells).
- `src/leaves_t1.py`: solves the leaf cases of `exact_dist.py` with Gurobi, CPLEX or SCIP.
- `src/verify_t1_leaves.py`: definition-level check of the K(1,B,1) leaf list.
- `src/classify.py`, `src/recheck.py`: the symmetry exclusion behind `paper/note.tex`
  (HiGHS first pass, SCIP/CBC rechecks).
- `src/rainbow.py`: rainbow domination numbers of small hypercubes, computed directly.
- `src/search.py`, `src/orbit_ilp.py`, `src/constructions.py`: earlier searches for new upper
  bounds (none beat the tables).
- `paper/`: the four notes (LaTeX and PDF) and the table generator `make_table.py`.
- `results/`: logs, `*.jsonl` result files, found codes (`codes/`) and record codes (`records/`).
- `sources/`: copies of Kéri's tables and bibliography (accessed 2026-09-29).
- `tests/`: unit tests for the checker and the conjugacy-class enumeration.

## Reproducing

Python 3 with numpy, scipy, numba, ortools and pyscipopt; gurobipy and cplex for the Gurobi and
CPLEX checks (their pip packages include size-limited licences that cover these models).

```
python -m pytest tests -q                              # unit tests
python src/checker.py FILE t b R                       # verify a code file (exit 0 = valid)

# K(2,5,1) >= 34
python src/split_t2.py 5 32                            # SCIP: no 32-word code
python src/split_t2.py 5 33                            # SCIP: no 33-word code
python src/recheck_t2_gurobi.py 5 33                   # Gurobi recheck
python src/verify_t2_distributions.py 5 33             # distribution step from the definition

# K(1,7,1) >= 45  (repeat for 42 and 43 words)
python src/leaves_t1.py 7 44 --leaf 2 --solver gurobi
python src/leaves_t1.py 7 44 --leaf 2 --solver cplex
python src/verify_t1_leaves.py 7 44
```

Code files contain one codeword per line, ternary digits first, then binary digits; lines
starting with `#` are comments.
