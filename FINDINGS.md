# Findings

No published upper bound has been beaten. Every code this project has produced that passes
`src/checker.py` has size equal to or larger than the table value.

## Round 1–2 (2026-09-29): tabu search on raw codewords (`src/search.py`)

- Reproduced K(4,4,2) <= 24, K_3(7,2) <= 34, K_3(6,1) <= 73, K_3(8,3) <= 27. Failed to
  reproduce K_3(7,1) <= 186 (reached 192), so the search is weaker than the published ones.
- At record − 1 the best results were 2 uncovered words for K(4,3,1) at M=71 and K(5,1,1) at
  M=53, and 5 uncovered for K(4,5,2) at M=47 (two 10-minute runs gave 6 and 5).
- Full run list: `results/runs.csv`.

## Round 3 (2026-10-13): group-invariant codes by exact ILP (`src/orbit_ilp.py`)

Method: pick a random group G of isometries (coordinate permutations within each block,
x -> 2x on ternary coordinates, translations), restrict to codes that are unions of
G-orbits, and solve the resulting weighted set cover exactly with HiGHS, with the extra
constraint "size <= record − 1". "Infeasible" is then a proof that no G-invariant code beats
the record. Groups left undecided at 12 s were rerun as a plain minimisation for 20–45 s.

Validation: a random group of order 24 gave an 81-word code for K_3(8,2), the published
record, verified by the checker (`results/codes/K_8_0_2_M81_orbit.txt`).

| Target | Record | Groups tried | Proved no code below record | Undecided | Best size found |
|--------|--------|--------------|-----------------------------|-----------|-----------------|
| K_3(8,2) | 81  | 184 | 179 | 5 | 81 (optimal for 2 groups) |
| K(8,1,2) | 162 | 67  | 65  | 2 | 162 (optimal for 4 groups) |
| K(4,5,2) | 48  | 71  | 64  | 7 | 48 (optimal for 1 group) |
| K(7,2,2) | 108 | 360 | 359 | 1 | 120 |
| K(6,4,2) | 144 | 192 | 191 | 1 | 204 |

Every group, its generators, and the solver outcome is in `results/orbit_runs.jsonl`.

What this does and does not show:

- Most "infeasible" results are cheap: a group whose orbits are all large cannot give a small
  code, and the solver sees that in milliseconds. They say little about the bound.
- The informative cases are the groups where the solver had real work to do. There the
  optimum was exactly the record (81, 162, 48) in 7 groups and never below it.
- The groups were sampled at random, not enumerated up to conjugacy, so this is not a
  classification. The undecided groups were not resolved.
- Symmetric-code searches of this kind were already run in the 1990s (Östergård and
  co-authors), with tabu search rather than exact solving.

## Round 4 (2026-10-13): which symmetries can an improving code have? (`src/classify.py`)

Change of question: instead of looking for a smaller code, prove what a smaller code would
have to look like. A code with a non-trivial automorphism is invariant under some isometry of
prime order, so it is enough to take one representative of every conjugacy class of
prime-order isometries and show that no invariant code of size <= record − 1 exists. The
class representatives are checked against brute-force enumeration of the whole isometry group
on five small spaces (`tests/test_classify.py`).

| Target | Record | Classes (order 2 / odd) | Odd-order classes excluded | Order-2 classes excluded |
|--------|--------|-------------------------|----------------------------|--------------------------|
| K(1,7,1) | 48 | 39 / 7  | 7 of 7  | 39 of 39 |
| K(3,4,1) | 48 | 53 / 9  | 9 of 9  | 47 of 53 (1 timed out at 240 s, 5 not yet run) |
| K(5,1,1) | 54 | 23 / 9  | 9 of 9  | 2 of 23 (3 timed out at 300 s, 18 not yet run) |
| K(4,3,1) | 72 | 53 / 13 | 2 of 13 (2 timed out) | not started |

State when the runs were stopped on 2026-10-13 (they resume from `results/classify.jsonl`):

- K(1,7,1): every class excluded. A code with fewer than 48 words has trivial automorphism
  group. First-pass result; recheck 11 of 46 done.
- K(3,4,1) and K(5,1,1): a code below the record has no automorphism of odd order.
- The hard order-2 classes were solved by minimising with no size cap and reading HiGHS's
  proven lower bound (`solve_min`), which closes far faster than the capped feasibility form.

Recheck (`src/recheck.py`, `results/recheck.jsonl`): rebuilds each excluded class with
independent code and a random conjugate of the generator. With HiGHS, 11 of the 46 K(1,7,1)
classes are reconfirmed so far, none contradicted. OR-tools CP-SAT confirmed the 3 easiest and
timed out at 150 s on the next two, so it is not usable as the second solver.

CBC second-solver recheck (started 2026-10-15). CBC 2.9 from `pip install pulp[cbc]` (the
`cbcbox` package). PuLP 4 changed its API completely, so `recheck.py` writes an LP file and
runs `cbc.exe` directly, minimising with no size cap and comparing the proven optimum with
47. Positive control first: on two order-2 classes CBC found a 48-word invariant code (the
record) and proved 48 optimal, matching HiGHS's bound. CBC is roughly 30x slower than HiGHS
on these problems; a size cutoff of 47.5 made it slower still. Result: 14 of 46 confirmed, 5
timed out at 30 min, then stopped in favour of SCIP. Logs: `results/logs/`.

**SCIP second-solver recheck, complete (2026-10-15). All 46 K(1,7,1) classes confirmed.**
SCIP 10.0 (`pip install pyscipopt`), `recheck.py --solver scip`, minimising with no size cap
on a random conjugate of each generator (seed 2026), rebuilt by `recheck.py`'s own code. Same
positive control first: with 48 allowed it finds a 48-word invariant code; with 47 it proves
none exists. Every class: proven optimum or proven lower bound above 47; longest class 626 s.
So the K(1,7,1) statement now rests on two independent solvers built by different teams (HiGHS
first pass, SCIP recheck), plus CBC on 14 classes and HiGHS-on-a-conjugate on 11, all agreeing:

> Every covering code of radius 1 in F_3 x F_2^7 with at most 47 words has trivial
> automorphism group.

Still a floating-point branch-and-bound verdict from each solver, not an exact certificate.

Caveats: "excluded" means HiGHS reported the integer program infeasible. That is a
floating-point branch-and-bound verdict, not a certified proof; a publishable statement needs
an independent second solver or an exact check. Classes "tried" but not excluded hit the
time limit and are open, not counterexamples. Per-class log: `results/classify.jsonl`.

### Novelty check for round 4 (2026-10-13, web search, not exhaustive)

- No published statement was found that codes below the record for K(1,7,1) or K(3,4,1) must
  be asymmetric, and no classification of codes for these entries.
- The method is not new. Östergård and Weakley (Des. Codes Cryptogr. 16, 1999) search for
  covering codes with a prescribed automorphism group, but heuristically (tabu), to build
  codes rather than to exclude them. Exhaustive exclusion of prime-order automorphisms by
  integer programming is standard in neighbouring areas (self-dual codes, subspace codes).
- No update of the mixed binary/ternary table after Kéri's (2011) was found, so 42–48 and
  44–48 appear to be current. Gijswijt–Polak (arXiv:2504.01932) does not treat mixed codes.
- Not read (paywalled): Bertolo–Östergård–Weakley (J. Combin. Des. 12, 2004) and
  Östergård–Hämäläinen (1997). Either could contain a remark that overlaps.
- Related and recent: arXiv:2608.19872 (August 2026) improves 84 entries of Kéri's q-ary
  tables for alphabet sizes 5–21 by local search and semidefinite bounds.

## Assessment

Three different attacks (plain annealing, tabu from the record, exact search under symmetry)
have all stopped at the table values. The evidence so far is that these bounds are hard to
improve with laptop-scale search, not that they are neglected.

## Round 5 (2026-10-16): exact attempt at K(1,7,1) = 48 (`src/exact_dist.py`)

Prompted by P. Östergård's reply (2026-10-16): symmetry-exclusion results like round 4 are
not published unless the problem is of very wide interest; he expects K(1,7,1) = 48, provable
the way Östergård and Blass proved K(9,1) = 62 (IEEE Trans. Inform. Theory 47, 2001).

Method: their codeword-distribution refinement adapted to one ternary coordinate (fix the
ternary coordinate, then binary coordinates one at a time; cell inequality
(1+f) M_x + sum over neighbouring cells M_y >= 2^f; dedupe under S_3 x hyperoctahedral group),
then exact completion of each surviving distribution with CP-SAT. Integer arithmetic only.

Validation (all correct, codes confirmed by `src/checker.py`):
K(1,4,1): no code at 7, code at 8. K(1,5,1): none at 15, code at 16. K(1,6,1): none at 23,
code at 24. Each in under 10 s.

K(1,7,1), M = 47 (measured, not run to the end):

| level | inequivalent distributions | exact completion of a sample |
|---|---|---|
| 0 | 16 | |
| 1 | 429 | |
| 2 | 21,396 | 8 of 8 timed out at 40 s |
| 3 | ~15 million raw (sampled; a few million after dedupe) | 10 of 15 infeasible in 0.04-14 s, 5 timed out at 20 s |
| 4 | splitting 20 level-3 cases did not finish in ~8 min | not measured |

Estimate: months of CPU time on this laptop at level 3; level 4 multiplies the count faster
than it eases each case. Not feasible here without a new idea (much stronger pruning than the
cell inequality) or the full orderly-generation engineering of the 2001 paper plus far more
compute. Stopped before any long run.

## Round 6 (2026-10-16): new lower-bound ideas, and a bridge to rainbow domination

- Column observation: a codeword covers all 3 words of its column, so K(1,b,1) is the
  3-rainbow domination number of Q_b (Brešar–Henning–Rall 2008: gamma_rk(G) = gamma(G □ K_k)).
- Colour-blind relaxation (g(v) in 0..3, empty column needs neighbour total >= 3): exact for
  b = 4, 6 but 14 < 16 for b = 5, and has a weight-47 solution for b = 7 (not colourable). Too weak.
- Symmetry-reduced level-1 SDP / Lovász–Schrijver bound (16 pair classes): 39.10 for K(1,7,1),
  below the known 42. (A first version with a sign error gave invalid bounds above known values;
  caught by validating on b = 3..6.) The strong Gijswijt–Polak hierarchy, not yet done for mixed
  codes, would be a major project.
- Bridge confirmed by direct computation (`src/rainbow.py`, graph definition, SCIP):
  gamma_r2(Q_n) = 2, 2, 4, 7, 12, 16 for n = 1..6 (= K(n+1,1)) and
  gamma_r3(Q_n) = 2, 3, 6, 8, 16 for n = 1..5 (= K(1,n,1)). No published hypercube values of
  rainbow domination were found (web search, 2026-10-16), so the coding tables give values that
  literature lacks: gamma_r2(Q_n) exact for n <= 8 (gamma_r2(Q_8) = K(9,1) = 62), gamma_r3(Q_n)
  exact for n <= 6 (gamma_r3(Q_6) = 24), 42 <= gamma_r3(Q_7) <= 48.
