# Targets

Notation: K(t, b, R) is the smallest code with t ternary and b binary coordinates and
covering radius R. K_3(n, R) = K(n, 0, R). "Space" = 3^t 2^b.

Bounds are from Kéri's tables (mixed_tables.pdf, 3_tables.pdf, accessed 2026-09-29; the
index page's last listed update is 2011-11-21). Ternary lower bounds marked (GP) are the
improved values from Gijswijt and Polak, arXiv:2504.01932 (v2, June 2026), which compares
against Kéri's tables and lists the same upper bounds, so the upper bounds below are, as far
as I could find, still current.

Upper-bound key letters (mixed table): a = the generic recursion K(t,b,R) <= 2K(t,b-1,R)
or <= 3K(t-1,b,R); b = Hämäläinen–Rankinen (1991). Ternary table: u = Hämäläinen–Rankinen
(1991).

| # | Entry | Lower | Upper | Gap | Upper-bound provenance | Space | Why chosen |
|---|-------|-------|-------|-----|------------------------|-------|------------|
| 1 | K(4,5,2) | 32 | 48 | 16 | a: 2·K(4,4,2) = 2·24 (no dedicated search listed) | 2592 | Biggest relative gap (1.5x) among small mixed R=2 entries; bound is a plain doubling |
| 2 | K(1,7,1) | 42 | 48 | 6 | a: 3·K(7,1) = 3·16 (Hamming code x Z3) | 384 | Tiny space, generic upper bound |
| 3 | K(3,4,1) | 44 | 48 | 4 | a: 2·K(3,3,1) = 2·24 | 432 | Tiny space, generic upper bound |
| 4 | K(4,3,1) | 60 | 72 | 12 | a: 2·K(4,2,1) = 2·36 | 648 | Generic upper bound, gap 12 |
| 5 | K(5,1,1) | 45 | 54 | 9 | a: 2·K_3(5,1) = 2·27 | 486 | Generic upper bound |
| 6 | K(7,2,2) | 72 | 108 | 36 | a: 2·K(7,1,2) = 2·54 | 8748 | Large gap, generic upper bound |
| 7 | K(6,4,2) | 94 | 144 | 50 | a: 2·K(6,3,2) = 2·72 | 11664 | Large gap, generic upper bound |
| 8 | K(8,1,2) | 99 | 162 | 63 | a: 2·K_3(8,2) = 2·81 | 13122 | Large gap, generic upper bound |
| 9 | K_3(8,3) | 16 (GP; table 14) | 27 | 11 | u: Hämäläinen–Rankinen 1991 (a [8,3]_3 linear code works) | 6561 | Old (1991) bound, pure ternary |
| 10 | K_3(8,2) | 58 (GP; table 54) | 81 | 23 | u: Hämäläinen–Rankinen 1991 (a [8,4]_3 linear code works) | 6561 | Old bound, large gap, pure ternary |

Caveat on "less-worked": the R=1 mixed entries (2–5) have tiny spaces and were the subject of
the heavy computations behind Bertolo–Östergård–Weakley (2004) and later lower-bound work;
an 'a' upper bound surviving there most likely means nobody could beat it, not that nobody
tried. They were included because they are cheap to test.

Reproduction targets (known records, used to validate the search): K(4,4,2) <= 24
(Hämäläinen–Rankinen 1991), K_3(7,2) <= 34 (Östergård–Hämäläinen 1993), K_3(6,1) <= 73
(Laarhoven–Aarts–van Lint–Wille 1989), K_3(8,3) <= 27, and (failed) K_3(7,1) <= 186.
