# Covering array records (strength 4)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23741408.svg)](https://doi.org/10.5281/zenodo.23741408)

Explicit covering arrays that improve the best known sizes for ten
parameter sets in the
live covering array table
[chrispbsmolen/ca-tools](https://github.com/chrispbsmolen/ca-tools)
(table checked on 2026-10-28; version 1.1 adds three arrays from the
same day; version 1.2 adds CA(825; 4, 7, 5) and CA(1050; 4, 8, 5), table
re-checked on 2026-10-29).

| Parameters | Previous best | Previous construction | New N | File |
|---|---|---|---|---|
| CA(N; 4, 100, 3) | 645 | Density Mult Cyclic (Colbourn) | **633** | `arrays/CA_633_4_100_3.txt` |
| CA(N; 4, 105, 3) | 657 | Density Mult Cyclic (Colbourn) | **644** | `arrays/CA_644_4_105_3.txt` |
| CA(N; 4, 7, 5) | 910 | cross-sum (CKRS) | **825** | `arrays/CA_825_4_7_5.txt` |
| CA(N; 4, 8, 5) | 1125 | derived from strength 5 (table entry k = 11) | **1050** | `arrays/CA_1050_4_8_5.txt` |
| CA(N; 4, 7, 6) | 1844 | orthogonal array, fuse, postop (NCK) | **1838** | `arrays/CA_1838_4_7_6.txt` |
| CA(N; 4, 8, 6) | 1975 | orthogonal array, fuse, postop (NCK) | **1966** | `arrays/CA_1966_4_8_6.txt` |
| CA(N; 4, 12, 6) | 3944 | Density (Bryce-Colbourn), postop (NCK) | **3942** | `arrays/CA_3942_4_12_6.txt` |
| CA(N; 4, 26, 6) | 6456 | Density Cyclic (Colbourn) | **6453** | `arrays/CA_6453_4_26_6.txt` |
| CA(N; 4, 27, 6) | 6606 | Density Cyclic (Colbourn) | **6581** | `arrays/CA_6581_4_27_6.txt` |

`arrays/CA_636_4_100_3.txt` (version 1.0) is kept; CA(633; 4, 100, 3)
was obtained from it by a second search. `arrays/CA_908_4_7_5.txt`
(version 1.0) is kept as well; it is superseded by CA(825; 4, 7, 5).

By deleting columns, CA(633; 4, 100, 3) also improves k = 95 to 99
(previous best 637 or 639) and CA(644; 4, 105, 3) improves k = 101 to
104 (previous best 657).

## Format

Plain text, one run per line, symbols `0` to `v-1` separated by
spaces. File names follow the table convention `CA_N_t_k_v.txt`.

## Verification

Each array is a certificate of its own bound. Check it exhaustively
with `ca_verify.c` from ca-tools:

    cc -O3 -pthread -o ca_verify ca_verify.c
    ./verify.sh ./ca_verify

Every array here was checked by `ca_verify` and by a second,
independently written checker.

## Method

### Coset unions (version 1.2: CA(825; 4, 7, 5), CA(1050; 4, 8, 5))

Take a linear code C of length k and dimension d over GF(5) whose
generator matrix has rank d on every set of 4 columns, choose m shift
vectors s_1, ..., s_m in GF(5)^k, and take the union of the cosets
s_i + C as the rows (N = m * 5^d), then remove redundant rows greedily.
For a 4-set S of columns, the projection of s_i + C onto S is a coset of
a d-dimensional subspace of GF(5)^4; the union is a CA(N; 4, k, 5)
exactly when, for every S, the shifts hit all 5^(4-d) cosets of that
subspace. The search chooses C and the shifts by local search on the
number of uncovered cosets.

CA(825; 4, 7, 5) uses d = 2: the search found 34 cosets of a [7, 2] code
(850 rows) and greedy pruning removed one whole coset, leaving 33.
CA(1050; 4, 8, 5) uses d = 2 and 42 cosets of an [8, 2] code.
Unions of cosets of linear codes belong to the known family of
group-based and hash-family constructions; the new part is the
parameters reached.

`generators/CA_825_4_7_5.c` was written by Claude Opus inside an
automated search loop and produced the 825-row array in one run of one
hour on one CPU core. `generators/CA_1050_4_8_5.c` is the same program
with k and the search budgets made parameters; it produced the
1050-row array with no further LLM involvement:

    gcc -O2 -DK=8 -DNT=70 -o gen8 generators/CA_1050_4_8_5.c
    M3=8 M2=44 P1=300 P2=900 P3=0 ./gen8 > out.txt

Checked on 2026-10-29: Colbourn's table (2019 and 2024-11 versions)
gives 910 and 1125, the NIST table gives 1280 and 1460, and no paper
found on arXiv gives a smaller array for these parameters.

### Post-optimisation (versions 1.0 and 1.1)

Post-optimisation of the published record array. Each search starts
from the explicit record array stored in ca-tools
(`table/arrays/CA_<previous N>_4_k_v.txt`), deletes a run whose
uniquely covered 4-tuples are few, and repairs the resulting holes by
local search over non-essential entries (simulated-annealing style
reassignment), repeating while time allows. The C search programs in
`generators/` were written by Claude Opus (Anthropic) inside an
automated search loop; each array was produced by the program with
the same name in one run of at most one hour on one CPU core.

To rerun a generator:

    gcc -O3 -march=native -std=gnu11 -o gen generators/CA_636_4_100_3.c -lm
    BEST_CANDIDATE=path/to/ca-tools/table/arrays/CA_645_4_100_3.txt ./gen > out.txt

The search for CA(633; 4, 100, 3) started from `arrays/CA_636_4_100_3.txt`
instead of the table array.

The searches are randomised and time-limited, so a rerun gives a valid
array of similar size, not necessarily the identical file. The arrays
in `arrays/` are the exact files that were verified.

## Citation

T. N. Sitjar (tns-research), *Covering array records: strength-4 covering arrays improving the best known sizes*, 2026, doi:10.5281/zenodo.23741408. See also `CITATION.cff`.

## Licence

Arrays: CC BY 4.0. Code: MIT. See `LICENSE`.
