# Covering array records (strengths 3 and 4)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23741408.svg)](https://doi.org/10.5281/zenodo.23741408)

Explicit covering arrays that improve the best known sizes for 31
parameter sets. All 39 arrays in the table below improved the live covering
array table
[chrispbsmolen/ca-tools](https://github.com/chrispbsmolen/ca-tools)
at commit 392d8c9; the eight marked as overtaken no longer hold the best known
size (see below). The note describing them is at
[doi:10.5281/zenodo.23842218](https://doi.org/10.5281/zenodo.23842218).
Versions 1.0 to 1.3 hold strength-4 arrays. Version 1.4 adds 28 parameter
sets (13 at strength 3, 15 at strength 4) and improves three of our own
entries: CA(1006; 4, 8, 5), CA(1957; 4, 8, 6) and CA(3240; 4, 12, 6). All
39 were re-checked on 2026-10-30 against ca-tools commit 392d8c9.

| Parameters | Previous best | Previous construction | New N | File |
|---|---|---|---|---|
| CA(N; 3, 32, 5) | 324 | Raaphorst-Moura-Stevens extension (Bergeron 2026): m columns added to the RMS CA(2q^3 - 1; 3, q^2 + q + 1, q), q a prime power, using L = max(q - 2, ceil(m (q - 2) / q)) layers of q^2 rows; N = 2q^3 - 1 + q^2 L for m = 1 and N = 2q^3 + q^2 L for m >= 2 | **323** | `arrays/CA_323_3_32_5.txt` |
| CA(N; 3, 58, 7) | 930 | Raaphorst-Moura-Stevens extension (Bergeron 2026): m columns added to the RMS CA(2q^3 - 1; 3, q^2 + q + 1, q), q a prime power, using L = max(q - 2, ceil(m (q - 2) / q)) layers of q^2 rows; N = 2q^3 - 1 + q^2 L for m = 1 and N = 2q^3 + q^2 L for m >= 2 | **929** | `arrays/CA_929_3_58_7.txt` |
| CA(N; 3, 74, 8) | 1407 | Raaphorst-Moura-Stevens extension (Bergeron 2026): m columns added to the RMS CA(2q^3 - 1; 3, q^2 + q + 1, q), q a prime power, using L = max(q - 2, ceil(m (q - 2) / q)) layers of q^2 rows; N = 2q^3 - 1 + q^2 L for m = 1 and N = 2q^3 + q^2 L for m >= 2 | **1406** | `arrays/CA_1406_3_74_8.txt` |
| CA(N; 3, 92, 9) | 2024 | Raaphorst-Moura-Stevens extension (Bergeron 2026): m columns added to the RMS CA(2q^3 - 1; 3, q^2 + q + 1, q), q a prime power, using L = max(q - 2, ceil(m (q - 2) / q)) layers of q^2 rows; N = 2q^3 - 1 + q^2 L for m = 1 and N = 2q^3 + q^2 L for m >= 2 | **2023** | `arrays/CA_2023_3_92_9.txt` |
| CA(N; 3, 8, 10) | 1180 | ordered design construction of Cohen, Colbourn and Ling (2008, Construction 2), strength 3, with a binary ingredient CA(12; 3, 8, 2) with two constant rows, as implemented in CAs (ODbasedCA), the route chosen by CAs::bestCA | **1170** | `arrays/CA_1170_3_8_10.txt` |
| CA(N; 3, 16, 10) | 2215 | Chateauneuf-Kreher doubling (Chateauneuf and Kreher 2002), strength 3, as implemented in CAs (CK_doublingCA), the route chosen by CAs::bestCA | **2205** | `arrays/CA_2205_3_16_10.txt` |
| CA(N; 3, 134, 11) | 3750 | Raaphorst-Moura-Stevens extension (Bergeron 2026): m columns added to the RMS CA(2q^3 - 1; 3, q^2 + q + 1, q), q a prime power, using L = max(q - 2, ceil(m (q - 2) / q)) layers of q^2 rows; N = 2q^3 - 1 + q^2 L for m = 1 and N = 2q^3 + q^2 L for m >= 2 | **3749** | `arrays/CA_3749_3_134_11.txt` |
| CA(N; 3, 8, 12) | 1992 | ordered design construction of Cohen, Colbourn and Ling (2008, Construction 2), strength 3, with a binary ingredient CA(12; 3, 8, 2) with two constant rows, as implemented in CAs (ODbasedCA), the route chosen by CAs::bestCA | **1980** | `arrays/CA_1980_3_8_12.txt` |
| CA(N; 3, 184, 13) | 6252 | Raaphorst-Moura-Stevens extension (Bergeron 2026): m columns added to the RMS CA(2q^3 - 1; 3, q^2 + q + 1, q), q a prime power, using L = max(q - 2, ceil(m (q - 2) / q)) layers of q^2 rows; N = 2q^3 - 1 + q^2 L for m = 1 and N = 2q^3 + q^2 L for m >= 2 | **6251** | `arrays/CA_6251_3_184_13.txt` |
| CA(N; 3, 8, 14) | 3108 | ordered design construction of Cohen, Colbourn and Ling (2008, Construction 2), strength 3, with a binary ingredient CA(12; 3, 8, 2) with two constant rows, as implemented in CAs (ODbasedCA), the route chosen by CAs::bestCA | **3094** | `arrays/CA_3094_3_8_14.txt` |
| CA(N; 3, 8, 18) | 6444 | ordered design construction of Cohen, Colbourn and Ling (2008, Construction 2), strength 3, with a binary ingredient CA(12; 3, 8, 2) with two constant rows, as implemented in CAs (ODbasedCA), the route chosen by CAs::bestCA | **6426** | `arrays/CA_6426_3_8_18.txt` |
| CA(N; 3, 8, 20) | 8760 | ordered design construction of Cohen, Colbourn and Ling (2008, Construction 2), strength 3, with a binary ingredient CA(12; 3, 8, 2) with two constant rows, as implemented in CAs (ODbasedCA), the route chosen by CAs::bestCA | **8740** | `arrays/CA_8740_3_8_20.txt` |
| CA(N; 3, 8, 24) | 14928 | ordered design construction of Cohen, Colbourn and Ling (2008, Construction 2), strength 3, with a binary ingredient CA(12; 3, 8, 2) with two constant rows, as implemented in CAs (ODbasedCA), the route chosen by CAs::bestCA | **14904** | `arrays/CA_14904_3_8_24.txt` |
| CA(N; 4, 100, 3) | 645 | Density Mult Cyclic (Colbourn) | **633** | `arrays/CA_633_4_100_3.txt` |
| CA(N; 4, 105, 3) | 657 | Density Mult Cyclic (Colbourn) | **644** | `arrays/CA_644_4_105_3.txt` |
| CA(N; 4, 7, 5) | 910 | cross-sum (CKRS) | **825** | `arrays/CA_825_4_7_5.txt` |
| CA(N; 4, 8, 5) | 1125 | derived from strength 5 (table entry k = 11) | **1006** | `arrays/CA_1006_4_8_5.txt` |
| CA(N; 4, 9, 5) | 1125 | Derive from strength 5 (table entry k = 11) | **1008** | `arrays/CA_1008_4_9_5.txt` |
| CA(N; 4, 7, 6) | 1844 | orthogonal array, fuse, postop (NCK) | **1832** | `arrays/CA_1832_4_7_6.txt` |
| CA(N; 4, 8, 6) | 1975 | orthogonal array, fuse, postop (NCK) | **1957** | `arrays/CA_1957_4_8_6.txt` |
| CA(N; 4, 9, 6) | 2768 | Density (Bryce-Colbourn) postop NCK | **2518** | `arrays/CA_2518_4_9_6.txt` |
| CA(N; 4, 10, 6) | 3208 | Density (Bryce-Colbourn) postop NCK | **2736** | `arrays/CA_2736_4_10_6.txt` |
| CA(N; 4, 11, 6) | 3671 | Density Cyclic (Colbourn) postop NCK | **3000** | `arrays/CA_3000_4_11_6.txt` |
| CA(N; 4, 12, 6) | 3944 | Density (Bryce-Colbourn), postop (NCK) | **3240** | `arrays/CA_3240_4_12_6.txt` |
| CA(N; 4, 13, 6) | 4027 | SCPHF Tabu Search (Walker-Colbourn) fuse postop NCK | **3606** | `arrays/CA_3606_4_13_6.txt` |
| CA(N; 4, 14, 6) | 4080 | SCPHF Tabu Search (Walker-Colbourn) fuse postop NCK | **3846** | `arrays/CA_3846_4_14_6.txt` |
| CA(N; 4, 15, 6) | 4132 | SCPHF Tabu Search (Walker-Colbourn) fuse postop NCK | 4116 (overtaken: 4101) | `arrays/CA_4116_4_15_6.txt` |
| CA(N; 4, 16, 6) | 4203 | SCPHF Tabu Search (Walker-Colbourn), fuse, postop (NCK) | 4202 (overtaken: 4111) | `arrays/CA_4202_4_16_6.txt` |
| CA(N; 4, 19, 6) | 5422 | Density Cyclic (Colbourn) postop NCK | 5046 (overtaken: 4293) | `arrays/CA_5046_4_19_6.txt` |
| CA(N; 4, 20, 6) | 5584 | Density Cyclic (Colbourn) postop NCK | 5346 (overtaken: 4414) | `arrays/CA_5346_4_20_6.txt` |
| CA(N; 4, 21, 6) | 5743 | Density Cyclic (Colbourn) postop NCK | **5487** | `arrays/CA_5487_4_21_6.txt` |
| CA(N; 4, 22, 6) | 5895 | Density Cyclic (Colbourn) postop NCK | **5678** | `arrays/CA_5678_4_22_6.txt` |
| CA(N; 4, 23, 6) | 6062 | Density Cyclic (Colbourn) postop NCK | 6051 (overtaken: 6023) | `arrays/CA_6051_4_23_6.txt` |
| CA(N; 4, 26, 6) | 6456 | Density Cyclic (Colbourn) | **6445** | `arrays/CA_6445_4_26_6.txt` |
| CA(N; 4, 27, 6) | 6606 | Density Cyclic (Colbourn) | 6581 (overtaken: 6568) | `arrays/CA_6581_4_27_6.txt` |
| CA(N; 4, 29, 6) | 6851 | 3-Restricted SCPHF RE (CL) fuse (table entry k = 30) | 6846 (overtaken: 6807) | `arrays/CA_6846_4_29_6.txt` |
| CA(N; 4, 9, 7) | 3569 | orthogonal array, fuse, postop (NCK) | 3568 (overtaken: 3557) | `arrays/CA_3568_4_9_7.txt` |
| CA(N; 4, 6, 10) | 11957 | Augment OA | **11710** | `arrays/CA_11710_4_6_10.txt` |
| CA(N; 4, 6, 12) | 24201 | Augment OA | **23772** | `arrays/CA_23772_4_6_12.txt` |

`arrays/CA_636_4_100_3.txt` (version 1.0) is kept; CA(633; 4, 100, 3)
was obtained from it by a second search. `arrays/CA_908_4_7_5.txt`
(version 1.0) is kept as well; it is superseded by CA(825; 4, 7, 5).
Likewise the earlier arrays CA(1838; 4, 7, 6), CA(1833; 4, 7, 6),
CA(1966; 4, 8, 6), CA(1963; 4, 8, 6), CA(6453; 4, 26, 6) and
CA(6451; 4, 26, 6) are kept, each being the starting array of the next
search (see Method). CA(1050; 4, 8, 5) and CA(3942; 4, 12, 6) are kept and
superseded by CA(1006; 4, 8, 5) and CA(3240; 4, 12, 6).

Eight arrays are overtaken by the deposit of Bai XunDi, published on 30 October
2026 ([doi:10.5281/zenodo.23787055](https://doi.org/10.5281/zenodo.23787055)).
We missed it when we first claimed them. Its eight smaller arrays pass both of
our checkers. Our files stay here for reproducibility and are no longer claimed
as improvements.

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

### Version 1.4

**Orbit search (14 arrays, v = 5, 6, 10, 12).** Rows are unions of orbits of
PGL(2, q) acting on the q + 1 symbols, chosen by local search over starter
rows, then a descent through subgroups of PGL(2, q) and a row-level repair.
This starter-vector idea is known: Chateauneuf, Colbourn and Kreher at
strength 3, Meagher and Stevens at strength 2, and Maity, Akhtar,
Chandrasekharan and Colbourn (arXiv:1509.03547) at strength 4, with bounds
improved for 3 symbols. The new part here is the parameters reached, with the subgroup descent and the
repair.

`generators/CA_2736_4_10_6.c` and `generators/CA_3240_4_12_6.c` come from the
automated search loop (programs written by Claude Opus). `generators/orbit.c`
is the same search with q and k as parameters; it produced CA(11710; 4, 6, 10),
CA(23772; 4, 6, 12), CA(1006; 4, 8, 5), CA(1008; 4, 9, 5) and CA(N; 4, k, 6) for
k = 13, 14, 15, 19, 20, 21, 22, 23, with no LLM at run time:

    gcc -O3 -march=native -mcmodel=medium -DQ=5 -DK=19 -DNS=3876 -DC3=816 -DMAXN=8000 -w -o orbit generators/orbit.c -lm
    OUT=out.txt TABLE_N=5422 SEED=1 REPAIR=300 SCALE=1 ./orbit

`NS` is C(k, 4), `C3` is C(k-1, 3), `MAXN` is the larger of 8000 and 1.3 times
the table value, `SCALE` is (v/12)^4 rounded to one decimal
with a minimum of 1. Arrays finished before commit 4c01e16 of the search used
a version of `orbit.c` that differs only by the size of the group table.

**Post-optimisation by annealing (3 arrays).** `generators/CA_1957_4_8_6.c`,
`generators/CA_2518_4_9_6.c` and `generators/CA_3000_4_11_6.c` delete a run and repair
the holes by simulated annealing, in the same automated loop as below (program written
by Claude Opus). CA(1957; 4, 8, 6) continues from CA(1963; 4, 8, 6) through 1960, 1959
and 1958 rows, each step with its own program and array here. The program behind
CA(2518; 4, 9, 6) and CA(3000; 4, 11, 6) was written for k = 12 and transposed to these
sizes with no LLM at run time; they start from our orbit arrays CA(2520; 4, 9, 6) and
CA(3031; 4, 11, 6), kept in `arrays/` with their programs in `generators/`.

**Classic post-optimisation (14 arrays, all 13 at strength 3 and CA(6846; 4, 29, 6)).**
`generators/postop.c` with no LLM at run time, on the record array of the table:

    ./postop CA_1180_3_8_10.txt 3 10 300 out.txt

CA(6846; 4, 29, 6) starts from CA(6852; 4, 29, 6) of `aadwyer/CA_Database`, with the
table value 6851 coming from its k = 30 entry.

Prior art, checked on 2026-10-30: the ca-tools table, `aadwyer/CA_Database`,
the 100 most recent arXiv papers and the 25 most recent Zenodo records matching
"covering array", plus a targeted search on PGL and strength-4 constructions. None
gives these sizes. Colbourn's online table was compared through its 2024-11 snapshot
held in ca-tools, which is the source of most previous values above. A smaller
unpublished or very recent array remains possible; an issue on this repository is welcome.

### Version 1.3

Three arrays continue the post-optimisation searches below, each
search starting from our previous best array instead of the table
array (`BEST_CANDIDATE` set to it):

- CA(N; 4, 7, 6): 1838 (version 1.0) -> 1833 (`generators/CA_1833_4_7_6.c`) -> 1832 (`generators/CA_1832_4_7_6.c`)
- CA(N; 4, 8, 6): 1966 (version 1.0) -> 1963 (`generators/CA_1963_4_8_6.c`)
- CA(N; 4, 26, 6): 6453 (version 1.1) -> 6451 (`generators/CA_6451_4_26_6.c`) -> 6445 (`generators/CA_6445_4_26_6.c`)

These programs were written by Claude Opus inside the same automated
search loop; each ran for at most one hour on one CPU core.

CA(3568; 4, 9, 7) and CA(4202; 4, 16, 6) involve no LLM at run time.
They come from `generators/postop.c`, a classic delete-and-patch
post-optimiser used as a baseline: pick a run with few uniquely
covered 4-tuples, delete it, and re-home each lost tuple by
reassigning entries of another run, a change being kept only if it
uncovers nothing. Each array came out of a 5-minute run on the record
array of ca-tools (from `aadwyer/CA_Database`):

    gcc -O3 -march=native -o postop generators/postop.c -lm
    ./postop CA_3569_4_9_7.txt 4 7 300 out.txt

Checked on 2026-10-30: the live ca-tools table still gives 1844, 1975,
3569, 4203 and 6456 for these five parameter sets.

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
