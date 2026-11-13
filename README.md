# Covering array records (strengths 3 to 6)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23741408.svg)](https://doi.org/10.5281/zenodo.23741408)

Explicit covering arrays that improve the best known sizes for 104
parameter sets. All 109 parameter sets in the table below improved the live
covering array table
[chrispbsmolen/ca-tools](https://github.com/chrispbsmolen/ca-tools)
when they were added (commit 392d8c9 for the first 39, commit 6c1e11c for the
eight of version 1.6, commit 7b76f80 for version 1.7, commit 84933d2 for version 1.8, commit 8a30577
for version 1.9, commit 69c708e for version 1.10); the five marked as overtaken no
longer hold the best known size (see below).
Version 1.10 adds three parameter sets, CA(559529; 4, 48, 20), CA(559570; 4, 48, 21)
and CA(559613; 4, 48, 22), and improves eight of our own entries (see Method). The
largest gains are CA(291; 3, 32, 5) and CA(884; 3, 58, 7) (our previous 323 and 929,
previous best 324 and 930), from search programs written by Claude Opus, and
CA(77835; 4, 31, 12), CA(36435; 4, 29, 10) and CA(42436; 5, 21, 6) by post-optimisation
(0.5 % to 0.8 % below our previous entries).
Version 1.9 adds five parameter sets and improves 41 of our own entries (see Method).
The largest gains are on the line CA(N; 4, k, 6): CA(6871; 4, 38, 6) and
CA(7195; 4, 39, 6) (previous best 7812 and 7902), from search programs written by
Claude Opus, and k = 26 to 37 by 5.4 % to 11.3 % over the previous best (for example
CA(6335; 4, 33, 6), previous best 7145), from post-optimisation of our CA(6445; 4, 37, 6),
CA(6446; 4, 37, 6) and CA(6254; 4, 29, 6). Among the five new parameter sets, CA(36615; 4, 29, 10) improves
41007 by 10.7 %.
Version 1.8 adds 21 parameter sets, among them the first three at strength 6, and
improves 28 of our own entries (see Method). The largest gains are CA(3516; 4, 14, 6)
and CA(3786; 4, 15, 6) (previous best 4080 and 4132), from the PGL(2, 5) family of
version 1.7 with one fixed column, CA(6254; 4, k, 6) for k = 26 to 29 (previous best
6456 to 6851) and CA(35983; 4, 28, 10) (39954). The note at
[doi:10.5281/zenodo.23842218](https://doi.org/10.5281/zenodo.23842218)
describes the arrays of version 1.4.
Version 1.7 adds 33 parameter sets and improves 13 of our own entries (see Method).
The largest gains are on the line CA(N; 4, k, 6): k = 28 to 39 by 3.4 % to 10.5 %
(CA(6445; 4, 37, 6), previous best 7199), and k = 11, 12, 15, 16 and 21 to 25 by a
family of arrays developed under PGL(2, 5) and a cyclic shift of the columns.
Version 1.6 adds eight parameter sets: CA(954; 3, 12, 8), CA(1813; 3, 13, 10),
CA(3057; 3, 15, 12), CA(5522; 3, 19, 14), CA(5870; 3, 19, 15), CA(6947; 4, 35, 6),
CA(7037; 4, 36, 6) and CA(2646; 5, 13, 4).
Version 1.5 improves four entries: CA(2616; 4, 10, 6), CA(3216; 4, 12, 6),
CA(6559; 4, 27, 6) and CA(3556; 4, 9, 7).
Versions 1.0 to 1.3 hold strength-4 arrays. Version 1.4 adds 28 parameter
sets (13 at strength 3, 15 at strength 4) and improves three of our own
entries: CA(1006; 4, 8, 5), CA(1957; 4, 8, 6) and CA(3240; 4, 12, 6). All
39 were re-checked on 2026-10-30 against ca-tools commit 392d8c9.

| Parameters | Previous best | Previous construction | New N | File |
|---|---|---|---|---|
| CA(N; 3, 34, 3) | 86 | Bai XunDi (doi:10.5281/zenodo.23787055): column extension + tabu repair and row-deletion post-optimisation starting from CA(87;3,39,3) = Kianoosh Shokri, Lucia Moura, Thm 5.2.8 (q-plication), q=3 (table entry k = 39) | **82** | `arrays/CA_82_3_34_3.txt` |
| CA(N; 3, 32, 5) | 324 | Raaphorst-Moura-Stevens extension (Bergeron 2026): m columns added to the RMS CA(2q^3 - 1; 3, q^2 + q + 1, q), q a prime power, using L = max(q - 2, ceil(m (q - 2) / q)) layers of q^2 rows; N = 2q^3 - 1 + q^2 L for m = 1 and N = 2q^3 + q^2 L for m >= 2 | **291** | `arrays/CA_291_3_32_5.txt` |
| CA(N; 3, 11, 6) | 418 | Grey (doi:10.5281/zenodo.23877721): Iterative row deletion with Nurmela-style tabu single-cell post-optimization from a 426-row PSL(2,5) group-orbit union on P^1(F_5) (1 constant 6-row orbit + 5 generic 60-row orbits + 4 binary-stabilizer {0,inf}^11 30-row orbits, repaired to a 426-row covering array via single-cell tabu moves) | 416 (overtaken: 411) | `arrays/CA_416_3_11_6.txt` |
| CA(N; 3, 58, 7) | 930 | Raaphorst-Moura-Stevens extension (Bergeron 2026): m columns added to the RMS CA(2q^3 - 1; 3, q^2 + q + 1, q), q a prime power, using L = max(q - 2, ceil(m (q - 2) / q)) layers of q^2 rows; N = 2q^3 - 1 + q^2 L for m = 1 and N = 2q^3 + q^2 L for m >= 2 | **884** | `arrays/CA_884_3_58_7.txt` |
| CA(N; 3, 12, 8) | 955 | Bai XunDi (doi:10.5281/zenodo.23787055): OA(q^3; 3, q+1, q) with evolved extension columns, greedy patch rows, row deletion and tabu | **954** | `arrays/CA_954_3_12_8.txt` |
| CA(N; 3, 74, 8) | 1407 | Raaphorst-Moura-Stevens extension (Bergeron 2026): m columns added to the RMS CA(2q^3 - 1; 3, q^2 + q + 1, q), q a prime power, using L = max(q - 2, ceil(m (q - 2) / q)) layers of q^2 rows; N = 2q^3 - 1 + q^2 L for m = 1 and N = 2q^3 + q^2 L for m >= 2 | **1406** | `arrays/CA_1406_3_74_8.txt` |
| CA(N; 3, 92, 9) | 2024 | Raaphorst-Moura-Stevens extension (Bergeron 2026): m columns added to the RMS CA(2q^3 - 1; 3, q^2 + q + 1, q), q a prime power, using L = max(q - 2, ceil(m (q - 2) / q)) layers of q^2 rows; N = 2q^3 - 1 + q^2 L for m = 1 and N = 2q^3 + q^2 L for m >= 2 | **2023** | `arrays/CA_2023_3_92_9.txt` |
| CA(N; 3, 8, 10) | 1180 | ordered design construction of Cohen, Colbourn and Ling (2008, Construction 2), strength 3, with a binary ingredient CA(12; 3, 8, 2) with two constant rows, as implemented in CAs (ODbasedCA), the route chosen by CAs::bestCA | **1170** | `arrays/CA_1170_3_8_10.txt` |
| CA(N; 3, 13, 10) | 1848 | Bai XunDi (doi:10.5281/zenodo.23787055): poly2 + pencil cover, q = 11 | **1807** | `arrays/CA_1807_3_13_10.txt` |
| CA(N; 3, 16, 10) | 2215 | Chateauneuf-Kreher doubling (Chateauneuf and Kreher 2002), strength 3, as implemented in CAs (CK_doublingCA), the route chosen by CAs::bestCA | 2205 (overtaken: 2178) | `arrays/CA_2205_3_16_10.txt` |
| CA(N; 3, 134, 11) | 3750 | Raaphorst-Moura-Stevens extension (Bergeron 2026): m columns added to the RMS CA(2q^3 - 1; 3, q^2 + q + 1, q), q a prime power, using L = max(q - 2, ceil(m (q - 2) / q)) layers of q^2 rows; N = 2q^3 - 1 + q^2 L for m = 1 and N = 2q^3 + q^2 L for m >= 2 | **3749** | `arrays/CA_3749_3_134_11.txt` |
| CA(N; 3, 8, 12) | 1992 | ordered design construction of Cohen, Colbourn and Ling (2008, Construction 2), strength 3, with a binary ingredient CA(12; 3, 8, 2) with two constant rows, as implemented in CAs (ODbasedCA), the route chosen by CAs::bestCA | **1980** | `arrays/CA_1980_3_8_12.txt` |
| CA(N; 3, 15, 12) | 3098 | Bai XunDi (doi:10.5281/zenodo.23787055): poly2 + pencil cover, q = 13 | **2955** | `arrays/CA_2955_3_15_12.txt` |
| CA(N; 3, 184, 13) | 6252 | Raaphorst-Moura-Stevens extension (Bergeron 2026): m columns added to the RMS CA(2q^3 - 1; 3, q^2 + q + 1, q), q a prime power, using L = max(q - 2, ceil(m (q - 2) / q)) layers of q^2 rows; N = 2q^3 - 1 + q^2 L for m = 1 and N = 2q^3 + q^2 L for m >= 2 | **6251** | `arrays/CA_6251_3_184_13.txt` |
| CA(N; 3, 8, 14) | 3108 | ordered design construction of Cohen, Colbourn and Ling (2008, Construction 2), strength 3, with a binary ingredient CA(12; 3, 8, 2) with two constant rows, as implemented in CAs (ODbasedCA), the route chosen by CAs::bestCA | **3094** | `arrays/CA_3094_3_8_14.txt` |
| CA(N; 3, 19, 14) | 5568 | hyperoval-poly2, alphabet-restricted, q = 16, as entered in the table with the arrays of Bai XunDi (doi:10.5281/zenodo.23787055) | **5255** | `arrays/CA_5255_3_19_14.txt` |
| CA(N; 3, 19, 15) | 5879 | Bai XunDi (doi:10.5281/zenodo.23787055): hyperoval-poly2 with alphabet restriction and fusion | **5704** | `arrays/CA_5704_3_19_15.txt` |
| CA(N; 3, 21, 17) | 9476 | Bai XunDi (doi:10.5281/zenodo.23787055): poly2 + pencil cover (paper A Thm 5.1, Prop 5.2) with Lemma F' on rows (a)+(c) | **9415** | `arrays/CA_9415_3_21_17.txt` |
| CA(N; 3, 8, 18) | 6444 | ordered design construction of Cohen, Colbourn and Ling (2008, Construction 2), strength 3, with a binary ingredient CA(12; 3, 8, 2) with two constant rows, as implemented in CAs (ODbasedCA), the route chosen by CAs::bestCA | **6426** | `arrays/CA_6426_3_8_18.txt` |
| CA(N; 3, 21, 18) | 9884 | Bai XunDi (doi:10.5281/zenodo.23787055): poly2 + pencil cover (paper A Thm 5.1, Prop 5.2) with Lemma F' on rows (a)+(c) | **9869** | `arrays/CA_9869_3_21_18.txt` |
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
| CA(N; 4, 10, 6) | 3208 | Density (Bryce-Colbourn) postop NCK | **2616** | `arrays/CA_2616_4_10_6.txt` |
| CA(N; 4, 11, 6) | 3671 | Density Cyclic (Colbourn) postop NCK | **2976** | `arrays/CA_2976_4_11_6.txt` |
| CA(N; 4, 12, 6) | 3944 | Density (Bryce-Colbourn), postop (NCK) | **2976** | `arrays/CA_2976_4_12_6.txt` |
| CA(N; 4, 13, 6) | 4027 | SCPHF Tabu Search (Walker-Colbourn) fuse postop NCK | **3510** | `arrays/CA_3510_4_13_6.txt` |
| CA(N; 4, 14, 6) | 4080 | SCPHF Tabu Search (Walker-Colbourn) fuse postop NCK | **3516** | `arrays/CA_3516_4_14_6.txt` |
| CA(N; 4, 15, 6) | 4132 | SCPHF Tabu Search (Walker-Colbourn) fuse postop NCK | **3786** | `arrays/CA_3786_4_15_6.txt` |
| CA(N; 4, 16, 6) | 4203 | SCPHF Tabu Search (Walker-Colbourn), fuse, postop (NCK) | **4043** | `arrays/CA_4043_4_16_6.txt` |
| CA(N; 4, 19, 6) | 5422 | Density Cyclic (Colbourn) postop NCK | 5046 (overtaken: 4293) | `arrays/CA_5046_4_19_6.txt` |
| CA(N; 4, 20, 6) | 5584 | Density Cyclic (Colbourn) postop NCK | 5346 (overtaken: 4414) | `arrays/CA_5346_4_20_6.txt` |
| CA(N; 4, 21, 6) | 5743 | Density Cyclic (Colbourn) postop NCK | **5136** | `arrays/CA_5136_4_21_6.txt` |
| CA(N; 4, 22, 6) | 5895 | Density Cyclic (Colbourn) postop NCK | **5406** | `arrays/CA_5406_4_22_6.txt` |
| CA(N; 4, 23, 6) | 6062 | Density Cyclic (Colbourn) postop NCK | **5676** | `arrays/CA_5676_4_23_6.txt` |
| CA(N; 4, 24, 6) | 6147 | Bai XunDi (doi:10.5281/zenodo.23787055): row deletion + Nurmela-style tabu repair post-optimisation of CA(6176;4,24,6) | **5945** | `arrays/CA_5945_4_24_6.txt` |
| CA(N; 4, 25, 6) | 6302 | Bai XunDi (doi:10.5281/zenodo.23787055): row deletion + Nurmela-style tabu repair post-optimisation of CA(6325;4,25,6) | **6108** | `arrays/CA_6108_4_25_6.txt` |
| CA(N; 4, 26, 6) | 6456 | Density Cyclic (Colbourn) | **6108** | `arrays/CA_6108_4_26_6.txt` |
| CA(N; 4, 27, 6) | 6606 | Density Cyclic (Colbourn) | **6152** | `arrays/CA_6152_4_27_6.txt` |
| CA(N; 4, 28, 6) | 6684 | Bai XunDi (doi:10.5281/zenodo.23787055): row deletion + Nurmela-style tabu repair post-optimisation of CA(6714;4,28,6) | **6207** | `arrays/CA_6207_4_28_6.txt` |
| CA(N; 4, 29, 6) | 6851 | 3-Restricted SCPHF RE (CL) fuse (table entry k = 30) | **6246** | `arrays/CA_6246_4_29_6.txt` |
| CA(N; 4, 30, 6) | 6851 | 3-Restricted SCPHF RE (CL) fuse | **6275** | `arrays/CA_6275_4_30_6.txt` |
| CA(N; 4, 31, 6) | 6851 | symbol fusion of a 7-symbol entry (ca-tools derived.csv) | **6306** | `arrays/CA_6306_4_31_6.txt` |
| CA(N; 4, 32, 6) | 7023 | SCPHF Conditional Expectation (CLS) postop NCK | **6306** | `arrays/CA_6306_4_32_6.txt` |
| CA(N; 4, 33, 6) | 7145 | 2-Restricted SCPHF RE (CL) fuse | **6335** | `arrays/CA_6335_4_33_6.txt` |
| CA(N; 4, 34, 6) | 7185 | Bai XunDi (doi:10.5281/zenodo.23787055): near-CPHF over GF(7), fusion 7 -> 6, tabu repair, row deletion (table entry k = 36) | **6374** | `arrays/CA_6374_4_34_6.txt` |
| CA(N; 4, 35, 6) | 7185 | Bai XunDi (doi:10.5281/zenodo.23787055): near-CPHF over GF(7), fusion 7 -> 6, tabu repair, row deletion (table entry k = 36) | **6401** | `arrays/CA_6401_4_35_6.txt` |
| CA(N; 4, 36, 6) | 7185 | Bai XunDi (doi:10.5281/zenodo.23787055): near-CPHF over GF(7), fusion 7 -> 6, tabu repair, row deletion | **6429** | `arrays/CA_6429_4_36_6.txt` |
| CA(N; 4, 37, 6) | 7199 | Bai XunDi (doi:10.5281/zenodo.23787055): near-CPHF over GF(7) -> CA, symbol fusion 7->6, tabu repair, row-deletion post-optimisation | **6442** | `arrays/CA_6442_4_37_6.txt` |
| CA(N; 4, 38, 6) | 7812 | Density Cyclic (Colbourn) | **6871** | `arrays/CA_6871_4_38_6.txt` |
| CA(N; 4, 39, 6) | 7902 | Density Cyclic (Colbourn) | **7195** | `arrays/CA_7195_4_39_6.txt` |
| CA(N; 4, 9, 7) | 3569 | orthogonal array, fuse, postop (NCK) | **3556** | `arrays/CA_3556_4_9_7.txt` |
| CA(N; 4, 6, 10) | 11957 | Augment OA | **11710** | `arrays/CA_11710_4_6_10.txt` |
| CA(N; 4, 13, 10) | 20555 | Bai XunDi (doi:10.5281/zenodo.23787055): degree-2 polynomial construction (poly2) with a pencil cover of AG(2,q), CA(20562; 4, 13, 10); then dedupe + optimal hitting-set pruning (CP-SAT) of CA(20562; 4, 13, 10) (R=6 redundant rows) | **20395** | `arrays/CA_20395_4_13_10.txt` |
| CA(N; 4, 21, 10) | 28951 | Bai XunDi (doi:10.5281/zenodo.23787055): Sherwood CPHF(2;k,q,4) -> CA, symbol fusion, wildcard row-deletion post-optimisation (table entry k = 25) | **28754** | `arrays/CA_28754_4_21_10.txt` |
| CA(N; 4, 22, 10) | 28951 | Bai XunDi (doi:10.5281/zenodo.23787055): Sherwood CPHF(2;k,q,4) -> CA, symbol fusion, wildcard row-deletion post-optimisation (table entry k = 25) | **28824** | `arrays/CA_28824_4_22_10.txt` |
| CA(N; 4, 23, 10) | 28951 | Bai XunDi (doi:10.5281/zenodo.23787055): Sherwood CPHF(2;k,q,4) -> CA, symbol fusion, wildcard row-deletion post-optimisation (table entry k = 25) | **28877** | `arrays/CA_28877_4_23_10.txt` |
| CA(N; 4, 24, 10) | 28951 | Bai XunDi (doi:10.5281/zenodo.23787055): Sherwood CPHF(2;k,q,4) -> CA, symbol fusion, wildcard row-deletion post-optimisation (table entry k = 25) | **28914** | `arrays/CA_28914_4_24_10.txt` |
| CA(N; 4, 25, 10) | 28951 | Bai XunDi (doi:10.5281/zenodo.23787055): Sherwood CPHF(2;k,q,4) -> CA, symbol fusion, wildcard row-deletion post-optimisation | **28935** | `arrays/CA_28935_4_25_10.txt` |
| CA(N; 4, 26, 10) | 29269 | Bai XunDi (doi:10.5281/zenodo.23787055): Sherwood CPHF(2;k,q,4) (Singer-cycle seeds + tabu column extension), CPHF->CA expansion | **29210** | `arrays/CA_29210_4_26_10.txt` |
| CA(N; 4, 27, 10) | 36860 | Bai XunDi (doi:10.5281/zenodo.23787055): CPHF tabu + residual block CA(41371; 4, 29, 11); symbol fusion CKRS 2010 Lemma 3.1 (Colbourn, Keri, Rivas Soriano, Schlage-Puchta), CAN(t,k,v-1) <= CAN(t,k,v) - 2; column drop; drop columns [27, 28] + | **35188** | `arrays/CA_35188_4_27_10.txt` |
| CA(N; 4, 28, 10) | 39954 | Bai XunDi (doi:10.5281/zenodo.23787055): CPHF tabu + residual block CA(41371; 4, 29, 11); symbol fusion CKRS 2010 Lemma 3.1 (Colbourn, Keri, Rivas Soriano, Schlage-Puchta), CAN(t,k,v-1) <= CAN(t,k,v) - 2; column drop; drop columns [27] + opt | **35963** | `arrays/CA_35963_4_28_10.txt` |
| CA(N; 4, 29, 10) | 41007 | Bai XunDi (doi:10.5281/zenodo.23787055): CPHF tabu + residual block construction CA(41371; 4, 29, 11); symbol fusion CKRS 2010 Lemma 3.1 (Colbourn, Keri, Rivas Soriano, Schlage-Puchta), CAN(t,k,v-1) <= CAN(t,k,v) - 2; prune -> CKRS fusion (fuse2, rows [12417, 11793]) x1 with optimal hitting-set pruning after each fusion, from CA(41371; 4, 29, 11) | **36435** | `arrays/CA_36435_4_29_10.txt` |
| CA(N; 4, 6, 12) | 24201 | Augment OA | **23772** | `arrays/CA_23772_4_6_12.txt` |
| CA(N; 4, 15, 12) | 40542 | Bai XunDi (doi:10.5281/zenodo.23787055): degree-2 polynomial construction (poly2) with a pencil cover of AG(2,q), CA(40550; 4, 15, 12); then dedupe + optimal hitting-set pruning (CP-SAT) of CA(40550; 4, 15, 12) (R=7 redundant rows) | **40476** | `arrays/CA_40476_4_15_12.txt` |
| CA(N; 4, 27, 12) | 56785 | Bai XunDi (doi:10.5281/zenodo.23787055): Sherwood CPHF -> CA, symbol fusion, wildcard row-deletion post-optimisation (table entry k = 29) | **56755** | `arrays/CA_56755_4_27_12.txt` |
| CA(N; 4, 28, 12) | 56785 | Bai XunDi (doi:10.5281/zenodo.23787055): Sherwood CPHF -> CA, symbol fusion, wildcard row-deletion post-optimisation (table entry k = 29) | **56768** | `arrays/CA_56768_4_28_12.txt` |
| CA(N; 4, 29, 12) | 56785 | Bai XunDi (doi:10.5281/zenodo.23787055): Sherwood CPHF -> CA, symbol fusion, wildcard row-deletion post-optimisation | **56776** | `arrays/CA_56776_4_29_12.txt` |
| CA(N; 4, 31, 12) | 79324 | Bai XunDi (doi:10.5281/zenodo.23787055): CPHF tabu + residual block CA(81445; 4, 32, 13); symbol fusion CKRS 2010 Lemma 3.1 (Colbourn, Keri, Rivas Soriano, Schlage-Puchta), CAN(t,k,v-1) <= CAN(t,k,v) - 2; column drop; drop columns [31] + optimal hitting-set pruning of CA(80989; 4, 32, 12) | **77835** | `arrays/CA_77835_4_31_12.txt` |
| CA(N; 4, 32, 12) | 80972 | Bai XunDi (doi:10.5281/zenodo.23787055): CPHF tabu + residual block construction CA(81445; 4, 32, 13); symbol fusion CKRS 2010 Lemma 3.1 (Colbourn, Keri, Rivas Soriano, Schlage-Puchta), CAN(t,k,v-1) <= CAN(t,k,v) - 2; prune -> CKRS fusion (fuse2, rows [18341, 67895]) x1 with optimal hitting-set pruning after each fusion, from CA(81445; 4, 32, 13) | **80780** | `arrays/CA_80780_4_32_12.txt` |
| CA(N; 4, 19, 14) | 106281 | Bai XunDi (doi:10.5281/zenodo.23787055): degree-2 polynomial construction (poly2) with a pencil cover of AG(2,q), CA(106293; 4, 19, 14); then dedupe + optimal hitting-set pruning (CP-SAT) of CA(106293; 4, 19, 14) (R=9 redundant rows) | **106023** | `arrays/CA_106023_4_19_14.txt` |
| CA(N; 4, 24, 14) | 126972 | 3-Restricted SCPHF RE (CL) fuse fuse | **126726** | `arrays/CA_126726_4_24_14.txt` |
| CA(N; 4, 25, 14) | 126972 | Bai XunDi (doi:10.5281/zenodo.23787055): elliptic two-block array, CAN(t,#E(F_q),q) <= 2q^t - q^(t-1), q=16; then symbol fusion CKRS 2010 Lemma 3.1 (Colbourn, Keri, Rivas Soriano, Schlage-Puchta), CAN(t,k,v-1) <= CAN(t,k,v) - 2, 2 step(s) | 126935 (overtaken: 126929) | `arrays/CA_126935_4_25_14.txt` |
| CA(N; 4, 19, 15) | 110386 | Bai XunDi (doi:10.5281/zenodo.23787055): degree-2 polynomial construction (poly2) with a pencil cover of AG(2,q), CA(110397; 4, 19, 15); then dedupe + optimal hitting-set pruning (CP-SAT) of CA(110397; 4, 19, 15) (R=9 redundant rows) | **110273** | `arrays/CA_110273_4_19_15.txt` |
| CA(N; 4, 24, 15) | 126974 | 3-Restricted SCPHF RE (CL) fuse | **126845** | `arrays/CA_126845_4_24_15.txt` |
| CA(N; 4, 25, 15) | 126974 | Bai XunDi (doi:10.5281/zenodo.23787055): elliptic two-block array, CAN(t,#E(F_q),q) <= 2q^t - q^(t-1), q=16; then symbol fusion CKRS 2010 Lemma 3.1 (Colbourn, Keri, Rivas Soriano, Schlage-Puchta), CAN(t,k,v-1) <= CAN(t,k,v) - 2, 1 step(s) | **126931** | `arrays/CA_126931_4_25_15.txt` |
| CA(N; 4, 19, 16) | 114495 | Bai XunDi (doi:10.5281/zenodo.23787055): degree-2 polynomial construction (poly2) with a pencil cover of AG(2,q), CA(114505; 4, 19, 16); then dedupe + optimal hitting-set pruning (CP-SAT) of CA(114505; 4, 19, 16) (R=9 redundant rows) | **114489** | `arrays/CA_114489_4_19_16.txt` |
| CA(N; 4, 21, 18) | 188415 | Bai XunDi (doi:10.5281/zenodo.23787055): degree-2 polynomial construction (poly2) with a pencil cover of AG(2,q), CA(188426; 4, 21, 18); then dedupe + optimal hitting-set pruning (CP-SAT) of CA(188426; 4, 21, 18) (R=10 redundant rows) | **188412** | `arrays/CA_188412_4_21_18.txt` |
| CA(N; 4, 40, 18) | 260621 | Bai XunDi (doi:10.5281/zenodo.23787055): Sherwood CPHF(2;k,q,4) (Singer-cycle seeds + tabu column extension), CPHF->CA expansion; symbol fusion of the own array CA(260623; 4, 40, 19) (1 step, 2 rows removed) | **260585** | `arrays/CA_260585_4_40_18.txt` |
| CA(N; 4, 48, 20) | 559653 | Bai XunDi (doi:10.5281/zenodo.23787055): cyclic / structured CPHF (2 rows) -> CA expansion; then symbol fusion CKRS 2010 Lemma 3.1, CAN(t,k,v-1) <= CAN(t,k,v) - 2, 3 steps | **559529** | `arrays/CA_559529_4_48_20.txt` |
| CA(N; 4, 48, 21) | 559655 | Bai XunDi (doi:10.5281/zenodo.23787055): cyclic / structured CPHF (2 rows) -> CA expansion; then symbol fusion CKRS 2010 Lemma 3.1, CAN(t,k,v-1) <= CAN(t,k,v) - 2, 2 steps | **559570** | `arrays/CA_559570_4_48_21.txt` |
| CA(N; 4, 48, 22) | 559657 | symbol fusion CKRS 2010 Lemma 3.1 (ca-tools derived entry) of the table's CA(559659; 4, 48, 23): two-row Sherwood CPHF(2; 2q+2, q, 4) from the Singer cycle of PG(3,q) (Bai XunDi, doi:10.5281/zenodo.23787053, Theorem 3.2 and Corollary 5.6), q = 23; affine Sherwood-Martirosyan-Colbourn expansion | **559613** | `arrays/CA_559613_4_48_22.txt` |
| CA(N; 5, 13, 4) | 2669 | Bai XunDi (doi:10.5281/zenodo.23787055): near-CPHF(2) + aligned linear residual block | **2595** | `arrays/CA_2595_5_13_4.txt` |
| CA(N; 5, 14, 5) | 8575 | Bai XunDi (doi:10.5281/zenodo.23787055): near-CPHF(2) + aligned linear residual block | **8353** | `arrays/CA_8353_5_14_5.txt` |
| CA(N; 5, 13, 6) | 28399 | Bai XunDi (doi:10.5281/zenodo.23787055): Sherwood CPHF -> CA, symbol fusion, wildcard row-deletion post-optimisation (table entry k = 15) | **27222** | `arrays/CA_27222_5_13_6.txt` |
| CA(N; 5, 14, 6) | 28399 | Bai XunDi (doi:10.5281/zenodo.23787055): Sherwood CPHF -> CA, symbol fusion, wildcard row-deletion post-optimisation (table entry k = 15) | **28082** | `arrays/CA_28082_5_14_6.txt` |
| CA(N; 5, 16, 6) | 30682 | Bai XunDi (doi:10.5281/zenodo.23787055): Sherwood CPHF -> CA, symbol fusion, wildcard row-deletion post-optimisation | **30264** | `arrays/CA_30264_5_16_6.txt` |
| CA(N; 5, 21, 6) | 43531 | Density Cyclic (Colbourn) postop NCK | **42436** | `arrays/CA_42436_5_21_6.txt` |
| CA(N; 5, 22, 6) | 43663 | Bai XunDi (doi:10.5281/zenodo.23787055): Sherwood CPHF -> CA, symbol fusion, wildcard row-deletion post-optimisation | **43440** | `arrays/CA_43440_5_22_6.txt` |
| CA(N; 5, 23, 6) | 44320 | Bai XunDi (doi:10.5281/zenodo.23787055): Sherwood CPHF -> CA, symbol fusion, wildcard row-deletion post-optimisation | **44155** | `arrays/CA_44155_5_23_6.txt` |
| CA(N; 5, 10, 8) | 49136 | Bai XunDi (doi:10.5281/zenodo.23787055): degree-2 polynomial construction + block cover / SA post-processing; then 16 redundant rows of CA(49152; 5, 10, 8) deleted | **49128** | `arrays/CA_49128_5_10_8.txt` |
| CA(N; 5, 13, 10) | 226208 | Bai XunDi (doi:10.5281/zenodo.23787055): degree-2 polynomial construction (poly2) with a pencil cover of AG(2,q), CA(226262; 5, 13, 10); then dedupe + optimal hitting-set pruning (CP-SAT) of CA(226262; 5, 13, 10) (R=53 redundant rows) | **224844** | `arrays/CA_224844_5_13_10.txt` |
| CA(N; 5, 15, 10) | 307433 | Bai XunDi (doi:10.5281/zenodo.23787055): elliptic two-block array, CAN(t,#E(F_q),q) <= 2q^t - q^(t-1), q=11; then symbol fusion CKRS 2010 Lemma 3.1 (Colbourn, Keri, Rivas Soriano, Schlage-Puchta), CAN(t,k,v-1) <= CAN(t,k,v) - 2, 1 step(s); t (table entry k = 18) | **306455** | `arrays/CA_306455_5_15_10.txt` |
| CA(N; 5, 16, 10) | 307433 | Bai XunDi (doi:10.5281/zenodo.23787055): elliptic two-block array, CAN(t,#E(F_q),q) <= 2q^t - q^(t-1), q=11; then symbol fusion CKRS 2010 Lemma 3.1 (Colbourn, Keri, Rivas Soriano, Schlage-Puchta), CAN(t,k,v-1) <= CAN(t,k,v) - 2, 1 step(s); t (table entry k = 18) | **307024** | `arrays/CA_307024_5_16_10.txt` |
| CA(N; 5, 17, 10) | 307433 | Bai XunDi (doi:10.5281/zenodo.23787055): elliptic two-block array, CAN(t,#E(F_q),q) <= 2q^t - q^(t-1), q=11; then symbol fusion CKRS 2010 Lemma 3.1 (Colbourn, Keri, Rivas Soriano, Schlage-Puchta), CAN(t,k,v-1) <= CAN(t,k,v) - 2, 1 step(s); t (table entry k = 18) | **307215** | `arrays/CA_307215_5_17_10.txt` |
| CA(N; 5, 18, 10) | 307433 | Bai XunDi (doi:10.5281/zenodo.23787055): elliptic two-block array, CAN(t,#E(F_q),q) <= 2q^t - q^(t-1), q=11; then symbol fusion CKRS 2010 Lemma 3.1 (Colbourn, Keri, Rivas Soriano, Schlage-Puchta), CAN(t,k,v-1) <= CAN(t,k,v) - 2, 1 step(s); then 26 redundant rows of CA(307459; 5, 18, 10) deleted | **307356** | `arrays/CA_307356_5_18_10.txt` |
| CA(N; 5, 15, 12) | 527262 | Bai XunDi (doi:10.5281/zenodo.23787055): degree-2 polynomial construction (poly2) with a pencil cover of AG(2,q), CA(527270; 5, 15, 12); then dedupe + optimal hitting-set pruning (CP-SAT) of CA(527270; 5, 15, 12) (R=7 redundant rows) | **527205** | `arrays/CA_527205_5_15_12.txt` |
| CA(N; 5, 16, 12) | 714023 | 4-Restricted SCPHF RE (CL) fuse | **713888** | `arrays/CA_713888_5_16_12.txt` |
| CA(N; 6, 14, 6) | 223223 | Bai XunDi (doi:10.5281/zenodo.23787055): Singer CPHF CA(235291; 6, 16, 7); symbol fusion CKRS 2010 Lemma 3.1 (Colbourn, Keri, Rivas Soriano, Schlage-Puchta), CAN(t,k,v-1) <= CAN(t,k,v) - 2; column drop; drop columns [0, 6] + optimal hitting-set pruning of CA(232255; 6, 16, 6) | **222949** | `arrays/CA_222949_6_14_6.txt` |
| CA(N; 6, 15, 6) | 229013 | Bai XunDi (doi:10.5281/zenodo.23787055): Singer CPHF CA(235291; 6, 16, 7); symbol fusion CKRS 2010 Lemma 3.1 (Colbourn, Keri, Rivas Soriano, Schlage-Puchta), CAN(t,k,v-1) <= CAN(t,k,v) - 2; column drop; drop columns [0] + optimal hitting-set pruning of CA(232255; 6, 16, 6) | **228773** | `arrays/CA_228773_6_15_6.txt` |
| CA(N; 6, 16, 6) | 232255 | Bai XunDi (doi:10.5281/zenodo.23787055): Singer CPHF array CA(235291; 6, 16, 7), symbol fusion (CKRS 2010 Lemma 3.1 (Colbourn, Keri, Rivas Soriano, Schlage-Puchta)) to CA(235289; 6, 16, 6); then dedupe + optimal hitting-set pruning (CP-SAT) of CA(235289; 6, 16, 6) (R=3040 redundant rows) | **232101** | `arrays/CA_232101_6_16_6.txt` |
| CA(N; 6, 10, 8) | 393088 | Bai XunDi (doi:10.5281/zenodo.23787055): OA-jump: RS-OA(q^t;t,q+1,q) + extra column(s) + patch blocks + residual rows (SA); then 128 redundant rows of CA(393216; 6, 10, 8) deleted | **393018** | `arrays/CA_393018_6_10_8.txt` |

`arrays/CA_636_4_100_3.txt` (version 1.0) is kept; CA(633; 4, 100, 3)
was obtained from it by a second search. `arrays/CA_908_4_7_5.txt`
(version 1.0) is kept as well; it is superseded by CA(825; 4, 7, 5).
Likewise the earlier arrays CA(1838; 4, 7, 6), CA(1833; 4, 7, 6),
CA(1966; 4, 8, 6), CA(1963; 4, 8, 6), CA(6453; 4, 26, 6) and
CA(6451; 4, 26, 6) are kept, each being the starting array of the next
search (see Method). CA(1050; 4, 8, 5), CA(3942; 4, 12, 6) and CA(3240; 4, 12, 6)
are kept and superseded by CA(1006; 4, 8, 5) and CA(3216; 4, 12, 6). CA(2736; 4, 10, 6),
CA(6581; 4, 27, 6) and CA(3568; 4, 9, 7) are kept and superseded by the arrays of
version 1.5.

Version 1.10 replaces earlier arrays of ours for eight parameter sets: CA(N; 3, 32, 5),
CA(N; 3, 58, 7), CA(N; 4, 29, 10), CA(N; 4, 31, 12), CA(N; 5, 13, 10), CA(N; 5, 16, 6),
CA(N; 5, 21, 6) and CA(N; 5, 23, 6); the earlier files stay in `arrays/`.

Version 1.9 replaces earlier arrays of ours for 41 parameter sets: CA(N; 3, 15, 12),
CA(N; 3, 19, 14), CA(N; 3, 19, 15), CA(N; 3, 21, 18), CA(N; 3, 34, 3), CA(N; 4, 16, 6),
CA(N; 4, k, 6) for k = 24 to 39, CA(N; 4, 25, 10), CA(N; 4, 27, 10), CA(N; 4, 28, 10),
CA(N; 4, k, 12) for k = 27, 28, 31 and 32, CA(N; 5, 13, 4), CA(N; 5, 13, 6),
CA(N; 5, 14, 5), CA(N; 5, 14, 6), CA(N; 5, 13, 10), CA(N; 5, k, 10) for k = 15 to 18,
CA(N; 5, 15, 12), CA(N; 6, 15, 6) and CA(N; 6, 16, 6); the earlier files stay in
`arrays/`. The files for CA(6108; 4, 25, 6) and CA(6306; 4, 31, 6) are the first k
columns of CA(6108; 4, 26, 6) and CA(6306; 4, 32, 6).

Version 1.8 replaces earlier arrays of ours for 28 parameter sets: CA(N; 3, 19, 14),
CA(N; 3, 21, 17), CA(N; 4, 13, 10), CA(N; 4, k, 10) for k = 21 to 24, 27 and 28,
CA(N; 4, 24, 14), CA(N; 4, 24, 15), CA(N; 4, k, 6) for k = 14, 15, 26 to 30, 38 and 39,
CA(N; 5, 13, 4), CA(N; 5, 14, 5), CA(N; 5, 13, 6), CA(N; 5, 14, 6), CA(N; 5, 16, 6) and
CA(N; 5, k, 10) for k = 15 to 17; the earlier files stay in `arrays/`. The files for
CA(6254; 4, k, 6) with k = 26, 27 and 28 are the first k columns of CA(6254; 4, 29, 6).

Version 1.7 replaces earlier arrays of ours for CA(N; 3, 13, 10), CA(N; 3, 15, 12),
CA(N; 4, k, 6) with k = 11, 12, 13, 21, 22, 35 and 36, and takes back k = 15, 16, 23 and 29,
which the deposit of Bai XunDi had overtaken; the earlier files stay in `arrays/`.
By deleting columns, CA(6445; 4, 37, 6) also gives CA(6445; 4, 27, 6), below our
CA(6559; 4, 27, 6) of version 1.5.

Six arrays were overtaken by the deposit of Bai XunDi, published on 30 October
2026 ([doi:10.5281/zenodo.23787055](https://doi.org/10.5281/zenodo.23787055)).
We missed it when we first claimed them. Its smaller arrays pass both of
our checkers. Our files stay here for reproducibility and are no longer claimed
as improvements. CA(416; 3, 11, 6) of version 1.8 is overtaken by CA(411; 3, 11, 6)
of Kevin Grey in the ca-tools table (doi:10.5281/zenodo.24096230, built from his
CA(418; 3, 11, 6)). Our CA(126935; 4, 25, 14) is overtaken by the derived entry
CA(126929; 4, 25, 14) of the table, obtained by symbol fusion (CKRS 2010, Lemma 3.1)
from our CA(126931; 4, 25, 15).
Two more were overtaken by the deposit of Bai XunDi until version 1.5:
its CA(6568; 4, 27, 6) and CA(3557; 4, 9, 7) are now above our CA(6559; 4, 27, 6)
and CA(3556; 4, 9, 7).

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

### Version 1.10

Checked on 2026-11-13, after the table of ca-tools at commit 69c708e was read
again: the 11 files added in this version pass `ca_verify` of ca-tools, the second
checker and a negative control (the array with one run removed is rejected).
The three post-optimised arrays CA(N; 4, 48, v), v = 20, 21, 22, each had one
redundant run, their first; the files here are those arrays without it, and they
pass all three checks themselves.

**Search programs written by Claude Opus (2 arrays).** CA(291; 3, 32, 5) and
CA(884; 3, 58, 7) come from programs written by Claude Opus in the automated search
loop (delete a run and repair, as in version 1.9). They are in `generators/v1.10/`,
named after the array they produced, and ran with no LLM at run time. Each reads its
starting array from the environment variable `BEST_CANDIDATE`: the smallest array the
lab held for that size, our CA(323; 3, 32, 5) and CA(929; 3, 58, 7) of earlier versions.

**Post-optimisation (9 arrays).** `generators/postop.c`, with no LLM at run time, as
in versions 1.8 and 1.9. The chains start from:
- our earlier entries CA(36615; 4, 29, 10), CA(78458; 4, 31, 12), CA(225184; 5, 13, 10),
  CA(30386; 5, 16, 6), CA(42655; 5, 21, 6) and CA(44281; 5, 23, 6), which descend from
  the record arrays of Bai XunDi (doi:10.5281/zenodo.23787055) (for CA(N; 5, 21, 6), his
  CA(43663; 5, 22, 6) with one column deleted);
- the arrays CA(559653; 4, 48, 20) and CA(559655; 4, 48, 21) of Bai XunDi
  (doi:10.5281/zenodo.23787055) for v = 20 and 21;
- for CA(N; 4, 48, 22), the table's CA(559659; 4, 48, 23) (two-row Sherwood CPHF from
  the Singer cycle of PG(3, 23), Bai XunDi, doi:10.5281/zenodo.23787053), fused to
  22 symbols (CKRS 2010, Lemma 3.1).
The gains on CA(N; 4, 48, v) are small (0.022 %, 0.015 % and 0.008 % for v = 20, 21 and 22);
we list these values as improved sizes only.

Prior art for the 11, checked on 2026-11-13: the ca-tools table at commit 69c708e with
stair-step, symbol fusion and derived entries, its open pull requests, the database of
Dwyer (`aadwyer/CA_Database`), and the Zenodo, GitHub and arXiv records matching
"covering array" or "covering perfect hash family". None gives a size at or below ours.

### Version 1.9

Checked on 2026-11-11, after the table of ca-tools at commit 8a30577 was read
again: the 46 files added in this version pass `ca_verify` of ca-tools, the second
checker and a negative control (the array with its first run removed is rejected).
The post-optimised CA(N; 4, 40, 18) had one redundant run, its first; the file here is
that array without it, 260585 runs, and it passes all three checks itself.

**Search programs written by Claude Opus (7 arrays).** CA(2955; 3, 15, 12),
CA(5704; 3, 19, 15), CA(82; 3, 34, 3), CA(4043; 4, 16, 6), CA(6871; 4, 38, 6),
CA(7195; 4, 39, 6) and CA(2595; 5, 13, 4) come from programs written by Claude Opus in
the automated search loop, one per parameter set, each a revision of the previous
program for that size (delete a run and repair, by local search with clause weighting,
simulated annealing or a merge of two runs). They are in `generators/v1.9/`, named
after the array they produced, and ran with no LLM at run time. Each reads its starting array from the environment variable `BEST_CANDIDATE`:
the smallest array the lab held for that size, which descends from our earlier entry
in the table (for CA(N; 4, 38, 6), from a one-column extension of
CA(6445; 4, 37, 6)).

**Post-optimisation (37 arrays).** `generators/postop.c`, with no LLM at run time, as
in version 1.8. The chains start from:
- the record arrays of Bai XunDi (doi:10.5281/zenodo.23787055), through our earlier
  entries where we had one, for 24 parameter sets: CA(N; 3, 21, 18), CA(N; 4, 21, 18),
  CA(N; 4, 24, 6), CA(N; 4, 25, 10), CA(N; 4, k, 10) for k = 27 to 29, CA(N; 4, k, 12)
  for k = 27, 28, 31 and 32, CA(N; 4, 40, 18), CA(N; 5, 13, 6), CA(N; 5, 14, 5),
  CA(N; 5, 14, 6), CA(N; 5, 13, 10), CA(N; 5, k, 10) for k = 15 to 18, CA(N; 5, 15, 12)
  and CA(N; 6, k, 6) for k = 14 to 16 (with columns deleted where the table entry has a
  larger k);
- the table's CA(714023; 5, 16, 12) (construction "4-Restricted SCPHF RE (CL) fuse")
  for CA(N; 5, 16, 12);
- our CA(6445; 4, 37, 6) of version 1.7 with columns deleted for k = 26, 27 and 32 to 37,
  our CA(6254; 4, 29, 6) of version 1.8 for k = 28 and 29, and our CA(6446; 4, 37, 6)
  with columns deleted for k = 30;
- for CA(N; 3, 19, 14), the 5261-run output of an Opus search program of the lab,
  which starts from our CA(5518; 3, 19, 14).
For the five new parameter sets the gains are small except CA(36615; 4, 29, 10)
(0.002 % to 0.12 % for the other four); we list these values as improved sizes only.

Prior art for the 46, checked on 2026-11-11: the ca-tools table at commit 8a30577 with
stair-step, symbol fusion and derived entries, its open pull requests, the database of
Dwyer (`aadwyer/CA_Database`), and the Zenodo, GitHub and arXiv records matching
"covering array" or "covering perfect hash family". None gives a size at or below ours.

### Version 1.8

Checked on 2026-11-08, after the table of ca-tools at commit 84933d2 was read
again: every array here passes `ca_verify` of ca-tools and a negative control (the
array with its first run removed is rejected), and the 50 files added in this version
pass the second checker as well.

**CA(3516; 4, 14, 6) and CA(3786; 4, 15, 6).** The family of version 1.7 with one
fixed column (f = 1): the 6 constant rows and the orbits of three starter rows under
PGL(2, 5) on the symbols times Z_(k-1) on the first k - 1 columns, 270 (k - 1) + 6 rows.
The starters were found by `generators/v1.8/family-c/starter.c`, a C port of the
search of version 1.7 (local search with clause weighting over the three starters),
written by Claude Opus and run with no LLM at run time: `run.py C4 14 10800 101` found
the first in 15 minutes and `run.py C4 15 10800 102` the second in 2 h 47 min, on one
core. The starters are in `generators/v1.8/family-c/*_starters.json`, and `common.py`
(the file of version 1.7) develops them into the arrays here.

**CA(6254; 4, 29, 6) and k = 26 to 28.** `generators/CA_6254_4_29_6.c`, written by
Claude Opus in the automated search loop, ran once for 31 minutes on one core. It
first takes the first 29 columns of a wider array held in the lab, here the 6446-run
array of the version 1.7 chain towards CA(6445; 4, 37, 6) (`arrays/CA_6446_4_37_6.txt`,
added in this version as the starting array), then deletes runs and
repairs as its parent program did. For k = 26 to 28 the files are the first k columns
of CA(6254; 4, 29, 6).

**Post-optimisation (43 arrays).** `generators/postop.c`, with no LLM at run time, run
repeatedly in passes of 5, 30 and 120 minutes on one core, each pass starting from the
smallest verified array we held for that parameter set or for a larger k with columns
deleted. For the 21 parameter sets that are new here, the chains start from the record
arrays of the table, credited above: those of Bai XunDi (doi:10.5281/zenodo.23787055)
for 19 of them, CA(418; 3, 11, 6) of Grey (doi:10.5281/zenodo.23877721) and, for
CA(N; 5, 21, 6), the CA(43663; 5, 22, 6) of Bai XunDi with a column deleted, the table
value 43531 being a cited entry with no array. The other chains continue from our
arrays of versions 1.6 and 1.7; CA(N; 4, 24, 14) descends from the CA(126972; 4, 25, 14)
of Bai XunDi with a column deleted. The gains over the previous best are small for the
new parameter sets (0.003 % to 2.3 %); we list these values as improved sizes only.

Prior art for the 49, checked on 2026-11-08: the ca-tools table at commit 84933d2 with
stair-step, symbol fusion and derived entries, its open pull requests, the database of
Dwyer (`aadwyer/CA_Database`), and the Zenodo, GitHub and arXiv records matching
"covering array" or "covering perfect hash family", among them the deposits of
Bai XunDi, Smolen (doi:10.5281/zenodo.23848283, whose strength-6 arrays have at most 5
symbols), Canli (doi:10.5281/zenodo.23839216) and Grey. None gives a size at or below ours.

### Version 1.7

Checked on 2026-11-04, after the table of ca-tools at commit 7b76f80 was read
again: every array here passes `ca_verify` of ca-tools (version 2) and a negative
control (the array with its first run removed is rejected).

**A family for CA(N; 4, k, 6) (k = 11, 12, 15, 16, 21 to 25).** Symbols are the
points of PG(1, 5) and the group is PGL(2, 5) (order 120) on the symbols times the
cyclic group Z_m on the first m = k - f columns, the last f columns being fixed. The
rows are the 6 constant rows and the orbits of three starter rows (two with at
least three distinct symbols, one binary). Each array has 270 m + 6 rows or
fewer: 2976 rows for k = 11 and 12, 4056 for k = 15 and 16, and 270 (k - 2) + 6 for
k = 21 to 25. The starters were found by a search written by Claude Opus
(`generators/v1.7/family/`, Python with NumPy, `c3.py`, `c4.py` and `c5.py` for no, one
and two fixed columns); the search found no starters for the other k it tried, which
decides nothing. This continues the starter-vector constructions of Maity, Akhtar,
Chandrasekharan and Colbourn (arXiv:1509.03547), which develop strength-4 arrays under
a group and give tables for three symbols; the structure was read off our CA(3510; 4, 13, 6),
invariant under a cyclic shift of its columns, and CA(3216; 4, 12, 6), invariant under
PGL(2, 5) on its symbols.

**CA(6445; 4, 37, 6) and k = 28 to 36.** `generators/CA_6445_4_37_6.c`, written by
Claude Opus in the automated search loop, deletes a run and repairs: before a run is
deleted, the 4-tuples it alone covers are moved to other runs where that is free, and
what remains is repaired by ejection chains. The chain of searches starts from the
post-optimisation of the record array CA(7199; 4, 37, 6) of Bai XunDi
(doi:10.5281/zenodo.23787055) and goes 7027, 6830, ..., 6446, 6445, each step a new
program in the same loop. For k = 28 to 36 the files are the first k columns of
CA(6445; 4, 37, 6).

**CA(7350; 4, 38, 6) and CA(7783; 4, 39, 6).** `generators/v1.7/extend/`, no LLM at run
time: `extend.c` adds one column to the best array for k - 1 (local search on the new
column, then greedy runs for what is still uncovered), and the program of an earlier
search for CA(N; 4, 34, 6), recompiled with K = k (`gen.c`), shrinks the result for 12
minutes. The chain starts from our CA(6830; 4, 37, 6). From k = 40 the same chain stays
above the table.

**CA(1807; 3, 13, 10), CA(2964; 3, 15, 12) and CA(3510; 4, 13, 6).** Programs written by
Claude Opus in the same loop (`generators/CA_1807_3_13_10.c`, `generators/CA_2964_3_15_12.c`,
`generators/CA_3510_4_13_6.c`), each run for at most one hour on one CPU core, starting from
our arrays of version 1.6 and 1.5; the first two descend from the arrays of Bai XunDi.

**Post-optimisation of the table arrays (22 parameter sets, v = 5, 6, 10, 12, 14, 15, 17, 18).**
`generators/postop.c` with no LLM at run time, on the record arrays of the table, run
repeatedly: CA(N; 3, 21, 17), CA(N; 3, 21, 18), CA(N; 4, 13, 10), CA(N; 4, k, 10) for k = 21 to 28,
CA(N; 4, 27, 12), CA(N; 4, 28, 12), CA(N; 4, 24, 14), CA(N; 4, 24, 15), CA(N; 5, 13, 6),
CA(N; 5, 14, 5), CA(N; 5, 14, 6), CA(N; 5, 16, 6) and CA(N; 5, k, 10) for k = 15 to 17. Most
starting arrays are those of Bai XunDi (doi:10.5281/zenodo.23787055), credited in the table
above; when the table gives a parameter set the value of a larger k, the start is that array
with columns deleted. The gains are small (0.002 % to 4.1 %), and post-optimisation of
published arrays is standard; we list these values as improved sizes only.

Prior art for the 46, checked on 2026-11-04: the ca-tools table at commit 7b76f80 with
stair-step, symbol fusion and derived entries (`derived.csv`), its open pull requests, the
database of Dwyer (`aadwyer/CA_Database`), the arXiv and Zenodo records matching "covering
array", and the deposits of Bai XunDi, Smolen (doi:10.5281/zenodo.23848283, whose alphabets
stop at 9 symbols) and Canli (doi:10.5281/zenodo.23839216, whose improved cells start at
k = 100 for strength 4 and k = 241 for strength 3). None gives a size at or below ours.

### Version 1.6

All eight come from `generators/postop.c` (delete a run and repair), with no LLM at
run time, applied to the record arrays of the table. Seven of these starting arrays
are from the deposit of Bai XunDi (doi:10.5281/zenodo.23787055), the eighth,
CA(5568; 3, 19, 14), from the table entry built the same way; each improvement is a
post-optimisation of his array, and we credit the construction to him. The table
gives CA(4, 35, 6) the value of its k = 36 entry, 7185, and the gains there are the
largest of this version (3.3 % and 2.1 %). CA(1813; 3, 13, 10) continues from our own
intermediate CA(1814; 3, 13, 10), itself obtained from his CA(1848; 3, 13, 10).

CA(2205; 3, 16, 10) of version 1.4 is overtaken: the table now derives
CA(2178; 3, 16, 10) by Chateauneuf-Kreher doubling.

Prior art for the eight, checked on 2026-11-02: the ca-tools table at commit 6c1e11c
with stair-step, symbol fusion and derived entries, its open pull requests (including
those of Arhan Canli up to number 29), and the deposits of Bai XunDi, Smolen
(doi:10.5281/zenodo.23848283), Canli (doi:10.5281/zenodo.23839216) and Grey
(doi:10.5281/zenodo.23877721). None gives a size at or below ours.

### Version 1.5

All four come from the automated search loop (programs written by Claude Opus,
one hour on one CPU core at most, each program in `generators/`).

- **CA(2616; 4, 10, 6)** (from CA(2736; 4, 10, 6)): orbit search under PGL(2, 5),
  with a warm start from the previous array and a search over the mix of generic,
  two-symbol and constant orbits instead of a descent one starter at a time.
- **CA(3216; 4, 12, 6)**: annealing over row orbits under PGL(2, 5) on the symbols,
  the short orbits being frozen as a prefix, then a descent through subgroups.
- **CA(6559; 4, 27, 6)** (6581 -> 6565 -> 6560 -> 6559, the intermediate arrays are
  here): delete a run and repair, with multi-cell moves and large-neighbourhood
  search on a small subset of rows.
- **CA(3556; 4, 9, 7)**: a union of 508 cosets of the all-ones code of length 9 over
  Z_7 (508 x 7 runs), which is the case d = 1 of the coset unions of version 1.2,
  with no completion or pruning needed.

Prior art, checked on 2026-11-01 for all 33 current entries, with stair-step and
symbol fusion taken into account: the ca-tools table at commit e8afc9a and its open
pull requests, and the papers of Bai XunDi (doi:10.5281/zenodo.23787051,
doi:10.5281/zenodo.23787053 and the data at doi:10.5281/zenodo.23787055), Smolen
(doi:10.5281/zenodo.23848283) and Canli (doi:10.5281/zenodo.23839216). None gives
a size at or below ours. A smaller unpublished or very recent array remains
possible; an issue on this repository is welcome.

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

T. N. Sitjar (tns-research), *Covering array records: covering arrays of strengths 3 to 6 improving the best known sizes*, 2026, doi:10.5281/zenodo.23741408. See also `CITATION.cff`.

## Licence

Arrays: CC BY 4.0. Code: MIT. See `LICENSE`.
