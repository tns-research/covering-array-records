# Arrays from the constructions of the thesis of K. Shokri

Explicit arrays built from published constructions, submitted to the ca-tools table on behalf of their authors.
They are not findings of this repository and are not part of its releases.

| Array | Construction | Command |
|---|---|---|
| CA(1016; 3, 66, 8) | C. J. Colbourn, C. Ingalls, J. Jedwab, M. Saaltink, K. W. Smith, B. Stevens, Combinatorial Theory 4 (2024), doi:10.5070/C64163845, restated as Theorem 3.3.7 of the thesis | `python3 build.py 3.3.7 8` |
| CA(8176; 3, 258, 16) | same | `python3 build.py 3.3.7 16` |
| CA(769; 3, 775, 5) | K. Shokri, PhD thesis, University of Ottawa, 2026, doi:10.20381/ruor-31663, Theorem 5.2.11 | `python3 build.py 5.2.11 5` |
| CA(364; 3, 336, 4) | same, Corollary 5.2.12 | `python3 build.py 5.2.12 4` |
| CA(765; 3, 651, 5) | same, Theorem 5.2.13 | `python3 build.py 5.2.13 5` |

Chapter 5 of the thesis states that its results are published in K. Shokri and L. Moura, J. Combin. Des. 33 (2025) 156-171.
Deterministic, no search, no seed; needs Python 3 and numpy. Each command writes `out/CA_N_3_k_q.txt`; check it with
`ca_verify FILE 3 q`. Code written with Claude Opus; arrays checked by the exhaustive verifier of the ca-tools table.
