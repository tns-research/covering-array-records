#!/usr/bin/env python3
"""Rebuild the five strength-3 covering arrays taken from the thesis of K. Shokri (2026, doi:10.20381/ruor-31663).

usage: python3 build.py THEOREM q [OUTDIR]
  THEOREM q  -> array
  3.3.7   8  -> CA(1016; 3, 66, 8)     construction of Colbourn, Ingalls, Jedwab, Saaltink, Smith, Stevens (2024),
  3.3.7  16  -> CA(8176; 3, 258, 16)   restated as Theorem 3.3.7 of the thesis
  5.2.11  5  -> CA(769; 3, 775, 5)     Theorem 5.2.11 of the thesis
  5.2.12  4  -> CA(364; 3, 336, 4)     Corollary 5.2.12 of the thesis
  5.2.13  5  -> CA(765; 3, 651, 5)     Theorem 5.2.13 of the thesis
Deterministic: no search, no seed. The choices the thesis makes are used as stated (default column order of the last
block, scaling step as stated, the primitive cubic returned by lib.primitive_cubic). Needs numpy. The array is written
to OUTDIR (default ./out) as CA_N_3_k_q.txt; check it with ca_verify (see ../../verify.sh). --check also runs the
built-in exhaustive triple check (slow for q = 16).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

args = [a for a in sys.argv[1:] if not a.startswith("--")]
if len(args) < 2:
    sys.exit(__doc__)
thm, q = args[0], int(args[1])
if len(args) > 2:
    lib.OUT = args[2]
if thm == "3.3.7":
    A, info = lib.build337(q)
    N, k = lib.formula("3.3.7", q)
else:
    A, info = lib.build5(thm, q)
    N, k = lib.formula(thm, q)
assert A.shape == (N, k), (A.shape, N, k)
if "--check" in sys.argv:
    f = lib.check3(A, q)
    assert f is None, f
print(lib.write_array(A, q))
