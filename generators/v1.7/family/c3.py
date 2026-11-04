"""C3: two rich cyclic starters + one binary cyclic starter, no fixed column (f = 0).

Statement: for k in {9, 11, 14, 15} there are u, v in PG(1,5)^k and b in {0,1}^k such that
  {6 constant rows} U {g(rho^r w) : g in PGL(2,5), r in Z_k, w in {u, v, b}}
is a CA(4,k,6); it has at most 270k+6 rows (2436, 2976, 3786, 4056), below the live bars 2518, 3000, 3846, 4101.
(At k = 13 such an array, 3516 rows, was found by run 20261102-192344; this extends it.)
Check: annealing over (u, v, b); every hit is expanded and verified by the exact CA check.
Existential: a miss within budget decides nothing, so this script prints "survives" with the undecided k listed.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import family_check  # noqa: E402

STATEMENT = ("For k in {9,11,14,15}: exist u,v in PG(1,5)^k, b in {0,1}^k with constants U (Z_k x PGL(2,5))-development "
             "of {u,v,b} a CA(4,k,6) with <= 270k+6 rows")

if __name__ == "__main__":
    family_check("C3", STATEMENT, [9, 11, 14, 15], f=0, n_rich=2, n_bin=1,
                 total_seconds=float(os.environ.get("C3_SECONDS", "2000")))
