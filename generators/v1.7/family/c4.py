"""C4: two rich starters + one binary starter, Z_{k-1} on columns 0..k-2, column k-1 fixed (f = 1).

Statement: for k in {10,...,16} U {21, 22, 23} there are u, v in PG(1,5)^k and b in {0,1}^k such that
  {6 constant rows} U {g(rho^r w) : g in PGL(2,5), r in Z_{k-1}, w in {u, v, b}}
(rho rotates columns 0..k-2 and fixes column k-1) is a CA(4,k,6) with at most 270(k-1)+6 rows.
Values: 10:2436 11:2706 12:2976 13:3246 14:3516 15:3786 16:4056 21:5406 22:5676 23:5946;
live bars: 2616 3000 3216 3606 3846 4101 4111 5487 5678 6023 (lab: 13:3510, 21:5477, 22:5671).
This is the Meagher-Stevens / Maity et al. section 3.1 fixed-point variant, at g = 6 with a binary starter.
Existential: a miss within budget decides nothing.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import family_check  # noqa: E402

STATEMENT = ("For k in {10..16, 21, 22, 23}: exist u,v in PG(1,5)^k, b in {0,1}^k with constants U (Z_{k-1} x PGL(2,5))-"
             "development of {u,v,b} (one fixed column) a CA(4,k,6) with <= 270(k-1)+6 rows")

if __name__ == "__main__":
    family_check("C4", STATEMENT, [10, 11, 12, 13, 14, 15, 16, 21, 22, 23], f=1, n_rich=2, n_bin=1,
                 total_seconds=float(os.environ.get("C4_SECONDS", "4200")))
