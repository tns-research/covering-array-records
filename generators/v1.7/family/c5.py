"""C5: two rich starters + one binary starter, Z_{k-2} on columns 0..k-3, columns k-2, k-1 fixed (f = 2).

Statement: for k in {13,...,17} U {21,...,25} there are u, v in PG(1,5)^k and b in {0,1}^k such that
  {6 constant rows} U {g(rho^r w) : g in PGL(2,5), r in Z_{k-2}, w in {u, v, b}}
is a CA(4,k,6) with at most 270(k-2)+6 rows.
Values: 13:2976 14:3246 15:3516 16:3786 17:4056 21:5136 22:5406 23:5676 24:5946 25:6216;
live bars: 3606 3846 4101 4111 4148 5487 5678 6023 6147 6302 (lab: 13:3510, 21:5477, 22:5671).
Deleting a fixed column gives C4 at k-1 and C3 at k-2, so C5 is the strongest of the three.
(f = 3 is impossible with two rich starters: the 4-sets made of the 3 fixed columns and one cyclic column
need fixed triples of shapes abc, aab, aba, abb, one per rich starter.)
Existential: a miss within budget decides nothing.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import family_check  # noqa: E402

STATEMENT = ("For k in {13..17, 21..25}: exist u,v in PG(1,5)^k, b in {0,1}^k with constants U (Z_{k-2} x PGL(2,5))-"
             "development of {u,v,b} (two fixed columns) a CA(4,k,6) with <= 270(k-2)+6 rows")

if __name__ == "__main__":
    family_check("C5", STATEMENT, [13, 14, 15, 16, 17, 21, 22, 23, 24, 25], f=2, n_rich=2, n_bin=1,
                 total_seconds=float(os.environ.get("C5_SECONDS", "4200")))
