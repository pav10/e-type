#!/usr/bin/env python3
"""Survey the genera of Hermitian unimodular O_K-lattices of rank 2, 3, 4 by
random Kneser walks; prints the largest Hermitian minimum found per field.
These are the possible small 'bricks' A in  A (x)_{O_K} B."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
fields = [(1, 1), (2, 1), (3, 1), (2, 0), (4, 1), (5, 1), (5, 0), (6, 1), (6, 0), (8, 1), (9, 1), (10, 1)]
for K in fields:
    for m in (2, 3, 4):
        s = brick_walk(K, m, steps=400 if m < 4 else 1500)
        print("d=%3d rank %d: max Hermitian min %d   (%d invariant classes seen)" % (disc(K), m, max(s)[0] // 2, len(s)), flush=True)
