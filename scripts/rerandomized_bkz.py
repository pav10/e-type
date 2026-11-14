#!/usr/bin/env python3
"""Heuristic kill test with rerandomised bases: repeat (random unimodular transform -> LLL -> BKZ)
and keep the shortest vector.  usage: rerandomized_bkz.py SURVIVOR.json [tries] [block] [target]"""
import sys, os, random, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
G = load(os.path.basename(sys.argv[1]))[0]["G"]
tries = int(sys.argv[2]) if len(sys.argv) > 2 else 8
block = int(sys.argv[3]) if len(sys.argv) > 3 else 40
target = int(sys.argv[4]) if len(sys.argv) > 4 else 8
n = len(G); best = None
for k in range(tries):
    random.seed(k)
    U = pari.matid(n)
    for _ in range(4 * n):                                   # random elementary operations
        i, j = random.sample(range(n), 2)
        U[i, j] = U[i, j] + random.choice([-1, 1])
    Gr = pari.mattranspose(U) * G * U
    t0 = time.time()
    nr, x = bkz_shortest(Gr, block=block, loops=4)
    best = nr if best is None else min(best, nr)
    print("try %d: shortest found %d (%.0fs)" % (k, nr, time.time() - t0), flush=True)
    if best <= target:
        break
print("BEST", best)
