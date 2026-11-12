#!/usr/bin/env python3
"""Rank-96 candidates  A (x)_{O_K} Leech,  A Hermitian unimodular of rank 4.

Split vectors have norm  hmin(A) * min(Leech) = 4 * hmin(A), so hmin(A) >= 3 is
necessary (first possible for d = 19).  Every candidate is even unimodular of
rank 96; we look for a vector of norm <= 8 (kill).  A survivor of BKZ-60 would
be a serious candidate for an extremal lattice (min 10) and goes to certification.

usage:  rank96_tensor.py N t [block] [max number of Leech structures]      e.g.  rank96_tensor.py 5 1   (d = 19)
"""
import sys, time, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *

K = (int(sys.argv[1]), int(sys.argv[2]))
block = int(sys.argv[3]) if len(sys.argv) > 3 else 50
d = disc(K)
print("K = Q(sqrt(-%d))" % d, flush=True)
bricks = {k: v for k, v in brick_walk(K, 4, steps=2500).items() if k[0] >= 6}
print("rank-4 bricks with Hermitian min >= 3 (Zmin, kissing):", sorted(bricks), flush=True)
cache = "leech_d%d.json" % d
if os.path.exists(os.path.join(RES, cache)):
    leech = load(cache)
else:
    leech = leech_structures(K)
    save(cache, leech)
random.seed(0); random.shuffle(leech)
limit = int(sys.argv[4]) if len(sys.argv) > 4 else len(leech)
leech = leech[:limit]
print("%d Leech structures over O_K in use" % len(leech), flush=True)
rows = []
for (mn, kiss), (A, JA) in sorted(bricks.items()):
    for i, P in enumerate(leech):
        t0 = time.time()
        G, J = tensor_OK(A, JA, P["G"], P["J"], K)
        assert len(G) == 96 and pari.matdet(G) == 1 and is_even(G) and is_structure(G, J, K)
        dead, nr = kill(G, 8, block=block)
        rows.append({"d": d, "brick": [mn, kiss], "leech": i, "m": P["m"], "shortest_found": int(nr)})
        print("brick %s (x) Leech#%d (m=%d): shortest found %d  %s  (%.0fs)" %
              ((mn, kiss), i, P["m"], nr, "DEAD" if dead else "*** SURVIVOR ***", time.time() - t0), flush=True)
        if not dead:
            save("SURVIVOR_rank96_d%d_%d_%d.json" % (d, mn, i), [{"G": G, "J": J}])
save("rank96_tensor_d%d.json" % d, rows)
