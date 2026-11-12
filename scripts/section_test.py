#!/usr/bin/env python3
"""Rank-2 sections test: for every rank-2 section A2 of the brick A spanned by
two short vectors, build the 48-dimensional sublattice  A2 (x) P  of  A (x) P
and search it hard (BKZ up to block 48 ~ near-exact SVP in dimension 48).
Covers ALL tensors a1(x)p1 + a2(x)p2 with a1, a2 in A2 and p1, p2 arbitrary.

usage: section_test.py LEECH_INDEX [normA]"""
import sys, os, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
from etype.search import _kconj, _kmul
K = (5, 1)
idx = int(sys.argv[1]); normA = int(sys.argv[2]) if len(sys.argv) > 2 else 8
A, JA = brick_walk(K, 4, steps=2500)[(6, 20)]
leech = load("leech_d19.json"); random.seed(0); random.shuffle(leech)
P, JP = leech[idx]["G"], leech[idx]["J"]
types = pair_types(A, JA, K, short_vectors(A, normA))
best = None
for (h1, h2, u, v) in sorted(types, key=lambda x: (x[0] * x[1] - _kmul((x[2], x[3]), _kconj((x[2], x[3]), K), K)[0])):
    dA = h1 * h2 - _kmul((u, v), _kconj((u, v), K), K)[0]
    if dA <= 0:
        continue
    A2, J2 = herm_to_z([[(h1, 0), (u, v)], [_kconj((u, v), K), (h2, 0)]], K)
    t = time.time()
    G, J = tensor_OK(A2, J2, P, JP, K)
    nr, x = bkz_shortest(G, block=48, loops=4)
    best = nr if best is None else min(best, nr)
    print("section (%d,%d,%s+%s eta) det %s: shortest found %d   (%.0fs)" % (h1, h2, u, v, dA, nr, time.time() - t), flush=True)
print("Leech#%d: smallest norm found on rank-2 sections: %d" % (idx, best))
