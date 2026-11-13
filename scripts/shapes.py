#!/usr/bin/env python3
"""Enumerate the admissible Gram shapes S(W, T) for the d=19 brick (6,20).
usage: shapes.py T part nparts      (T = Hermitian norm bound: 4 <-> Z-norm 8)"""
import sys, os, time, pickle
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
from etype.certify import *
K = (5, 1)
T = float(sys.argv[1]); part = (int(sys.argv[2]), int(sys.argv[3]))
A, JA = brick_walk(K, 4, steps=2500)[(6, 20)]
vecs, ex, W = ok_basis(A, JA, K)
t = time.time()
R = enumerate_X(W, K, T, part=part, verbose=True)
print("%d numeric candidates (%.0fs)" % (len(R), time.time() - t), flush=True)
good = []
for Tv, dg, off in R:
    G, J = X_lattice(dg, off, K)
    if count_short(G, 2) == 0:
        good.append((round(Tv, 6), dg, off, int(pari.matdet(G))))
print("%d admissible shapes" % len(good), flush=True)
pickle.dump((vecs, ex, good), open(os.path.join(RES, "X_d19_T%g_part%d_%d.pkl" % (T, part[0], part[1])), "wb"))
