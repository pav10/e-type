#!/usr/bin/env python3
"""Recompute the ideal-lattice Leech structures over Q(zeta_m, eta) together with
the action of zeta_m, and store those matching given indices of the shuffled
list leech_d19.json.   usage: leech_with_zeta.py m idx [idx ...]"""
import sys, os, random, pickle
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
K = (5, 1); N, t = K
m = int(sys.argv[1]); want = list(map(int, sys.argv[2:]))
leech = load("leech_d19.json"); random.seed(0); random.shuffle(leech)
pari("P = polredbest(polcompositum(polcyclo(%d), x^2 - %d*x + %d)[1])" % (m, t, N))
pari("P = subst(P, x, t); bnf = bnfinit(P, 1); et = Mod(nfroots(bnf.nf, x^2 - %d*x + %d)[1], P); ze = Mod(nfroots(bnf.nf, polcyclo(%d))[1], P);" % (t, N, m))
R = pari("ideal_lattices(bnf, 1, [et, ze], 4)")
found = {}
for i in range(len(R)):
    G = R[i][2]
    U = pari.qflllgram(G, 1)
    G2 = pari.mattranspose(U) * G * U
    for idx in want:
        if idx not in found and G2 == leech[idx]["G"]:
            J2 = U ** -1 * R[i][3][0] * U; Z2 = U ** -1 * R[i][3][1] * U
            assert J2 == leech[idx]["J"] and pari.mattranspose(Z2) * G2 * Z2 == G2 and Z2 * J2 == J2 * Z2
            found[idx] = {"G": str(G2), "J": str(J2), "Z": str(Z2), "m": m}
print("matched:", sorted(found))
for idx, v in found.items():
    pickle.dump(v, open(os.path.join(RES, "leech_d19_%d_zeta.pkl" % idx), "wb"))
