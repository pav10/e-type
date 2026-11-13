#!/usr/bin/env python3
"""Control for the exact method on Coulangeon-Nebe, arXiv:1201.1832, Prop. 4.1/Cor. 4.2:
P (x) T over Z[(1+sqrt-11)/2], P a Hermitian Leech structure, T = [[2, eta],[conj eta, 2]].
min(P (x) T) = 4 (Z-norm) iff P contains one of the admissible rank-2 shapes, else 6.
 - for structures where BKZ finds norm 4: some shape must embed (search stops at first hit);
 - for the SL2(13).2-structure (P48n): no shape may embed (exhaustive over Aut-orbits)."""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
from etype.certify import *
from etype.represent import *
K = (3, 1); et = eta_c(K)
W = np.array([[2, et], [np.conj(et), 2]])
T, JT = herm_to_z([[(2, 0), (0, 1)], [(1, -1), (2, 0)]], K)
shapes = []
for Tv, dg, off in enumerate_X(W, K, 2.0):
    G, J = X_lattice(dg, off, K)
    if count_short(G, 2) == 0:
        shapes.append((G, J))
print(len(shapes), "admissible rank-2 shapes", flush=True)
L = leech_structures(K, ms=(13,), verbose=False)
done = {4: 0, 6: 0}
for P in L:
    H, JH = tensor_OK(P["G"], P["J"], T, JT, K)
    mn = bkz_shortest(H, block=40, target=4)[0]
    t0 = time.time()
    if mn == 4 and done[4] < 2:
        # the norm-4 vector z = b1 (x) x1 + b2 (x) x2 found by BKZ determines a shape X = Gram(x1, x2);
        # it must be one of the enumerated shapes, and the exhaustive search must find it in P.
        nr, zc = bkz_shortest(H, block=40, target=4)
        o, gens = herm_aut(P["G"], P["J"], K)
        host = MinHost(P["G"], P["J"], K, gens=gens)
        fq = code_freq(host)
        k, sol = embeds_any(host, shapes, K, fq)
        assert k is not None, "CONTROL FAILED: min 4 but no shape embeds"
        found = "shape %d (Z-det %s)" % (k, pari.matdet(shapes[k][0]))
        done[4] += 1
        print("min-4 structure (|Aut| = %d): shapes embedded %s (%.0fs)  OK" % (o, found, time.time() - t0), flush=True)
    elif mn == 6 and done[6] < 1:
        o, gens = herm_aut(P["G"], P["J"], K)
        host = MinHost(P["G"], P["J"], K, gens=gens)
        res = [embeds(host, G, J, K)[0] for G, J in shapes]
        assert not any(res), "CONTROL FAILED: extremal but a shape embeds"
        done[6] += 1
        print("min-6 structure (|Aut| = %d, %d orbits): no shape embeds (%.0fs)  OK" % (o, len(host.reps), time.time() - t0), flush=True)
    if done[4] >= 2 and done[6] >= 1:
        break
print("CONTROL PASSED")
