#!/usr/bin/env python3
"""Validation harness -- run this first.  Every check is against known data.

 1. theta module: kissing numbers of E8, Leech, P48, Gamma72, D4, BW16, K12, ...
 2. tensor_OK: O(x)O, T(x)T over Z[(1+sqrt(-11))/2]
 3. ideal lattices: Leech with a Z[eta]-structure over Q(zeta_13, sqrt(-11))
 4. positive/negative control (Coulangeon-Nebe, arXiv:1201.1832, Cor. 4.2):
    P (x) T has rank 48 and minimum 4 or 6; the kill test must find norm 4 for
    some structures P and must NOT find anything below 6 for the others (P48n).
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
from etype.theta import target, extremal_theta

known = [(8, 1, 240), (24, 1, 196560), (48, 1, 52416000), (72, 1, 6218175600), (4, 2, 24), (16, 2, 4320),
         (2, 3, 6), (12, 3, 756), (14, 3, 756), (6, 7, 42), (8, 5, 120)]
for n, l, k in known:
    t = target(n, l)
    assert t["kissing"] == k, (n, l, t)
assert extremal_theta(12, 3, 11)[0] == [1, 0, 756, 4032, 20412, 60480, 139860, 326592, 652428, 1020096, 2000376]
print("theta: ok")

K = (3, 1)
T, JT = herm_to_z([[(2, 0), (0, 1)], [(1, -1), (2, 0)]], K)       # [[2, eta], [conj(eta), 2]]
assert pari.matdet(T) == 121 and minimum_exact(T) == 4 and is_structure(T, JT, K)
O, JO = standard(1, K)
G, J = tensor_OK(O, JO, O, JO, K)
assert pari.matdet(G) == 11 and is_structure(G, J, K)
G, J = tensor_OK(T, JT, T, JT, K)
assert pari.matdet(G) == 11**4 and minimum_exact(G) == 4          # identity tensor: h = 2
print("tensor: ok")

L = leech_structures(K, ms=(13,))
assert len(L) > 0
res = set()
for P in L[:16]:
    H, JH = tensor_OK(P["G"], P["J"], T, JT, K)
    assert len(H) == 48 and pari.matdet(H) == 1 and is_even(H) and is_structure(H, JH, K)
    res.add(bkz_shortest(H, block=40, target=4)[0])
print("rank-48 control: shortest vectors found", sorted(res))
assert res == {4, 6}, "expected both outcomes (min 4 and the extremal P48n with min 6)"
print("ALL OK")
