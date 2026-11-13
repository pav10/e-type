#!/usr/bin/env python3
"""Rank-m Hermitian unimodular O_K-lattices ("bricks") with large Hermitian minimum.

Starts: O_K^m (standard form) and ideal lattices over Q(zeta_k)K with [Q(zeta_k)K : K] = m,
alpha I conj(I) Diff_L = Diff_K  (Hermitian unimodular; reaches the EVEN genera that
O_K^m cannot reach by neighbours).  Random Kneser walks at split primes from each start,
classes up to O_K-linear isometry (qfisom).  Cross-check against Schiemann,
J. Symb. Comp. 26 (1998), Tables 1-2.

usage: bricks.py N t [m] [steps]      (eta^2 - t eta + N = 0)
"""
import sys, os, time, random, json, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
from etype.represent import herm_isom

N, t = int(sys.argv[1]), int(sys.argv[2]); K = (N, t); d = disc(K)
m = int(sys.argv[3]) if len(sys.argv) > 3 else 4
steps = int(sys.argv[4]) if len(sys.argv) > 4 else 4000

def parity(G):
    """'even' if every h(x,x) is even (all Z-norms = 2h(x,x) divisible by 4)."""
    n = len(G)
    ok = all(int(G[i, i]) % 4 == 0 for i in range(n)) and all(int(G[i, j]) % 2 == 0 for i in range(n) for j in range(n))
    return "even" if ok else "odd"

starts = [("standard",) + tuple(standard(m, K))]
cdiff = "2*et-%d" % t if t else "2*et"                          # generator of Diff_K (sqrt(-d) or 2 sqrt(-d))
for k in range(3, 80):
    if int(pari("eulerphi(%d)" % k)) == m and k % 4 != 2:
        try:
            fam = ideal_family(K, k, c=cdiff, maxsol=2)
        except Exception as e:
            print("  zeta_%d: %s" % (k, str(e)[:60])); continue
        fam = [(G, J) for G, J in fam if pari.matdet(G) == d ** m]
        for G, J in fam[:6]:
            starts.append(("ideal zeta_%d" % k, G, J))
print("d = %d, rank %d: %d starting lattices (%s)" % (d, m, len(starts), sorted(set(parity(s[1]) for s in starts))), flush=True)

classes = collections.defaultdict(list)                       # (parity, Zmin, kissing) -> [(G, J)]
def add(G, J):
    mn = minimum_exact(G); key = (parity(G), mn, count_short(G, mn))
    for G2, J2 in classes[key]:
        if herm_isom(G, J, G2, J2, K):
            return
    classes[key].append((G, J))

prs = split_primes(K, 60)[:4]
random.seed(1)
t0 = time.time()
for name, G, J in starts:
    assert is_structure(G, J, K) and pari.matdet(G) == d ** m
    add(G, J)
    for _ in range(steps // len(starts)):
        r = neighbour(G, J, K, random_vec(2 * m), random.choice(prs))
        if r is None:
            continue
        G, J = r
        add(G, J)
tot = collections.Counter()
for (par, mn, kiss), lst in classes.items():
    tot[par] += len(lst)
print("classes found: %s   (%.0fs)" % (dict(tot), time.time() - t0))
for key in sorted(classes, key=lambda k: (k[0], -k[1], k[2])):
    if key[1] >= 6:
        print("  %s  Hermitian min %d  kissing %d : %d class(es)" % (key[0], key[1] // 2, key[2], len(classes[key])))
best = [{"parity": k[0], "Zmin": k[1], "kissing": k[2], "G": str(G), "J": str(J)}
        for k, lst in classes.items() if k[1] >= 6 for G, J in lst]
json.dump(best, open(os.path.join(RES, "bricks_d%d_m%d.json" % (d, m)), "w"))
