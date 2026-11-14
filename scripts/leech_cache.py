#!/usr/bin/env python3
"""Leech lattices with O_K-structure from ideal lattices over Q(zeta_m)K, phi(m) = 12, cached per m.
usage: leech_cache.py N t m [m ...]      -> results/leech_d<d>_m<m>.json, merged into leech_d<d>.json"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
N, t = int(sys.argv[1]), int(sys.argv[2]); K = (N, t); d = disc(K)
for m in map(int, sys.argv[3:]):
    f = os.path.join(RES, "leech_d%d_m%d.json" % (d, m))
    if os.path.exists(f):
        print("m=%d cached" % m); continue
    t0 = time.time()
    L = leech_structures(K, ms=(m,))
    save(os.path.basename(f), L)
    print("m=%d: %d Leech structures (%.0fs)" % (m, len(L), time.time() - t0), flush=True)
allL = []
for m in (13, 21, 28, 36):
    f = os.path.join(RES, "leech_d%d_m%d.json" % (d, m))
    if os.path.exists(f):
        allL += load(os.path.basename(f))
save("leech_d%d.json" % d, allL)
print("leech_d%d.json: %d structures" % (d, len(allL)))
