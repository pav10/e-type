#!/usr/bin/env python3
"""Find and certify a short vector in an external Gram matrix (coordinates in ITS basis).
usage: srfp_kill.py FILE target maxblock OUT"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.core import *
G = pari(open(sys.argv[1]).read().strip()); target = int(sys.argv[2]); mb = int(sys.argv[3])
t = time.time()
nr, x = bkz_shortest(G, block=mb, loops=8, target=target, verbose=True, input_coords=True)
v = pari.matrix(len(G), 1, x); exact = (pari.mattranspose(v) * G * v)[0, 0]
print("vector found: norm %s (recomputed exactly in the input Gram matrix), %.0fs" % (exact, time.time() - t))
if exact <= target:
    open(sys.argv[4], "w").write(" ".join(map(str, x)) + "\n")
    print("saved to", sys.argv[4])
