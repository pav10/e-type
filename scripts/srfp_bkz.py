#!/usr/bin/env python3
"""Deeper BKZ search for short vectors in an external Gram matrix:  srfp_bkz.py FILE [block]"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.core import *
G = lll(pari(open(sys.argv[1]).read().strip())); blk = int(sys.argv[2]) if len(sys.argv) > 2 else 70
t = time.time(); nr, x = bkz_shortest(G, block=blk, loops=8, verbose=True)
v = pari.matrix(len(G), 1, x)
print("shortest found %d (exact %s) in %.0fs" % (nr, (pari.mattranspose(v) * G * v)[0, 0], time.time() - t))
