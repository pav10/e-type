#!/usr/bin/env python3
"""Stronger kill test on a saved candidate:  attack.py results/FILE.json [maxblock]"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
it = load(os.path.basename(sys.argv[1]))[0]
G = it["G"]; mb = int(sys.argv[2]) if len(sys.argv) > 2 else 70
print("rank", len(G), "det", pari.matdet(G), "even", is_even(G), flush=True)
t = time.time()
nr, x = bkz_shortest(G, block=mb, loops=6, verbose=True)
print("shortest found", nr, "in %.0fs" % (time.time() - t)); print(x)
