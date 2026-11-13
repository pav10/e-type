#!/usr/bin/env python3
"""Precompute sections + Gram shapes of a brick:  brick_shapes.py N t brick_index [nproc]"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
from etype.exact import brick_shapes
N, t, bi = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); K = (N, t); d = disc(K)
nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 1
b = json.load(open(os.path.join(RES, "bricks_d%d_m4.json" % d)))[bi]
r = brick_shapes(K, pari(b["G"]), pari(b["J"]), os.path.join(RES, "brick_shapes_d%d_b%d.pkl" % (d, bi)), nproc=nproc,
                 log=lambda s: print(s, flush=True))
print("done: %d shapes, %d non-free sections" % (len(r[0]), len(r[1])))
