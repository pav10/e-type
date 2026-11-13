#!/usr/bin/env python3
"""Exact decision for a rank-96 candidate  A (x)_{O_K} P  saved by rank96_tensor.py.

usage: exact96.py N t SURVIVOR_FILE.json [nproc]
The survivor file name encodes the brick key; the brick is re-read from
results/bricks_d<d>_m4.json and P from results/leech_d<d>.json (same shuffle)."""
import sys, os, re, json, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
from etype.exact import decide
N, t = int(sys.argv[1]), int(sys.argv[2]); K = (N, t); d = disc(K)
f = os.path.basename(sys.argv[3]); nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 3
m = re.match(r"SURVIVOR_rank96_d(\d+)_(\d+)_(\d+)_(\w+)_(\d+)_(\d+)\.json", f)
assert m and int(m.group(1)) == d, f
bidx, lidx = int(m.group(5)), int(m.group(6))
b = json.load(open(os.path.join(RES, "bricks_d%d_m4.json" % d)))[bidx]
A, JA = pari(b["G"]), pari(b["J"])
leech = load("leech_d%d.json" % d); random.seed(0); random.shuffle(leech)
P, JP = leech[lidx]["G"], leech[lidx]["J"]
G, J = tensor_OK(A, JA, P, JP, K)
print("candidate: brick #%d %s (x) Leech#%d, rank %d, det %s" % (bidx, (b["Zmin"], b["kissing"], b["parity"]), lidx, len(G), pari.matdet(G)), flush=True)
out = os.path.join(RES, "exact_" + f[:-5])
cache = os.path.join(RES, "brick_shapes_d%d_b%d.pkl" % (d, bidx))
res = decide(K, A, JA, P, JP, out, nproc=nproc, log=lambda s: print(s, flush=True), cache=cache)
print("RESULT:", res[0], res[1])
