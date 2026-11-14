#!/usr/bin/env python3
"""Fast exact KILL attempt (not a proof of extremality): embedding search over the rank-4
shapes available so far (finished parts of the brick-shape run), cheapest classes first;
the first embedding is turned into an explicit vector whose norm is computed exactly.
usage: kill_stream.py N t SURVIVOR_FILE.json [nproc] [maxcost]"""
import sys, os, re, json, random, glob, pickle, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
from etype.exact import classes, explicit_vector, _embed_worker, _G
from etype.represent import herm_aut, MinHost, code_freq, embeds_fast
from etype.certify import ok_basis
import multiprocessing as mp
import numpy as np
N, t = int(sys.argv[1]), int(sys.argv[2]); K = (N, t); d = disc(K)
f = os.path.basename(sys.argv[3]); nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 2
maxcost = float(sys.argv[5]) if len(sys.argv) > 5 else 1e5
m = re.match(r"SURVIVOR_rank96_d(\d+)_\d+_\d+_\w+?_(\d+)_(\d+)\.json", f)
bidx, lidx = int(m.group(2)), int(m.group(3))
b = json.load(open(os.path.join(RES, "bricks_d%d_m4.json" % d)))[bidx]
A, JA = pari(b["G"]), pari(b["J"])
vecsA = ok_basis(A, JA, K)[0]
leech = load("leech_d%d.json" % d); random.seed(0); random.shuffle(leech)
P, JP = leech[lidx]["G"], leech[lidx]["J"]
parts = sorted(glob.glob(os.path.join(RES, "brick_shapes_d%d_b%d.pkl.parts" % (d, bidx), "part_*.pkl")))
shapes = [(4,) + s for p in parts for s in pickle.load(open(p, "rb"))]
print("%s: %d parts, %d rank-4 shapes" % (f, len(parts), len(shapes)), flush=True)
o, gens = herm_aut(P, JP, K)
host = MinHost(P, JP, K, gens=gens); fq = code_freq(host)
allowed = set()
for i in host.reps:
    allowed |= set(np.unique(host.code(np.array([i]), np.arange(host.V.shape[1]))[0]).tolist())
print("|Aut_h(P)| = %d, %d orbits" % (o, len(host.reps)), flush=True)
cl = classes(shapes, K, allowed, log=lambda s: print(s, flush=True))
costs = []
for r, G, J in cl:
    res = embeds_fast(host, G, J, K, fq, prepare_only="cost")
    costs.append(None if isinstance(res[1], str) else res[0])
order = sorted((i for i in range(len(cl)) if costs[i] is not None and costs[i] <= maxcost), key=lambda i: costs[i])
print("searching %d classes" % len(order), flush=True)
_G.update(host=host, fq=fq, K=K, cl=cl)
hit = None; t0 = time.time(); n = 0
with mp.get_context("fork").Pool(nproc) as pool:
    for ci, ok in pool.imap_unordered(_embed_worker, order, chunksize=4):
        n += 1
        if ok:
            hit = ci; pool.terminate(); break
        if n % 500 == 0:
            print("   %d / %d, no embedding (%.0fs)" % (n, len(order), time.time() - t0), flush=True)
if hit is None:
    print("NO KILL among the searched classes (this is NOT a proof)"); sys.exit()
z, nz = explicit_vector(K, A, JA, P, JP, vecsA, host, fq, cl[hit])
print("KILLED: explicit vector of norm %d, integral coordinates in the lattice basis" % nz)
pickle.dump({"norm": nz, "z": z}, open(os.path.join(RES, "kill_" + f[:-5] + ".pkl"), "wb"))
