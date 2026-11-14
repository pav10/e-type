#!/usr/bin/env python3
"""Recover the Gram shape X of a verified kill vector z = sum b_i (x) x_i, then test the same
shape against other survivors with the same brick; every hit is turned into an explicit vector
whose norm is computed exactly.
usage: transfer_kill.py N t KILLED_SURVIVOR.json OTHER_SURVIVOR.json [...]"""
import sys, os, re, json, random, pickle
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
from etype.search import *
from etype.certify import ok_basis, herm_gram, X_lattice
from etype.represent import herm_aut, MinHost, code_freq, embeds_fast
from etype.exact import explicit_vector
from etype.core import _kron
N, t = int(sys.argv[1]), int(sys.argv[2]); K = (N, t); D = disc(K)
def parse(f):
    m = re.match(r"SURVIVOR_rank96_d(\d+)_\d+_\d+_\w+?_(\d+)_(\d+)\.json", os.path.basename(f))
    return int(m.group(2)), int(m.group(3))
b0, l0 = parse(sys.argv[3])
bricks = json.load(open(os.path.join(RES, "bricks_d%d_m4.json" % D)))
A, JA = pari(bricks[b0]["G"]), pari(bricks[b0]["J"])
vecs = ok_basis(A, JA, K)[0]
leech = load("leech_d%d.json" % D); random.seed(0); random.shuffle(leech)
kill = pickle.load(open(os.path.join(RES, "kill_" + os.path.basename(sys.argv[3])[:-5] + ".pkl"), "rb"))
P0, JP0 = leech[l0]["G"], leech[l0]["J"]
n1, n2 = 8, 24
def Pm_of(JP):
    return -(_kron(JA, pari.matid(n2)) + _kron(pari.matid(n1), JP) - t * pari.matid(n1 * n2)) * _kron(2 * JA - t * pari.matid(n1), pari.matid(n2))
Pm = Pm_of(JP0); B = pari.mathnf(Pm)
z = B * pari.matrix(96, 1, kill["z"]) / D
Phi = pari.matconcat([Pm * pari.matrix(192, 1, [vecs[k][a] * (1 if b == l else 0) for a in range(8) for b in range(24)]) / D
                      for k in range(4) for l in range(24)])
x = pari.matsolve(pari.mattranspose(Phi) * Phi, pari.mattranspose(Phi) * z)
assert Phi * x == z and pari.denominator(x) == 1
xs = [[int(x[24 * k + l, 0]) for l in range(24)] for k in range(4)]
exX, X = herm_gram(P0, JP0, K, xs)
print("kill vector of norm %d recovered; tensor rank %d; tr(conj W X) = %.4f" %
      (kill["norm"], np.linalg.matrix_rank(X, tol=1e-8), np.trace(np.conj(ok_basis(A, JA, K)[2]) @ X).real), flush=True)
# shape lattice M = (O_K^4, X)
dg = tuple(int(exX[k][k][0]) for k in range(4)); off = {}
for p in range(4):
    for q in range(p + 1, 4):
        u, v = exX[p][q]                               # h = u + v eta;  a + b eta = h * (2 eta - t)
        a = -t * u - 2 * N * v; b = 2 * u + t * v
        off[(p, q)] = (int(a), int(b))
GM, JM = X_lattice(dg, off, K)
assert pari.mattranspose(pari.matrix(1,1,[0])) is not None
for f in sys.argv[4:]:
    b, l = parse(f)
    assert b == b0, "different brick"
    P, JP = leech[l]["G"], leech[l]["J"]
    o, gens = herm_aut(P, JP, K); host = MinHost(P, JP, K, gens=gens); fq = code_freq(host)
    ok, sol = embeds_fast(host, GM, JM, K, fq)
    if not ok:
        print("%s: the killing shape does not embed (survivor not decided by this test)" % os.path.basename(f), flush=True)
        continue
    zc, nz = explicit_vector(K, A, JA, P, JP, vecs, host, fq, (4, GM, JM))
    print("%s: KILLED, explicit vector of norm %d" % (os.path.basename(f), nz), flush=True)
    pickle.dump({"norm": nz, "z": zc}, open(os.path.join(RES, "kill_" + os.path.basename(f)[:-5] + ".pkl"), "wb"))
