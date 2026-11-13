#!/usr/bin/env python3
"""Turn every recorded embedding into an explicit vector z of A (x) P and compute its
norm exactly in the 96-dimensional Gram matrix.   usage: verify_kill.py idx"""
import sys, os, random, pickle
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
from etype.certify import *
from etype.represent import *
from etype.core import _kron
K = (5, 1); N, t = K; d = 19
idx = int(sys.argv[1]); OUT = os.path.join(RES, "decide_d19_%d" % idx)
A, JA = brick_walk(K, 4, steps=2500)[(6, 20)]
vecs, ex, W = ok_basis(A, JA, K)
leech = load("leech_d19.json"); random.seed(0); random.shuffle(leech)
P, JP = leech[idx]["G"], leech[idx]["J"]
o, gens = herm_aut(P, JP, K); host = MinHost(P, JP, K, gens=gens); fq = code_freq(host)
cl = pickle.load(open(os.path.join(OUT, "classes.pkl"), "rb"))
hits = sorted({int(l.split()[0]) for f in os.listdir(OUT) if f.startswith("embed_log")
               for l in open(os.path.join(OUT, f)) if l.split()[2] == "True"})
Pm = -(_kron(JA, pari.matid(24)) + _kron(pari.matid(8), JP) - t * pari.matid(192)) * _kron(2 * JA - t * pari.matid(8), pari.matid(24))
GG = _kron(A, P)
# the 96-dim lattice itself, to check z is a lattice vector
B = pari.mathnf(Pm); G96 = pari.mattranspose(B) * GG * B / d**2
assert pari.matdet(G96) == 1
out = []
for ci in hits:
    r, Gs, Js = cl[ci]; Gs, Js = pari(Gs), pari(Js)
    F = embeds_fast(host, Gs, Js, K, fq, prepare_only="map")
    assert F is not None
    Fp = pari.matrix(24, 8, [int(x) for x in F.ravel()])
    assert pari.mattranspose(Fp) * P * Fp == Gs and Fp * Js == JP * Fp      # isometric, O_K-linear
    xs = [[int(F[l, 2 * i]) for l in range(24)] for i in range(4)]
    z = sum((Pm * pari.matrix(192, 1, [a * b for a in vecs[i] for b in xs[i]]) for i in range(4)), pari.matrix(192, 1)) / d
    c = pari.matsolve(B, z * d)                                              # coordinates in the lattice basis
    assert pari.denominator(c) == 1, "z not in the lattice"
    nz = (pari.mattranspose(c) * G96 * c)[0, 0]
    print("class %d: explicit z in A (x) P, coordinates integral, Z-norm = %s" % (ci, nz), flush=True)
    out.append((ci, int(nz), [int(c[k, 0]) for k in range(96)]))
    if len(out) >= 3:
        break
pickle.dump(out, open(os.path.join(OUT, "kill_vectors.pkl"), "wb"))
print("minimum norm among verified vectors:", min(o[1] for o in out))
