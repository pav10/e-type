#!/usr/bin/env python3
"""End-to-end positive control on d = 19 with the production brick (6,20):
for killed candidates (6,20) (x) Leech#i, decompose a short vector z = sum b_i (x) x_i
(b = the fixed O_K-basis of A used by the enumeration), compute X = Gram(x_i) and its
rank r.  For r = 4:  (b) X must occur EXACTLY in the production list of rank-4
shapes, and (c) the embedding search must find its class in P."""
import sys, os, random, time, pickle
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
from etype.certify import *
from etype.represent import *
from etype.core import _kron
K = (5, 1); N, t = K; d = 19
A, JA = brick_walk(K, 4, steps=2500)[(6, 20)]
vecs, ex, W = ok_basis(A, JA, K)
leech = load("leech_d19.json"); random.seed(0); random.shuffle(leech)
shapes = []
for part in range(3):
    shapes += pickle.load(open(os.path.join(RES, "X_d19_T4_part%d_3.pkl" % part), "rb"))[2]
v0 = pickle.load(open(os.path.join(RES, "X_d19_T4_part1_3.pkl"), "rb"))[0]
assert v0 == vecs, "basis of A differs from the one used by the enumeration"
index = {}
for Tv, dg, off, _det in shapes:
    index[(tuple(dg), tuple(sorted(off.items())))] = (Tv, dg, off)
print(len(index), "rank-4 shapes indexed", flush=True)
n1, n2 = 8, 24
Pm0 = lambda JP: -(_kron(JA, pari.matid(n2)) + _kron(pari.matid(n1), JP) - t * pari.matid(n1 * n2)) * _kron(2 * JA - t * pari.matid(n1), pari.matid(n2))
nfound = 0
for i in range(len(leech)):
    P, JP = leech[i]["G"], leech[i]["J"]
    Pm = Pm0(JP); B = pari.mathnf(Pm)
    G = pari.mattranspose(B) * _kron(A, P) * B / d**2
    U = pari.qflllgram(G, 1); Gr = pari.mattranspose(U) * G * U
    nr, c = bkz_shortest(Gr, block=30, target=8)
    if nr > 8:
        continue
    z = B * U * pari.matrix(96, 1, c) / d
    Phi = pari.matconcat([Pm * pari.matrix(192, 1, [vecs[k][a] * (1 if b == l else 0) for a in range(8) for b in range(24)]) / d
                          for k in range(4) for l in range(24)])
    x = pari.matsolve(pari.mattranspose(Phi) * Phi, pari.mattranspose(Phi) * z)
    assert Phi * x == z and pari.denominator(x) == 1
    xs = [[int(x[24 * k + l, 0]) for l in range(24)] for k in range(4)]
    exX, X = herm_gram(P, JP, K, xs)
    r = int(np.linalg.matrix_rank(X, tol=1e-8))
    print("Leech#%d: norm %d, tensor rank %d, tr(conj W X) = %.4f" % (i, nr, r, np.trace(np.conj(W) @ X).real), flush=True)
    if r < 4:
        continue
    # exact lookup: diag ints, off-diagonal (a,b) with h = (a + b eta)/sqrt(-19)
    dg = tuple(int(exX[k][k][0]) for k in range(4)); off = {}
    for p in range(4):
        for q in range(p + 1, 4):
            u, v = exX[p][q]                                  # h = u + v eta ;  a + b eta = h * (2 eta - 1)
            a = -u - 2 * N * v; b = 2 * u + v * (t + 0)        # (u + v eta)(2eta - 1) = -u + (2u - v) eta + 2v eta^2
            a = -u - 2 * N * v; b = 2 * u - v + 2 * v * t
            off[(p, q)] = (int(a), int(b))
    hit = index.get((dg, tuple(sorted(off.items()))))
    print("   (b) X in the production list:", hit is not None, flush=True)
    assert hit is not None, "CONTROL FAILED (b)"
    Gx, Jx = X_lattice(dg, off, K)
    o, gens = herm_aut(P, JP, K); host = MinHost(P, JP, K, gens=gens); fq = code_freq(host)
    ok, sol = embeds_fast(host, Gx, Jx, K, fq)
    print("   (c) embedding found:", ok, flush=True)
    assert ok, "CONTROL FAILED (c)"
    nfound += 1
    if nfound >= 2:
        break
print("CONTROL PASSED" if nfound else "no rank-4 killer found")
