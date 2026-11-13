#!/usr/bin/env python3
"""Exact decision of  min(A (x)_{O_K} P) <= 8  for the d = 19 candidates.

  A : rank-4 Hermitian unimodular brick (Z-min 6, kissing 20) over O_K = Z[(1+sqrt-19)/2]
  P : Hermitian Leech structure leech_d19.json[idx]  (shuffled order of rank96_tensor.py)

z in A (x) P of rank r lies in A_r (x) P, A_r a primitive rank-r section of A (O_K is a
PID).  h(z,z) = tr(W_r conj X) with W_r = Gram(A_r), X = Gram of r independent vectors of P.
  r = 1 : h >= hmin(A) hmin(P) = 6 > 4.
  r >= 2: Coulangeon (Prop 3.2): 4 >= h >= r (det A_r d_r(P))^(1/r), and d_r(P) >= Hermite
          bound  =>  det A_2 <= 9.50,  det A_3 <= 14.18,  det A_4 = 1.
For each section class and each Gram shape X with tr(W conj X) <= 4 (finite list,
etype.certify.enumerate_X) we test exhaustively whether P contains (O_K^r, X).

usage: decide.py idx stage      stage in {sections, shapes2, shapes3, embed}
"""
import sys, os, time, pickle, itertools, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from etype.search import *
from etype.certify import *
from etype.represent import *

K = (5, 1); N, t = K; d = 19
idx = int(sys.argv[1]); stage = sys.argv[2]
OUT = os.path.join(RES, "decide_d19_%d" % idx)
os.makedirs(OUT, exist_ok=True)
A, JA = brick_walk(K, 4, steps=2500)[(6, 20)]
An, JAn = np.array([[int(A[i, j]) for j in range(8)] for i in range(8)]), np.array([[int(JA[i, j]) for j in range(8)] for i in range(8)])
BOUND = {2: 9.50, 3: 14.18}
VSTACK = pari('(a,b)->matconcat([a;b])')


def hgram(vecs):
    ex, cx = herm_gram(A, JA, K, vecs)
    return ex, cx


def saturate(vecs):
    """Z-basis of A cap (K-span of vecs): kernel of the dual functionals."""
    V = np.array([w for v in vecs for w in (v, JAn @ v)]).T            # 8 x 2r
    M = pari.matrix(8, V.shape[1], [int(x) for x in V.ravel()])
    # A cap Q-span(M) = kernel of (orthogonal complement functionals)
    Kc = pari.matker(pari.mattranspose(M) * A)                          # complement (rational)
    F = pari.mattranspose(Kc) * A
    S = pari.matkerint(F)
    return S


def ok_basis_of(S, r):
    """O_K-basis (as Z-vectors in A) of the rank-r saturated sublattice with Z-basis S."""
    Gs = pari.mattranspose(S) * A * S
    Js = pari.matinverseimage(S, JA * S)
    vecs, ex, cx = ok_basis(lll(Gs, Js)[0], lll(Gs, Js)[1], K, bound=None)
    return Gs, Js, ex, cx


def sections():
    """Isometry classes of primitive rank-2 and rank-3 sections of A below the bounds."""
    from etype.search import short_vectors
    out = {2: [], 3: []}
    # rank 3:  A_3 = w^perp, h(w,w) = det A_3 <= 14
    W = short_vectors(A, int(2 * BOUND[3]))
    for c in W.T:
        w = pari.matrix(8, 1, list(map(int, c)))
        S = pari.matkerint(VSTACK(pari.mattranspose(A * w), pari.mattranspose(A * JA * w)))   # b(.,w) = b(.,eta w) = 0
        Gs = pari.mattranspose(S) * A * S; Js = pari.matinverseimage(S, JA * S)
        Gs, Js = lll(Gs, Js)
        _add(out[3], Gs, Js)
    # rank 2: spans of (w1, w2), h(w1) <= 9.5, h(w2) <= 9.5/3 + 1.32*9.5 ; keep det <= 9.5
    W1 = short_vectors(A, 19); W2 = short_vectors(A, 31)
    for a in W1.T:
        for b in W2.T:
            ex, cx = herm_gram(A, JA, K, [list(map(int, a)), list(map(int, b))])
            D = np.linalg.det(cx).real
            if D < 0.5 or D > BOUND[2] + 1e-6:
                continue
            U = pari.matrix(8, 4, [int(x) for x in np.array([a, JAn @ a, b, JAn @ b]).T.ravel()])
            S = pari.matkerint(pari.mattranspose(pari.matker(pari.mattranspose(U))))   # A cap Q U
            Gs = pari.mattranspose(S) * A * S; Js = pari.matinverseimage(S, JA * S)
            Gs, Js = lll(Gs, Js)
            if float(pari.matdet(Gs)) > d ** 2 * BOUND[2] ** 2 + 1e-6:
                continue
            _add(out[2], Gs, Js)
    return out


def _add(lst, G, J):
    for (G2, J2) in lst:
        if pari.matdet(G2) == pari.matdet(G) and herm_isom(G, J, G2, J2, K):
            return
    lst.append((G, J))


if stage == "sections":
    t0 = time.time()
    sec = sections()
    for r in (2, 3):
        print("rank-%d section classes:" % r, [round(float(pari.matdet(G)) ** 0.5 / d ** (r / 2), 3) for G, J in sec[r]], flush=True)
    pickle.dump({r: [(str(G), str(J)) for G, J in sec[r]] for r in sec}, open(os.path.join(OUT, "sections.pkl"), "wb"))
    print("%.0fs" % (time.time() - t0))

elif stage.startswith("shapes"):
    r = int(stage[-1])
    sec = pickle.load(open(os.path.join(OUT, "sections.pkl"), "rb"))[r]
    allshapes = []
    for (Gs, Js) in sec:
        Gs, Js = pari(Gs), pari(Js)
        vecs, ex, cx = ok_basis(Gs, Js, K, bound=int(minimum_exact(Gs)) + 12)
        t0 = time.time()
        R = enumerate_X(cx, K, 4.0)
        good = []
        for Tv, dg, off in R:
            G, J = X_lattice(dg, off, K)
            if count_short(G, 2) == 0:
                good.append((Tv, dg, off))
        print("section det %.3f: %d shapes (%.0fs)" % (np.linalg.det(cx).real, len(good), time.time() - t0), flush=True)
        allshapes += good
    pickle.dump(allshapes, open(os.path.join(OUT, "shapes%d.pkl" % r), "wb"))

elif stage == "classes":
    # all shapes -> pair-type filter (complete table of P's minimal vectors) -> O_K-isometry classes
    from etype.search import short_vectors
    leech = load("leech_d19.json"); import random; random.seed(0); random.shuffle(leech)
    P, JP = leech[idx]["G"], leech[idx]["J"]
    tf = os.path.join(RES, "types22_d19_%d.pkl" % idx)
    tp = pickle.load(open(tf, "rb"))
    allowed = set((int(u * d), int(v * d)) for (hx, hy, u, v) in tp)
    shapes = []
    for r in (2, 3):
        shapes += [(r,) + s for s in pickle.load(open(os.path.join(OUT, "shapes%d.pkl" % r), "rb"))]
    for part in range(3):
        shapes += [(4, s[0], s[1], s[2]) for s in pickle.load(open(os.path.join(RES, "X_d19_T4_part%d_3.pkl" % part), "rb"))[2]]
    print("%d shapes" % len(shapes), flush=True)
    t0 = time.time(); groups = collections.defaultdict(list); nf = 0
    for s in shapes:
        G, J = X_lattice(s[2], s[3], K)
        V = short_vectors(G, 4)
        if V.shape[1]:
            V = np.concatenate([V, -V], axis=1)
            Gn = np.array([[int(G[i, j]) for j in range(len(G))] for i in range(len(G))]); Jn = np.array([[int(J[i, j]) for j in range(len(G))] for i in range(len(G))])
            a_ = V.T @ Gn @ V; c_ = V.T @ Gn @ (Jn @ V)
            du, dv = 2 * N * a_ - t * c_, 2 * c_ - t * a_
            if not all((int(x), int(y)) in allowed for x, y in zip(du.ravel(), dv.ravel())):
                nf += 1; continue
        key = (s[0], int(pari.matdet(G))) + tuple(int(c) for c in pari.qfrep(G, 10, 1))
        groups[key].append((G, J))
    print("pair-type filter removes %d; %d invariant groups (%.0fs)" % (nf, len(groups), time.time() - t0), flush=True)
    classes = []
    for key, lst in groups.items():
        reps = []
        for G, J in lst:
            if not any(herm_isom(G, J, G2, J2, K) for G2, J2 in reps):
                reps.append((G, J))
        classes += [(key[0], str(G), str(J)) for G, J in reps]
    print("%d O_K-isometry classes (%.0fs)" % (len(classes), time.time() - t0), flush=True)
    print(collections.Counter(c[0] for c in classes))
    pickle.dump(classes, open(os.path.join(OUT, "classes.pkl"), "wb"))

elif stage == "embed":
    part, nparts = int(sys.argv[3]), int(sys.argv[4])
    leech = load("leech_d19.json"); import random; random.seed(0); random.shuffle(leech)
    P, JP = leech[idx]["G"], leech[idx]["J"]
    o, gens = herm_aut(P, JP, K)
    host = MinHost(P, JP, K, gens=gens)
    freq = code_freq(host)
    print("|Aut_h(P)| = %d, %d orbits on the 196560 minimal vectors" % (o, len(host.reps)), flush=True)
    allcl = pickle.load(open(os.path.join(OUT, "classes.pkl"), "rb"))
    maxcost = float(sys.argv[5]) if len(sys.argv) > 5 else float("inf")
    costs = pickle.load(open(os.path.join(OUT, "costs.pkl"), "rb"))
    order = sorted((i for i in range(len(allcl)) if costs[i] is not None and costs[i] <= maxcost), key=lambda i: costs[i])
    todo = order[part::nparts]
    logf = os.path.join(OUT, "embed_log_%d_%d.txt" % (part, nparts))
    done = set()
    if os.path.exists(logf):
        done = {int(l.split()[0]) for l in open(logf) if l.strip()}
    todo = [i for i in todo if i not in done]
    print("%d classes left in this part" % len(todo), flush=True)
    LOG = open(logf, "a")
    stats = collections.Counter(); kills = []; skipped = []
    t0 = time.time()
    for n, ci in enumerate(todo):
        r, G, J = allcl[ci]
        G, J = pari(G), pari(J)
        ok, sol = embeds_fast(host, G, J, K, freq)
        stats[(r, ok)] += 1
        LOG.write("%d %d %s\n" % (ci, r, ok)); LOG.flush()
        if ok:
            kills.append((r, str(G), str(J), [list(map(int, host.V[:, i])) for i in sol]))
            print("  EMBEDS: rank %d, Gram %s" % (r, G), flush=True)
        elif ok is None:
            skipped.append((r, str(G), str(J)))
        if n % 50 == 0:
            print("  %d / %d  %s  (%.0fs)" % (n, len(todo), dict(stats), time.time() - t0), flush=True)
    print("RESULT part %d:" % part, dict(stats), flush=True)
    pickle.dump((kills, skipped), open(os.path.join(OUT, "embed_%d_%d.pkl" % (part, nparts)), "wb"))
