"""Exact decision of  min(A (x)_{O_K} P) <= 8   (A rank-m brick, P Hermitian Leech structure).

z in A (x) P of tensor rank r lies in A_r (x) P for a primitive rank-r section A_r of A, and
     h(z,z) = tr(conj(W_r) X),   W_r = Gram of an O_K-basis of A_r,  X = Gram of x_1..x_r in P.
  r = 1 : h = h(a,a) h(p,p) >= 2 hmin(A) > 4  (needs hmin(A) >= 3).
  r >= 2: Coulangeon, Acta Arith. 92 (2000), Prop. 3.2:  4 >= r (det A_r * d_r(P))^(1/r);
          d_r(P) >= Hermite bound of the 2r-dimensional Z-sections (minimum 4) of P.
  r = m : A_m = A.
For every section and every Gram shape X with tr(conj(W) X) <= 4 (etype.certify.enumerate_X,
finite and exact) we search exhaustively for an O_K-linear isometric embedding (O_K^r, X) -> P.
Any embedding gives an explicit vector of norm <= 8, which is constructed and checked.

Requirements, checked at run time: A free over O_K and every section used is free
(automatic when h_K = 1); otherwise the section is reported, not silently skipped.
"""
import os, time, json, pickle, itertools, collections, math
import numpy as np
from .core import pari, _K, disc, minimum_exact, count_short, lll, _kron
from .certify import ok_basis, herm_gram, enumerate_X, X_lattice
from .represent import herm_isom, herm_aut, MinHost, code_freq, embeds_fast
from .search import short_vectors

VSTACK = pari('(a,b)->matconcat([a;b])')
GAMMA = {2: 2 ** 0.5, 3: (64 / 3) ** (1 / 6)}          # gamma_4, gamma_6 (exact values)


def constants(K, hminA, hminP=2, T=4.0):
    N, t = _K(K); D = disc(K)
    hK = int(pari("quadclassunit(%d)[1]" % (-D)))
    dP = {r: math.sqrt((2 * hminP / GAMMA[r]) ** (2 * r) / D ** r) for r in (2, 3)}
    B = {r: (T / r) ** r / dP[r] for r in (2, 3)}
    mu = (D + 1) ** 2 / (16 * D) if t == 1 else (1 + N) / 4      # covering radius^2 of O_K
    cK = 1.0 if hK == 1 else ideal_class_factor(K)              # min norm of an integral ideal, worst class
    assert hminA * hminP > T, "rank-1 tensors are not excluded"
    return dict(D=D, hK=hK, dP=dP, B=B, mu=mu, cK=cK)


def ideal_class_factor(K):
    """max over ideal classes c of  min{ N(b) : b integral ideal in c }.
    A rank-1 lattice a*w (a an ideal) of determinant D1 contains a vector of h <= c_K * D1."""
    N, t = _K(K); D = disc(K)
    bnf = pari.bnfinit(pari("x^2 - %d*x + %d" % (t, N)))
    best = {}
    bound = int(math.sqrt(D / 3)) + 2                            # every class has an ideal of norm <= sqrt(D/3)
    L = pari.ideallist(bnf, bound)
    for n in range(1, bound + 1):
        for I in L[n - 1]:
            c = tuple(int(x) for x in pari.bnfisprincipal(bnf, I, 0))
            best.setdefault(c, n)
    hK = int(bnf.bnf_get_no())
    assert len(best) == hK, "class representatives not found"
    return float(max(best.values()))


def _np(M):
    return np.array([[int(M[i, j]) for j in range(len(M))] for i in range(len(M))], dtype=np.int64)


def sections(A, JA, K, C, hminA):
    """Isometry classes of primitive rank-2, rank-3 sections of the rank-4 brick A with
    Hermitian determinant below the bounds C['B']."""
    D = C["D"]; An, JAn = _np(A), _np(JA)
    out = {2: [], 3: []}
    seen = collections.defaultdict(list)
    def add(lst, G, J):
        key = (len(G), int(pari.matdet(G))) + tuple(int(c) for c in pari.qfrep(G, 16, 1))
        for G2, J2 in seen[key]:
            if herm_isom(G, J, G2, J2, K):
                return
        seen[key].append((G, J)); lst.append((G, J))
    # rank 3: S = w^perp, det S = det(O_K-span of w, saturated) <= h(w,w) * cK
    # NB: loop bodies live in functions -- cypari2 frees PARI objects LIFO, and rebinding
    # loop variables at module/loop level makes the PARI stack climb without bound.
    def perp(c):
        w = pari.matrix(8, 1, [int(x) for x in c])
        S = pari.matkerint(pari.mattranspose(pari.matconcat([A * w, A * JA * w])))
        Gs, Js = lll(pari.mattranspose(S) * A * S, pari.matinverseimage(S, JA * S))
        if float(pari.matdet(Gs)) <= D ** 3 * C["B"][3] ** 2 + 1e-6:
            add(out[3], Gs, Js)
    W = short_vectors(A, int(2 * C["B"][3] * C["cK"]) + 1)
    for c in W.T:
        perp(c)
    # rank 2: S saturated; w1 minimal in S, w2 = first vector of S outside K w1.
    # Minkowski II on the 4-dim Z-lattice of S (det D^2 det(S)^2, minimum >= 2 hmin(A)):
    #   lam1 <= (gamma_4^4 det)^(1/4),  lam3 <= sqrt(gamma_4^4 det / lam1^2),  and w2 has norm <= lam3.
    detZ = 4 * D ** 2 * C["B"][2] ** 2
    l1 = detZ ** 0.25
    l3 = max(math.sqrt(detZ / (2 * hminA) ** 2),
             2 * (C["cK"] ** 2 * C["B"][2] / hminA + C["mu"] * l1 / 2))    # projection argument (Z-norm)
    W1 = short_vectors(A, int(l1) + 1).T; W2 = short_vectors(A, int(l3) + 1).T
    N_, t_ = _K(K)
    def hcoef(X, Y):                                   # h(x,y) for rows of X, Y: (u, v) with h = u + v eta, times D
        a = X @ An @ Y.T; c = X @ An @ (JAn @ Y.T)
        return 2 * N_ * a - t_ * c, 2 * c - t_ * a
    def sat2(a, b):
        U8 = pari.matrix(8, 4, [int(x) for x in np.array([a, JAn @ a, b, JAn @ b]).T.ravel()])
        S = pari.matkerint(pari.mattranspose(pari.matker(pari.mattranspose(U8))))
        Gs, Js = lll(pari.mattranspose(S) * A * S, pari.matinverseimage(S, JA * S))
        if float(pari.matdet(Gs)) <= D ** 2 * C["B"][2] ** 2 + 1e-6:
            add(out[2], Gs, Js)
    u11, _ = hcoef(W1, W1); h11 = np.diag(u11) * 1.0 / D            # h(w,w) real: u/D (v = 0)
    hh = lambda U, V: (U * U + t_ * U * V + N_ * V * V) / D ** 2      # |u + v eta|^2 / D^2 ... norm of h
    nW2 = np.einsum("ij,jk,ik->i", W2, An, W2) / 2.0
    for i, a in enumerate(W1):
        U, V = hcoef(a[None, :], W2)
        det = h11[i] * nW2 - (U[0] ** 2 + t_ * U[0] * V[0] + N_ * V[0] ** 2) / D ** 2
        # span(w1, w2) has index <= cK^2 in S when w1 is minimal in S and the projection of
        # w2 is minimal in the projection of S, so det(span) <= cK^2 * det(S) <= cK^2 * B.
        for j in np.nonzero((det > 1e-6) & (det <= C["B"][2] * C["cK"] ** 2 + 1e-6))[0]:
            sat2(a, W2[j])
    return out


def _admissible(K):
    def keep(Tv, dg, off):
        G, J = X_lattice(dg, off, K)
        return count_short(G, 2) == 0                  # no vector of Hermitian norm 1
    return keep


def _shape_worker(args):
    W, K, T, part, nparts = args[:5]
    out = args[5] if len(args) > 5 else None
    if out and os.path.exists(out):
        return pickle.load(open(out, "rb"))
    good = enumerate_X(W, K, T, part=(part, nparts), keep=_admissible(K))
    if out:
        pickle.dump(good, open(out + ".tmp", "wb")); os.replace(out + ".tmp", out)
    return good


def all_shapes(A, JA, K, secs, nproc=3, T=4.0, log=print, partdir=None):
    """[(r, Tv, diag, off)] for the rank-2/3 sections and the full brick; non-free sections
    are returned separately."""
    import multiprocessing as mp
    shapes, nonfree = [], []
    for r in (2, 3):
        for Gs, Js in secs[r]:
            try:
                v, ex, cx = ok_basis(Gs, Js, K, bound=int(minimum_exact(Gs)) + 12)
            except RuntimeError:
                cx, ideals = steinitz_gram(Gs, Js, K)                # non-free: Steinitz pseudo-basis
                nonfree.append((r, str(Gs), str(Js), ideals))
            shapes += [(r,) + s for s in _shape_worker((cx, K, T, 0, 1))]
    log("rank-2/3 sections: %d shapes, %d non-free sections" % (len(shapes), len(nonfree)))
    try:
        v, ex, W = ok_basis(A, JA, K)
    except RuntimeError:                                          # non-free brick: Steinitz pseudo-basis
        W, ideals = steinitz_gram(A, JA, K); v = None
        log("brick is not free: Steinitz ideals %s" % ideals)
    t0 = time.time()
    NP = 48
    pdir = (partdir or "/tmp") ; os.makedirs(pdir, exist_ok=True)
    with mp.get_context("fork").Pool(nproc, maxtasksperchild=1) as pool:
        parts = pool.map(_shape_worker, [(W, K, T, p, NP, os.path.join(pdir, "part_%02d.pkl" % p)) for p in range(NP)], chunksize=1)
    s4 = [s for p in parts for s in p]
    log("rank-4: %d shapes (%.0fs)" % (len(s4), time.time() - t0))
    shapes += [(4,) + s for s in s4]
    return shapes, nonfree, v


def classes(shapes, K, allowed, log=print):
    """pair-type filter against P, then O_K-isometry classes.  Memory-light: shapes stay plain
    Python data; PARI objects live only inside the scoped helpers; representatives are kept
    as strings.  Returns [(r, G, J)] with PARI G, J for the representatives."""
    N, t = _K(K)
    groups = collections.defaultdict(list); nf = [0]
    def key_of(s):
        G, J = X_lattice(s[2], s[3], K)
        V = short_vectors(G, 4)
        if V.shape[1]:
            V = np.concatenate([V, -V], axis=1)
            Gn, Jn = _np(G), _np(J)
            a = V.T @ Gn @ V; c = V.T @ Gn @ (Jn @ V)
            codes = set(((2 * N * a - t * c) * 4096 + (2 * c - t * a)).ravel().tolist())
            if not codes <= allowed:
                return None
        return (s[0], int(pari.matdet(G))) + tuple(int(c) for c in pari.qfrep(G, 10, 1))
    for s in shapes:
        k = key_of(s)
        if k is None:
            nf[0] += 1
        else:
            groups[k].append((s[2], s[3]))
    def reps_of(lst):
        reps = []                                   # strings
        for dg, off in lst:
            G, J = X_lattice(dg, off, K)
            if not any(herm_isom(G, J, pari(G2), pari(J2), K) for G2, J2 in reps):
                reps.append((str(G), str(J)))
        return reps
    out = []
    for key, lst in groups.items():
        out += [(key[0], G, J) for G, J in reps_of(lst)]
    log("pair-type filter removes %d; %d O_K-isometry classes %s" % (nf[0], len(out), dict(collections.Counter(c[0] for c in out))))
    return [(r, pari(G), pari(J)) for r, G, J in out]


_G = {}
def _embed_worker(ci):
    host, fq, K, cl = _G["host"], _G["fq"], _G["K"], _G["cl"]
    r, G, J = cl[ci]
    ok, sol = embeds_fast(host, G, J, K, fq)
    return ci, ok


def brick_shapes(K, A, JA, cache, nproc=3, log=print):
    """Sections and Gram shapes of a brick (independent of P), cached on disk."""
    if os.path.exists(cache):
        return pickle.load(open(cache, "rb"))
    hminA = minimum_exact(A) // 2
    C = constants(K, hminA)
    log("K = Q(sqrt(-%d)), h_K = %d, hmin(A) = %d; bounds det A_2 <= %.2f, det A_3 <= %.2f"
        % (C["D"], C["hK"], hminA, C["B"][2], C["B"][3]))
    t0 = time.time()
    sfile = cache + ".sections"
    if os.path.exists(sfile):
        secs = {r: [(pari(G), pari(J)) for G, J in v] for r, v in pickle.load(open(sfile, "rb")).items()}
    else:
        secs = sections(A, JA, K, C, hminA)
        pickle.dump({r: [(str(G), str(J)) for G, J in v] for r, v in secs.items()}, open(sfile, "wb"))
    log("sections: %d of rank 2, %d of rank 3 (%.0fs)" % (len(secs[2]), len(secs[3]), time.time() - t0))
    shapes, nonfree, vecsA = all_shapes(A, JA, K, secs, nproc=nproc, log=log, partdir=cache + ".parts")
    res = (shapes, nonfree, vecsA)
    pickle.dump(res, open(cache, "wb"))
    return res


def decide(K, A, JA, P, JP, outdir, nproc=3, maxcost=1e6, log=print, cache=None):
    """Returns ('KILLED', norm, z) / ('NO EMBEDDING', report) / ('INCOMPLETE', report)."""
    import multiprocessing as mp
    os.makedirs(outdir, exist_ok=True)
    N, t = _K(K)
    shapes, nonfree, vecsA = brick_shapes(K, A, JA, cache or os.path.join(outdir, "brick_shapes.pkl"), nproc=nproc, log=log)
    log("%d shapes, %d non-free sections" % (len(shapes), len(nonfree)))
    o, gens = herm_aut(P, JP, K)
    host = MinHost(P, JP, K, gens=gens); fq = code_freq(host)
    allowed = set()
    for i in host.reps:
        allowed |= set(np.unique(host.code(np.array([i]), np.arange(host.V.shape[1]))[0]).tolist())
    log("|Aut(P)| = %d, %d orbits on minimal vectors, %d product types" % (o, len(host.reps), len(allowed)))
    cl = classes(shapes, K, allowed, log=log)
    costs = []
    for r, G, J in cl:
        res = embeds_fast(host, G, J, K, fq, prepare_only="cost")
        costs.append(None if isinstance(res[1], str) else res[0])
    order = sorted((i for i in range(len(cl)) if costs[i] is not None and costs[i] <= maxcost), key=lambda i: costs[i])
    hard = [i for i in range(len(cl)) if costs[i] is not None and costs[i] > maxcost]
    nospan = [i for i in range(len(cl)) if costs[i] is None]
    log("embedding search: %d classes (cost <= %.0e), %d harder, %d not spanned by minimal vectors"
        % (len(order), maxcost, len(hard), len(nospan)))
    _G.update(host=host, fq=fq, K=K, cl=cl)
    hit = None; t0 = time.time(); n = 0
    with mp.get_context("fork").Pool(nproc) as pool:
        for ci, ok in pool.imap_unordered(_embed_worker, order, chunksize=4):
            n += 1
            if ok:
                hit = ci; pool.terminate(); break
            if n % 500 == 0:
                log("   %d / %d classes, no embedding (%.0fs)" % (n, len(order), time.time() - t0))
    if hit is not None:
        z, nz = explicit_vector(K, A, JA, P, JP, vecsA, host, fq, cl[hit])
        log("KILLED: class %d embeds; explicit vector of norm %d" % (hit, nz))
        pickle.dump({"class": hit, "norm": nz, "z": z}, open(os.path.join(outdir, "kill.pkl"), "wb"))
        return "KILLED", nz, z
    rep = dict(classes=len(cl), searched=len(order), hard=len(hard), nospan=len(nospan), nonfree=len(nonfree))
    pickle.dump({"report": rep, "hard": [cl[i] for i in hard], "nospan": [cl[i] for i in nospan], "nonfree": nonfree},
                open(os.path.join(outdir, "open.pkl"), "wb"))
    status = "NO EMBEDDING" if not (hard or nospan or nonfree) else "INCOMPLETE"
    log("%s: %s" % (status, rep))
    return status, rep


def explicit_vector(K, A, JA, P, JP, vecsA, host, fq, cls):
    """Build z = sum b_i (x) F(e_i) in A (x) P from an embedding of a rank-4 shape class and
    return (coordinates in the HNF basis of A (x) P, exact norm).  For rank < 4 classes the
    section basis would be needed; the embedding itself is then reported."""
    N, t = _K(K); D = disc(K)
    r, G, J = cls
    F = embeds_fast(host, G, J, K, fq, prepare_only="map")
    n1, n2 = len(A), len(P)
    Pm = -(_kron(JA, pari.matid(n2)) + _kron(pari.matid(n1), JP) - t * pari.matid(n1 * n2)) * _kron(2 * JA - t * pari.matid(n1), pari.matid(n2))
    B = pari.mathnf(Pm); G96 = pari.mattranspose(B) * _kron(A, P) * B / D ** 2
    if r != len(vecsA):
        return None, -1
    xs = [[int(F[l, 2 * i]) for l in range(n2)] for i in range(r)]
    z = sum((Pm * pari.matrix(n1 * n2, 1, [a * b for a in vecsA[i] for b in xs[i]]) for i in range(r)), pari.matrix(n1 * n2, 1)) / D
    c = pari.matsolve(B, z * D)
    assert pari.denominator(c) == 1
    return [int(c[k, 0]) for k in range(len(B))], int((pari.mattranspose(c) * G96 * c)[0, 0])


def steinitz_gram(Gs, Js, K):
    """Complex Gram W of a Steinitz pseudo-basis (b_1..b_r) of a (possibly non-free) O_K-lattice
    S given by Z-data (Gs, Js):  S = O_K b_1 + ... + O_K b_{r-1} + a b_r.
    Every z in S (x) P is sum b_i (x) x_i with x_i in P (x_r in aP), so the shapes for W
    contain all shapes of S (superset; kills are verified explicitly)."""
    N, t = _K(K)
    n = len(Gs); r = n // 2
    nf = pari.nfinit(pari("y^2 - %d*y + %d" % (t, N)))
    # choose a K-basis v_1..v_r of S among short vectors, express a Z-basis of S in it
    V = short_vectors(Gs, minimum_exact(Gs) + 40).T
    Jn = _np(Js)
    basis = []
    for v in V:
        cand = basis + [v]
        B = np.array([w for u in cand for w in (u, Jn @ u)], dtype=float)
        if np.linalg.matrix_rank(B, tol=1e-8) == 2 * len(cand):
            basis = cand
            if len(basis) == r:
                break
    Bv = np.array([w for u in basis for w in (u, Jn @ u)], dtype=float).T          # n x n, columns v_i, eta v_i
    C = np.linalg.solve(Bv, np.eye(n))                                               # coords of e_k
    # e_k = sum_i (C[2i,k] + C[2i+1,k] eta) v_i   ->  pseudo-matrix columns in K^r
    from fractions import Fraction
    cols = []
    for k in range(n):
        col = []
        for i in range(r):
            a = Fraction(C[2 * i, k]).limit_denominator(10**6); b = Fraction(C[2 * i + 1, k]).limit_denominator(10**6)
            col.append(pari("Mod(%s + (%s)*y, y^2 - %d*y + %d)" % (a, b, t, N)))
        cols.append(col)
    M = pari.matrix(r, n, [cols[k][i] for i in range(r) for k in range(n)])
    H = pari.nfhnf(nf, [M, [1] * n])
    St = pari.rnfsteinitz(nf, H)
    Bm, Id = St[0], St[1]
    # normalise the last ideal: a = lambda * c with c integral of minimal norm in its class,
    # b_r -> lambda b_r, so that x_r ranges over c P, a subset of P.
    bnf = pari.bnfinit(pari("y^2 - %d*y + %d" % (t, N)))
    a = Id[r - 1]
    cls = pari.bnfisprincipal(bnf, a, 0)
    D_ = 4 * N - t * t
    best = None
    L = pari.ideallist(bnf, int(math.sqrt(D_ / 3)) + 2)
    for n_ in range(1, len(L) + 1):
        for I in L[n_ - 1]:
            if pari.bnfisprincipal(bnf, I, 0) == cls:
                best = I; break
        if best is not None:
            break
    lam = pari.bnfisprincipal(bnf, pari.idealdiv(bnf, a, best), 3)[1]                 # a = lam * c (flag 1|2: generator)
    Id = list(Id[:r - 1]) + [best]
    assert pari.idealnorm(nf, best).denominator() == 1
    # pseudo-basis vectors b_j = sum_i Bm[i,j] v_i  (in K-coordinates); Hermitian Gram via v's
    ex, cx = herm_gram(Gs, Js, K, [list(map(int, u)) for u in basis])
    Hv = np.array(cx)
    et = complex(t / 2, math.sqrt(4 * N - t * t) / 2)
    def cval(z):                                        # integral basis of nf is [1, y], y = eta
        col = pari.nfalgtobasis(nf, z)
        return float(col[0]) + float(col[1]) * et
    Bc = np.array([[cval(Bm[i, j]) for j in range(r)] for i in range(r)])
    Bc[:, r - 1] *= cval(lam)                                     # b_r -> lambda b_r
    W = Bc.T @ Hv @ np.conj(Bc)                                   # h(b_j, b_k) = sum B_ij conj(B_lk) h(v_i, v_l)
    ideals = [str(Id[j]) for j in range(r)]
    return W, ideals
