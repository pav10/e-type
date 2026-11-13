"""Exact analysis of short vectors in  A (x)_{O_K} P   (A Hermitian unimodular).

Every z in A (x) P is  z = sum_i b_i (x) x_i  for a fixed O_K-basis (b_i) of a
primitive sublattice A' of A and x_i in P; if z has rank r we may take rank A' = r
and x_1..x_r independent.  With  W_ij = h_A(b_i,b_j),  X_ij = h_P(x_i,x_j):

        h(z,z) = sum_ij W_ij X_ij = tr(W conj(X)).

X is a positive definite Hermitian matrix with entries in Diff^-1 = O_K/sqrt(-d)
(P is a Z-unimodular lattice with O_K-structure), integer diagonal, and every
value of X on O_K^r is a Hermitian norm of P, hence >= hmin(P) (= 2 for Leech).

So:  min(A (x) P) <= T   iff   some X in the FINITE set
        S(W, T) = { X > 0 admissible :  tr(W conj X) <= T }
is represented by P.  `enumerate_X` lists S(W,T) (independent of P!);
the representation test is done against the Leech structure afterwards.
This is the method of Coulangeon (Acta Arith. 92, proof of Thm 4.2.1, rank 4)
and Coulangeon-Nebe (Thm 5.3), made systematic.
"""
import itertools, math
import numpy as np
from fractions import Fraction
from .core import pari, _K, disc, herm_to_z, minimum_exact, lll


def eta_c(K):
    N, t = _K(K)
    return complex(t / 2, math.sqrt(4 * N - t * t) / 2)


def _fp(G, m0, R, out, limit=10**7):
    """all integer m with (m-m0)^T G (m-m0) <= R  (Fincke-Pohst, small dimension)."""
    n = len(G)
    Lc = np.linalg.cholesky(G).T                       # G = Lc^T Lc, upper triangular
    m = [0] * n
    def rec(i, rem):
        if i < 0:
            out.append(tuple(m)); return
        s = sum(Lc[i, j] * (m[j] - m0[j]) for j in range(i + 1, n)) / Lc[i, i]
        c = m0[i] - s
        rad = math.sqrt(max(rem, 0.0)) / Lc[i, i]
        for v in range(math.ceil(c - rad - 1e-9), math.floor(c + rad + 1e-9) + 1):
            m[i] = v
            rec(i - 1, rem - (Lc[i, i] * (v - c)) ** 2)
            if len(out) > limit:
                raise RuntimeError("too many points")
    rec(n - 1, R + 1e-9)


def enumerate_X(W, K, T, hmin=2, part=(0, 1), verbose=False):
    """All positive definite X (r x r) with entries (a + b*eta)/sqrt(-d), integer
    diagonal >= hmin, value >= hmin on a set of small test vectors, and
    tr(W conj X) <= T.   W: complex Hermitian r x r (Gram of the basis of A').
    Returns a list of (T(X), diag (ints), off = {(i,j): (a,b)} for i<j)."""
    N, t = _K(K); d = disc(K); et = eta_c(K); delta = complex(0, math.sqrt(d))
    r = len(W)
    H = np.linalg.inv(np.conj(np.array(W, dtype=complex)))     # T = tr(H^-1 X)
    L = np.linalg.cholesky(H)                                  # H = L L^*
    # real coordinates of Diff^-1: z = (a + b eta)/delta
    e1, e2 = 1 / delta, et / delta
    small = [0, 1, -1, et, -et, np.conj(et), -np.conj(et), 1 + et, -1 - et, 1 - et, et - 1]
    res, stats = [], [0] * (r + 1)
    # leaf filter: values of X on all coefficient vectors with entries in a small set
    sm2 = [0, 1, -1, et, -et, 1 - et, et - 1]
    NT = np.array([n for n in itertools.product(sm2, repeat=r) if any(n)], dtype=complex).T    # r x M
    Winv_c = np.linalg.inv(H)
    counter = [0]

    def ok_min(X, k):
        # values n^T X conj(n) for n supported on <= 2 coordinates, last one = k-1
        for i in range(k - 1):
            for a in small[1:]:
                for b in small[1:3]:
                    v = (abs(a) ** 2 * X[i, i] + abs(b) ** 2 * X[k - 1, k - 1]
                         + 2 * (a * np.conj(b) * X[i, k - 1]).real)
                    if v < hmin - 1e-6:
                        return False
        return True

    def rec(k, X, diag, off):
        stats[k] += 1
        if k == r:
            vals = np.einsum("im,ij,jm->m", NT, X, NT.conj()).real
            if vals.min() < hmin - 1e-6:
                return
            Tv = float(np.trace(Winv_c @ X).real)
            res.append((Tv, tuple(diag), dict(off)))
            return
        if k == 2:                                   # split the search for parallel runs
            counter[0] += 1
            if counter[0] % part[1] != part[0]:
                return
        Lp = L[:k, :k]; l = np.conj(L[k, :k]); lam = L[k, k].real
        Y = np.linalg.solve(Lp, np.linalg.solve(Lp, X).conj().T).conj().T      # L'^-1 X L'^-*
        rho = T - np.trace(Y).real
        if rho < -1e-9:
            return
        Xi = np.linalg.inv(X)
        c0 = Lp @ (Y @ l)
        # real lattice coordinates m = (a_1, b_1, ..., a_k, b_k), c_i = a_i e1 + b_i e2
        B = np.zeros((k, 2 * k), dtype=complex)
        for i in range(k):
            B[i, 2 * i], B[i, 2 * i + 1] = e1, e2
        G = (B.conj().T @ Xi @ B).real
        Br = np.vstack([B.real, B.imag]); c0r = np.concatenate([c0.real, c0.imag])
        m0 = np.linalg.solve(Br, c0r)
        pts = []
        _fp(G, m0, rho * lam * lam, pts)
        lYl = (np.conj(l) @ Y @ l).real
        u = np.conj(l) @ np.linalg.inv(Lp)                                     # l^* L'^-1
        for m in pts:
            c = B @ np.array(m, dtype=float)
            dc = c - c0
            q = (np.conj(dc) @ Xi @ dc).real
            lin = 2 * (u @ c).real - lYl
            lo, hi = q + lin, rho * lam * lam + lin
            for x in range(max(hmin, math.ceil(lo - 1e-7)), math.floor(hi + 1e-7) + 1):
                if x - (np.conj(c) @ Xi @ c).real < 1e-7:                      # must stay definite
                    continue
                X2 = np.zeros((k + 1, k + 1), dtype=complex)
                X2[:k, :k] = X; X2[:k, k] = c; X2[k, :k] = np.conj(c); X2[k, k] = x
                if not ok_min(X2, k + 1):
                    continue
                off2 = dict(off)
                for i in range(k):
                    off2[(i, k)] = (m[2 * i], m[2 * i + 1])
                rec(k + 1, X2, diag + [x], off2)

    for x in range(hmin, math.floor(T * H[0, 0].real + 1e-9) + 1):
        rec(1, np.array([[x]], dtype=complex), [x], {})
    if verbose:
        print("   nodes per level:", stats[1:], flush=True)
    return res


def X_exact(diag, off, K):
    """Hermitian matrix as pairs (u, v) = u + v*eta (Fractions) from the integer data."""
    N, t = _K(K); d = disc(K)
    r = len(diag)
    def over_delta(a, b):            # (a + b eta)/sqrt(-d),  sqrt(-d) = 2 eta - t
        # 1/delta = -delta/d = (t - 2 eta)/d ;  (a + b eta)(t - 2 eta) = a t + (b t - 2a) eta - 2b eta^2
        #  eta^2 = t eta - N
        u = a * t + 2 * b * N
        v = b * t - 2 * a - 2 * b * t
        return (Fraction(u, d), Fraction(v, d))
    conj = lambda p: (p[0] + t * p[1], -p[1])
    M = [[None] * r for _ in range(r)]
    for i in range(r):
        M[i][i] = (Fraction(diag[i]), Fraction(0))
        for j in range(i + 1, r):
            M[i][j] = over_delta(*off[(i, j)])
            M[j][i] = conj(M[i][j])
    return M


def X_lattice(diag, off, K):
    """Z-Gram (Tr o X) and J of the Hermitian lattice (O_K^r, X); None if not an
    even integral lattice of minimum >= 4 (i.e. Hermitian minimum >= 2)."""
    G, J = herm_to_z(X_exact(diag, off, K), K)
    if pari.denominator(G) != 1:
        return None
    return G, J


def herm_gram(G, J, K, vecs):
    """Hermitian Gram h(v_i, v_j) of integer column vectors, as (u, v) Fractions and as complex."""
    N, t = _K(K); d = disc(K); et = eta_c(K)
    n = len(G)
    Gn = [[int(G[i, j]) for j in range(n)] for i in range(n)]
    Jn = [[int(J[i, j]) for j in range(n)] for i in range(n)]
    Gn, Jn = np.array(Gn, dtype=object), np.array(Jn, dtype=object)
    V = np.array(vecs, dtype=object).T                       # n x r
    a = V.T @ Gn @ V; c = V.T @ Gn @ (Jn @ V)                # b(x,y), b(x,Jy)
    r = V.shape[1]
    ex = [[(Fraction(2 * N * a[i, j] - t * c[i, j], d), Fraction(2 * c[i, j] - t * a[i, j], d)) for j in range(r)] for i in range(r)]
    cx = np.array([[float(p[0]) + float(p[1]) * et for p in row] for row in ex])
    return ex, cx


def ok_basis(G, J, K, bound=None):
    """An O_K-basis of a free O_K-lattice (Z-data G, J) made of short vectors:
    v_1..v_m with {v_i, eta v_i} a Z-basis.  Returns (vectors, exact Gram, complex Gram)."""
    from .search import short_vectors
    m = len(G) // 2
    mn = minimum_exact(G)
    Jn = np.array([[int(J[i, j]) for j in range(len(G))] for i in range(len(G))], dtype=np.int64)
    Gn = np.array([[int(G[i, j]) for j in range(len(G))] for i in range(len(G))], dtype=np.int64)
    for extra in (2, 6, 12, 20, 40):
        V = short_vectors(G, (bound or mn) + extra).T
        V = V[np.argsort(np.einsum("ki,ij,kj->k", V, Gn, V), kind="stable")][:60]
        for comb in itertools.combinations(range(len(V)), m):
            M = np.array([w for i in comb for w in (V[i], Jn @ V[i])], dtype=float)
            if abs(abs(np.linalg.det(M)) - 1) < 1e-6:
                vecs = [list(map(int, V[i])) for i in comb]
                ex, cx = herm_gram(G, J, K, vecs)
                return vecs, ex, cx
    raise RuntimeError("no O_K-basis among short vectors")
