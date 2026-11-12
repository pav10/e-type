"""Core lattice tools (Python + PARI via cypari2 + fpylll).

Conventions
-----------
* A Z-lattice is an integral symmetric positive definite Gram matrix G
  (PARI t_MAT).  min(L) = smallest value of x~ G x (so E8 has minimum 2).
* An O_K-structure, K = Q(sqrt(-d)), O_K = Z[eta], eta^2 - t*eta + N = 0
  (t = 1 for odd discriminant, t = 0 for even; d = 4N - t^2 = |disc K|),
  is an integer matrix J acting on coordinate columns with
        J^2 - t J + N = 0,     J~ G = G (t - J)      (adjoint = conjugate).
  Every function takes the pair (N, t) as its `K` argument, e.g. K = (3, 1)
  for Q(sqrt(-11)), (2, 1) for Q(sqrt(-7)), (1, 1) for Q(sqrt(-3)), (2, 0) for
  Q(sqrt(-2)), (5, 1) for Q(sqrt(-19)).
  The Hermitian form h (linear in the first variable) is recovered from
  G(x,y) = Tr_{K/Q} h(x,y);  C-N's "(P, Tr o h)" convention.
"""
import os, math, random
import cypari2

pari = cypari2.Pari()
pari.allocatemem(2 * 10**9, silent=True)
_GP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "gp", "ideal_lattices.gp")
pari('read("%s")' % _GP)
_STRAT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "fplll_default_strategies.json")


# ----------------------------------------------------------------------------
# basic invariants
# ----------------------------------------------------------------------------
def is_even(G):
    n = len(G)
    return all(int(G[i][i]) % 2 == 0 for i in range(n)) and pari.denominator(G) == 1

def level(G):
    """Smallest l with l*G^-1 integral and even."""
    Gi = G ** -1
    l = int(pari.denominator(Gi))
    if any(int(l * Gi[i][i]) % 2 for i in range(len(G))):
        l *= 2
    return l

def _K(K):
    return (K, 1) if isinstance(K, int) else tuple(K)

def disc(K):
    N, t = _K(K)
    return 4 * N - t * t

def is_structure(G, J, K):
    N, t = _K(K)
    one = pari.matid(len(G))
    return (J * J - t * J + N * one == 0) and (pari.mattranspose(J) * G == G * (t * one - J))

def lll(G, J=None):
    U = pari.qflllgram(G, 1)
    G2 = pari.mattranspose(U) * G * U
    if J is None:
        return G2
    return G2, U ** -1 * J * U

def count_short(G, B):
    """Number of nonzero vectors of norm <= B (exact, Fincke-Pohst)."""
    return int(pari.qfminim(G, B, 0)[0])

def minimum_exact(G):
    """Exact minimum; only for small rank (<= ~40)."""
    r = pari.qfminim(G)
    return int(r[1])


# ----------------------------------------------------------------------------
# Hermitian lattices over O_K = Z[eta], eta^2 - eta + N = 0
# ----------------------------------------------------------------------------
def herm_to_z(H, K):
    """H: m x m Hermitian matrix, entries (u, v) meaning u + v*eta (rationals).
    Returns (G, J), G = Tr o h, on the Z-basis e_1, eta e_1, e_2, eta e_2, ..."""
    N, t = _K(K)
    m = len(H)
    w = pari("Mod(y, y^2 - %d*y + %d)" % (t, N))
    wb = t - w
    one = pari(1) + 0 * w
    el = lambda p: pari(p[0]) + pari(p[1]) * w
    ent = []
    for i in range(m):
        for a in (one, w):
            for j in range(m):
                for b in (one, wb):
                    ent.append(pari.trace(a * el(H[i][j]) * b))
    G = pari.matrix(2 * m, 2 * m, ent)
    J = pari.matrix(2 * m, 2 * m)
    for i in range(m):
        J[2 * i + 1, 2 * i] = 1
        J[2 * i, 2 * i + 1] = -N
        J[2 * i + 1, 2 * i + 1] = t
    return G, J

def standard(m, K):
    """O_K^m with the standard form (Hermitian unimodular)."""
    H = [[(1 if i == j else 0, 0) for j in range(m)] for i in range(m)]
    return herm_to_z(H, K)

def _kron(A, B):
    return pari("kron")(A, B)

def tensor_OK(G1, J1, G2, J2, K):
    """Hermitian tensor product over O_K of (G1,J1) and (G2,J2).

    V1 (x)_Q V2 = W (+) W' with W the part where J1(x)1 = 1(x)J2; the lattice
    L1 (x)_{O_K} L2 is the projection of L1 (x)_Z L2 to W, and G1 (x) G2
    restricted to W is exactly Tr(h1 h2).  Returns (G, J) LLL-reduced."""
    N, t = _K(K)
    d = 4 * N - t * t
    n1, n2 = len(G1), len(G2)
    I1, I2 = pari.matid(n1), pari.matid(n2)
    # d * projection:  -(J1(x)1 + 1(x)J2 - t) * ((2 J1 - t) (x) 1),  (2J-t)^2 = -d
    P = -(_kron(J1, I2) + _kron(I1, J2) - t * pari.matid(n1 * n2)) * _kron(2 * J1 - t * I1, I2)
    B = pari.mathnf(P)                                   # Z-basis of d * pi(L1 (x) L2)
    assert len(B) * 2 == n1 * n2, "projection has wrong rank"
    G = pari.mattranspose(B) * _kron(G1, G2) * B / d**2
    Jt = _kron(J1, I2)
    J = pari.matinverseimage(B, Jt * B)
    return lll(G, J)


# ----------------------------------------------------------------------------
# short-vector search ("kill test")
# ----------------------------------------------------------------------------
def bkz_shortest(G, block=40, loops=8, target=None, verbose=False, seed=0):
    """Heuristic search for short vectors of the lattice with Gram G.

    Returns (norm, coefficient vector) of the shortest vector found; the norm is
    recomputed EXACTLY from G.  This is an upper bound for min(L): finding a
    vector below the extremal minimum kills a candidate rigorously; finding
    none proves nothing (see PLAN.md, 'certification')."""
    import numpy as np
    from fpylll import IntegerMatrix, GSO, LLL, BKZ
    from fpylll.algorithms.bkz2 import BKZReduction
    G = lll(G)
    n = len(G)
    Gf = np.array([[float(G[i, j]) for j in range(n)] for i in range(n)])
    R = np.linalg.cholesky(Gf).T                         # G = R^T R, columns = basis
    scale = 2.0 ** 40
    A = IntegerMatrix.from_matrix([[int(round(scale * R[i, j])) for i in range(n)] for j in range(n)])
    U = IntegerMatrix.identity(n)
    M = GSO.Mat(A, U=U, float_type="long double" if n < 160 else "mpfr")
    LLL.Reduction(M)()
    best = (None, None)
    def scan():
        nonlocal best
        Gi = [[int(G[i, j]) for j in range(n)] for i in range(n)]
        for r in range(n):
            x = [U[r, j] for j in range(n)]
            y = [sum(Gi[i][j] * x[j] for j in range(n)) for i in range(n)]
            nr = sum(x[i] * y[i] for i in range(n))
            if nr > 0 and (best[0] is None or nr < best[0]):
                best = (nr, x)
    scan()
    if target is not None and best[0] <= target:
        return best
    bkz = BKZReduction(M)
    betas = sorted(set([b for b in (20, 30, 40, 50, 60, 70, 80) if b < block] + [min(block, len(G))]))
    for b in betas:
        if os.path.exists(_STRAT):
            par = BKZ.Param(b, strategies=_STRAT, max_loops=loops,
                            flags=BKZ.MAX_LOOPS | BKZ.AUTO_ABORT | BKZ.GH_BND)
        else:                                            # no pruning: keep blocks small
            par = BKZ.Param(min(b, 30), max_loops=loops, flags=BKZ.MAX_LOOPS | BKZ.AUTO_ABORT)
        bkz(par)
        scan()
        if verbose:
            print("   BKZ-%d: shortest found %d" % (b, best[0]), flush=True)
        if target is not None and best[0] <= target:
            break
    return best


# ----------------------------------------------------------------------------
# Kneser neighbours of Hermitian lattices at a split prime of degree 1
# ----------------------------------------------------------------------------
def herm_coeffs(G, J, K, x):
    """h(x, .) = u + v*eta as rational row functionals (u, v)."""
    N, t = _K(K)
    d = 4 * N - t * t
    a = pari.mattranspose(x) * G            # b(x, y)       = Tr h
    c = a * J                               # b(x, eta y)   = Tr(conj(eta) h)
    return (2 * N * a - t * c) / d, (2 * c - t * a) / d

def split_primes(K, bound=50):
    """[(p, r)] with p prime, r a root of X^2 - tX + N mod p, p not dividing disc:
    the degree-one primes  P = (p, eta - r)."""
    N, t = _K(K)
    out = []
    for p in range(2, bound):
        if all(p % q for q in range(2, int(p ** .5) + 1)) and disc(K) % p:
            out += [(p, r) for r in range(p) if (r * r - t * r + N) % p == 0]
    return out

def neighbour(G, J, K, x, pr=None, scale=1):
    """P-neighbour of the O_K-lattice (G,J), P = (p, eta - r) a split prime:
         L' = { y in L : h(x,y) in P }  +  P^-1 conj(P) x .
    Index p in and out, same Hermitian genus.  `scale`: integer s such that
    s*h is O_K-valued (1 for Hermitian unimodular lattices).
    Returns (G', J') LLL-reduced, or None if h(x, L) lies in P."""
    N, t = _K(K)
    p, r = pr if pr else split_primes(K)[0]
    n = len(G)
    u, v = herm_coeffs(G, J, K, x)
    f = scale * (u + r * v)
    assert pari.denominator(f) == 1, "scale*h is not O_K-valued"
    if all(int(f[0, j]) % p == 0 for j in range(n)):
        return None
    Kr = pari.matkerint(pari.matconcat([f, pari.matrix(1, 1, [p])]))
    Kr = pari.matrix(n, len(Kr), [Kr[i, j] for i in range(n) for j in range(len(Kr))])
    C = (t - r) * pari.matid(n) - J         # conj(eta) - r
    xp = C * C * x                          # p * generator of P^-1 conj(P) x
    B = pari.mathnf(pari.matconcat([p * Kr, xp, J * xp])) / p
    assert len(B) == n
    return lll(pari.mattranspose(B) * G * B, B ** -1 * J * B)

def random_vec(n, bound=2):
    while True:
        v = [random.randint(-bound, bound) for _ in range(n)]
        if any(v):
            return pari.matrix(n, 1, v)

def dsum(*Ms):
    """Orthogonal direct sum of Gram matrices (or block-diagonal of operators)."""
    return pari("(v)->matconcat(matdiagonal(v))")(list(Ms))
