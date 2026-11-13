"""Candidate generators and drivers built on etype.core."""
import json, os, random, time
from .core import *
from .core import _K

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")

def save(name, items):
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, name), "w") as f:
        json.dump([{k: (str(v) if not isinstance(v, (int, str, list, dict)) else v) for k, v in it.items()}
                   for it in items], f)

def load(name):
    with open(os.path.join(RES, name)) as f:
        items = json.load(f)
    for it in items:
        for k in ("G", "J"):
            if k in it:
                it[k] = pari(it[k])
    return items


def ideal_family(K, m, c="1", reduce_pol=True, maxsol=4):
    """All ideal lattices (I, alpha) over L = Q(zeta_m, eta) with
           alpha * I * conj(I) * Diff_L = (c),
    c a GP expression in `et` (= eta), e.g. "1" (unimodular), "2*et-1" for
    sqrt(-d) when t = 1 (Hermitian unimodular over O_K), "et" (N-modular).
    Needs bnfinit(L): fine up to degree ~32, slow beyond.
    Returns [(G, J)] LLL-reduced, J = multiplication by eta."""
    N, t = _K(K)
    pari("P = polcompositum(polcyclo(%d), x^2 - %d*x + %d)[1]" % (m, t, N))
    if int(pari("poldegree(P)")) != 2 * int(pari("eulerphi(%d)" % m)):
        return []
    if reduce_pol:
        pari("P = polredbest(P)")
    pari("P = subst(P, x, t); bnf = bnfinit(P, 1); et = Mod(nfroots(bnf.nf, x^2 - %d*x + %d)[1], P);" % (t, N))
    R = pari("ideal_lattices(bnf, %s, [et], %d)" % (c, maxsol))
    return [lll(R[i][2], R[i][3][0]) for i in range(len(R))]


def leech_structures(K, ms=(13, 21, 28, 36), verbose=True):
    """Leech lattices with an O_K-structure, found as ideal lattices over
    Q(zeta_m, eta), phi(m) = 12.  Even unimodular rank 24 without roots = Leech."""
    out = []
    for m in ms:
        t0 = time.time()
        fam = ideal_family(K, m)
        good = [(G, J) for (G, J) in fam if pari.matdet(G) == 1 and is_even(G) and count_short(G, 2) == 0]
        for G, J in good:
            assert is_structure(G, J, K)
        if verbose:
            print("  d=%d m=%d: %d ideal lattices, %d Leech   (%.0fs)" % (disc(K), m, len(fam), len(good), time.time() - t0), flush=True)
        out += [{"K": list(_K(K)), "m": m, "G": G, "J": J} for G, J in good]
    return out


def brick_walk(K, m, steps=2000, seed=0, nprimes=4):
    """Random Kneser walk through the genus of Hermitian unimodular O_K-lattices
    of rank m, starting at O_K^m.  Returns {(Zmin, kissing): (G, J)}."""
    random.seed(seed)
    G, J = standard(m, K)
    prs = split_primes(K, 60)[:nprimes]
    seen = {}
    for _ in range(steps):
        r = neighbour(G, J, K, random_vec(2 * m), random.choice(prs))
        if r is None:
            continue
        G, J = r
        mn = minimum_exact(G)
        key = (mn, count_short(G, mn))
        if key not in seen:
            assert is_structure(G, J, K) and pari.matdet(G) == disc(K) ** m and is_even(G)
            seen[key] = (G, J)
    return seen


def kill(G, target, block=50, verbose=False):
    """True (+ witness norm) if a vector of norm <= target is found."""
    nr, x = bkz_shortest(G, block=block, target=target, verbose=verbose)
    return nr <= target, nr


# ----------------------------------------------------------------------------
# Structured kill: rank-2 tensors  a1(x)p1 + a2(x)p2   (Coulangeon-Nebe, Prop 3.3)
# ----------------------------------------------------------------------------
def _kmul(x, y, K):
    N, t = _K(K)
    (u1, v1), (u2, v2) = x, y
    return (u1 * u2 - N * v1 * v2, u1 * v2 + u2 * v1 + t * v1 * v2)

def _kconj(x, K):
    N, t = _K(K)
    return (x[0] + t * x[1], -x[1])

def short_vectors(G, B):
    """numpy int array (n x count), one representative of each +-pair, norm <= B."""
    import numpy as np
    M = pari.qfminim(G, B, None, 2)[2]
    if len(M) == 0:
        return np.zeros((len(G), 0), dtype=np.int64)
    import re
    return np.array(re.findall(r"-?\d+", str(pari.mattranspose(M))), dtype=np.int64).reshape(len(M), len(G)).T

def pair_types(G, J, K, V, W=None, chunk=400, sample=None, seed=0):
    """Hermitian inner products h(x,y) = u + v*eta for x in V (optionally a random
    sample of `sample` columns), y in W (default V).
    Returns {(h(x,x), h(y,y), u, v): count}; h(x,x) = norm/2."""
    import numpy as np
    from fractions import Fraction
    N, t = _K(K); d = disc(K); n = len(G)
    Gn = np.array([[int(G[i, j]) for j in range(n)] for i in range(n)], dtype=np.int64)
    Jn = np.array([[int(J[i, j]) for j in range(n)] for i in range(n)], dtype=np.int64)
    W = V if W is None else W
    if sample and sample < V.shape[1]:
        V = V[:, np.random.default_rng(seed).choice(V.shape[1], sample, replace=False)]
    hV = np.einsum("ik,ij,jk->k", V, Gn, V) // 2
    hW = np.einsum("ik,ij,jk->k", W, Gn, W) // 2
    res = {}
    for hv in np.unique(hV):
        Vg = V[:, hV == hv]
        for hw in np.unique(hW):
            Wg = W[:, hW == hw]
            GW = (Gn @ Wg).astype(np.float64); GJW = (Gn @ (Jn @ Wg)).astype(np.float64)
            B = int(4 * N * (hv * hw * d) ** 0.5) + 4 * N + 8      # |d*u|, |d*v| < B  (Cauchy-Schwarz)
            cnt = np.zeros((2 * B + 1) ** 2, dtype=np.int64)
            for s0 in range(0, Vg.shape[1], chunk):
                X = Vg[:, s0:s0 + chunk].T.astype(np.float64)
                a, c = X @ GW, X @ GJW                               # b(x,y), b(x, J y): exact in float64
                u = (2 * N * a - t * c).astype(np.int64); v = (2 * c - t * a).astype(np.int64)
                cnt += np.bincount(((u + B) * (2 * B + 1) + (v + B)).ravel(), minlength=cnt.size)
            for k in np.nonzero(cnt)[0].tolist():
                res[(int(hv), int(hw), Fraction(k // (2 * B + 1) - B, d), Fraction(k % (2 * B + 1) - B, d))] = int(cnt[k])
    return res

def rank2_kill(A, JA, P, JP, K, target, normA=None, normP=None, sample=None, verbose=True):
    """Search for vectors of norm <= target among a1(x)p1 + a2(x)p2, where
    <a1,a2> and <p1,p2> are generated by short vectors (norm <= normA / normP;
    default: minimal vectors).  For every pair of section types the minimum of
    the rank-4 O_K-lattice A2 (x) P2 is computed EXACTLY.  Returns the list of
    killing configurations (empty = this family contains no short vector)."""
    mA, mP = minimum_exact(A), minimum_exact(P)
    VA = short_vectors(A, normA or mA); VP = short_vectors(P, normP or mP)
    tA = pair_types(A, JA, K, VA); tP = pair_types(P, JP, K, VP, sample=sample)
    if verbose:
        print("   %d / %d short vectors, %d / %d section types" % (VA.shape[1], VP.shape[1], len(tA), len(tP)), flush=True)
    kills, best, ntest = [], None, 0
    for (h1, h2, u, v) in tA:
        gA = [[(h1, 0), (u, v)], [_kconj((u, v), K), (h2, 0)]]
        dA = h1 * h2 - _kmul((u, v), _kconj((u, v), K), K)[0]
        if dA <= 0:
            continue
        for (k1, k2, x, y) in tP:
            gP = [[(k1, 0), (x, y)], [_kconj((x, y), K), (k2, 0)]]
            dP = k1 * k2 - _kmul((x, y), _kconj((x, y), K), K)[0]
            if dP <= 0 or 16 * dA * dP > target * target:      # C-N (3.6): norm >= 4 sqrt(dA dP)
                continue
            ntest += 1
            H = [[_kmul(gA[i][j], gP[k][l], K) for j in range(2) for l in range(2)] for i in range(2) for k in range(2)]
            Gz, _ = herm_to_z(H, K)
            Gz = Gz * pari.denominator(Gz)
            sc = int(pari.denominator(herm_to_z(H, K)[0]))
            mn = minimum_exact(lll(Gz)) / sc
            if best is None or mn < best:
                best = mn
            if mn <= target:
                kills.append(((h1, h2, u, v), (k1, k2, x, y), mn))
    if verbose:
        print("   %d section pairs pass the determinant bound; smallest norm found: %s" % (ntest, best), flush=True)
    return kills, best
