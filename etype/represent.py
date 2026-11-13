"""Representation test: does the Hermitian lattice P (Z-data G, J) contain a
sublattice isometric to the Hermitian lattice M = (O_K^r, X)?

An O_K-linear isometric embedding F: M -> P is a Z-linear map with
F~ G_P F = G_M and F J_M = J_P F.  We embed a frame v_1..v_r of short vectors
of M (O_K-independent) by backtracking through the short vectors of P, with the
Hermitian inner products prescribed, and then check that F is integral on M.
"""
import math
import numpy as np
from .core import pari, _K, disc
from .search import short_vectors


def _np(M):
    n = len(M)
    return np.array([[int(M[i, j]) for j in range(n)] for i in range(n)], dtype=np.int64)


class Host:
    """Short vectors of P (both signs) and fast Hermitian products against them."""
    def __init__(self, G, J, K, max_norm=4):
        self.K = _K(K); self.d = disc(K)
        self.G, self.J = _np(G), _np(J)
        self.lists = {}
        V = short_vectors(G, max_norm)
        nr = np.einsum("ik,ij,jk->k", V, self.G, V)
        for h2 in np.unique(nr):
            W = V[:, nr == h2]
            self.lists[int(h2) // 2] = np.concatenate([W, -W], axis=1)
        self.f32 = {h: L.astype(np.float32) for h, L in self.lists.items()}
        self.GJ = self.G @ self.J

    def prods(self, y, h):
        """d*(u, v) of h(y, w) for all w in the list of Hermitian norm h."""
        N, t = self.K
        V = self.f32[h]
        a = np.rint((self.G @ y).astype(np.float32) @ V).astype(np.int32)
        c = np.rint((self.GJ.T @ y).astype(np.float32) @ V).astype(np.int32)
        return 2 * N * a - t * c, 2 * c - t * a


def frame(GM, JM, K, max_norm=4):
    """O_K-independent short vectors of M (norm <= max_norm), or None."""
    N, t = _K(K)
    G, J = _np(GM), _np(JM)
    V = short_vectors(GM, max_norm)
    V = np.concatenate([V, -V], axis=1)
    nr = np.einsum("ik,ij,jk->k", V, G, V)
    r = len(G) // 2
    chosen = []
    for idx in np.argsort(nr, kind="stable"):
        v = V[:, idx]
        cand = chosen + [v]
        B = np.array([w for u in cand for w in (u, J @ u)], dtype=float)      # Z-span of O_K-span
        if np.linalg.matrix_rank(B, tol=1e-8) == 2 * len(cand):
            chosen = cand
            if len(chosen) == r:
                return chosen, [int(v @ G @ v) // 2 for v in chosen]
    return None


def targets(GM, JM, K, fr):
    N, t = _K(K)
    G, J = _np(GM), _np(JM)
    r = len(fr)
    out = {}
    for i in range(r):
        for k in range(i + 1, r):
            a = int(fr[i] @ G @ fr[k]); c = int(fr[i] @ G @ (J @ fr[k]))
            out[(i, k)] = (2 * N * a - t * c, 2 * c - t * a)
    return out


def embed(host, GM, JM, K, max_first=None, start=0, step=1):
    """Search an embedding of M into the host.  Returns ("found", [y_1..y_r]),
    ("none", None) after an exhaustive search, or ("skip", reason)."""
    fr = frame(GM, JM, K, max_norm=2 * max(host.lists))
    if fr is None:
        return "skip", "no frame of short vectors"
    vs, hs = fr
    tg = targets(GM, JM, K, vs)
    r = len(vs)
    L1 = host.lists[hs[0]]
    n1 = L1.shape[1] // 2                              # first vector up to sign
    G, J = _np(GM), _np(JM)
    Vm = np.array([w for v in vs for w in (v, J @ v)], dtype=float).T
    for i1 in range(start, n1 if max_first is None else min(n1, max_first), step):
        ys = [L1[:, i1]]
        cands = {}
        for k in range(1, r):
            du, dv = host.prods(ys[0], hs[k])
            cands[k] = np.nonzero((du == tg[(0, k)][0]) & (dv == tg[(0, k)][1]))[0]
            if len(cands[k]) == 0:
                break
        else:
            sol = _extend(host, ys, cands, hs, tg, 1, r)
            if sol is not None:
                Ym = np.array([w for y in sol for w in (y, host.J @ y)], dtype=float).T
                F = Ym @ np.linalg.inv(Vm)
                if np.abs(F - np.rint(F)).max() < 1e-6:
                    return "found", sol
    return "none", None


def _extend(host, ys, cands, hs, tg, k, r):
    if k == r:
        return ys
    N, t = host.K
    Lk = host.lists[hs[k]]
    for idx in cands[k]:
        y = Lk[:, idx]
        new = {}
        ok = True
        for m in range(k + 1, r):
            W = host.lists[hs[m]][:, cands[m]]
            a = (host.G @ y) @ W; c = (host.GJ.T @ y) @ W
            du, dv = 2 * N * a - t * c, 2 * c - t * a
            new[m] = cands[m][(du == tg[(k, m)][0]) & (dv == tg[(k, m)][1])]
            if len(new[m]) == 0:
                ok = False; break
        if ok:
            sol = _extend(host, ys + [y], {**cands, **new}, hs, tg, k + 1, r)
            if sol is not None:
                return sol
    return None


# ----------------------------------------------------------------------------
# Hermitian automorphisms and orbits (PARI qfauto / qfisom)
# ----------------------------------------------------------------------------
def herm_forms(G, J, K):
    """[G, S] with S = G(2J - t) antisymmetric: a Z-map preserving both is an
    O_K-LINEAR isometry (anti-linear ones send S to -S)."""
    N, t = _K(K)
    return [G, G * (2 * J - t * pari.matid(len(G)))]

def herm_aut(G, J, K):
    A = pari.qfauto(herm_forms(G, J, K))
    return int(A[0]), [_np(g) for g in A[1]]

def herm_isom(G1, J1, G2, J2, K):
    """True iff the two Hermitian lattices are O_K-linearly isometric."""
    r = pari.qfisom(herm_forms(G1, J1, K), herm_forms(G2, J2, K))
    return r != 0

def orbit_reps(V, gens):
    """V: n x M integer array closed under the group <gens>; returns indices of
    orbit representatives and orbit sizes."""
    key = {tuple(V[:, i]): i for i in range(V.shape[1])}
    imgs = [np.array([key[tuple(c)] for c in (g @ V).T]) for g in gens]
    seen = np.zeros(V.shape[1], dtype=bool)
    reps, sizes = [], []
    for i in range(V.shape[1]):
        if seen[i]:
            continue
        orb = [i]; seen[i] = True; k = 0
        while k < len(orb):
            for im in imgs:
                j = im[orb[k]]
                if not seen[j]:
                    seen[j] = True; orb.append(j)
            k += 1
        reps.append(i); sizes.append(len(orb))
    return reps, sizes


# ----------------------------------------------------------------------------
# Exhaustive embedding test for frames of minimal vectors (all norms 2)
# ----------------------------------------------------------------------------
class MinHost:
    """Minimal vectors of P (Hermitian norm hmin), their Aut-orbit representatives,
    and Hermitian products coded as integers  code = du * 4096 + dv."""
    def __init__(self, G, J, K, hnorm=2, gens=None):
        self.K = _K(K); self.G = _np(G); self.J = _np(J)
        V = short_vectors(G, 2 * hnorm)
        nr = np.einsum("ik,ij,jk->k", V, self.G, V)
        V = V[:, nr == 2 * hnorm]
        self.V = np.concatenate([V, -V], axis=1)
        self.Vf = self.V.astype(np.float32)
        self.GV = (self.G @ self.V).astype(np.float32)               # b(., v)
        self.GJV = (self.G @ (self.J @ self.V)).astype(np.float32)   # b(., J v)
        self.reps = list(range(self.V.shape[1])) if gens is None else orbit_reps(self.V, gens)[0]

    def code(self, I, Jdx):
        """codes of h(v_i, v_j) for i in I, j in Jdx (index arrays)."""
        N, t = self.K
        a = np.rint(self.Vf[:, I].T @ self.GV[:, Jdx]).astype(np.int64)
        c = np.rint(self.Vf[:, I].T @ self.GJV[:, Jdx]).astype(np.int64)
        return (2 * N * a - t * c) * 4096 + (2 * c - t * a)


def _code(du, dv):
    return du * 4096 + dv


def _row(host, i1):
    cache = host.__dict__.setdefault("_rows", {})
    if i1 not in cache:
        cache[i1] = host.code(np.array([i1]), np.arange(host.V.shape[1]))[0]
    return cache[i1]


def embed_frame(host, tg, r, check=None, reps=None):
    """Search y_1..y_r among minimal vectors with codes h(y_i, y_j) = tg[(i,j)]
    (i < j), y_1 running over orbit representatives (or the subset `reps`).
    check(ys) -> bool is the final integrality test.  Returns ys or None (exhaustive)."""
    allidx = np.arange(host.V.shape[1])
    for i1 in (host.reps if reps is None else reps):
        row = _row(host, i1)
        C = [None] + [allidx[row == tg[(0, k)]] for k in range(1, r)]
        if any(len(c) == 0 for c in C[1:]):
            continue
        sol = _ext(host, [i1], C, tg, 1, r, check)
        if sol is not None:
            return sol
    return None


def _ext(host, ys, C, tg, k, r, check):
    if k == r:
        return ys if (check is None or check(ys)) else None
    Ck = C[k]
    if k + 1 < r:
        rest = [host.code(Ck, C[m]) == tg[(k, m)] for m in range(k + 1, r)]
    for a, ik in enumerate(Ck):
        if k + 1 < r:
            C2 = C[:k + 1] + [C[m][rest[m - k - 1][a]] for m in range(k + 1, r)]
            if any(len(c) == 0 for c in C2[k + 1:]):
                continue
        else:
            C2 = C
        sol = _ext(host, ys + [ik], C2, tg, k + 1, r, check)
        if sol is not None:
            return sol
    return None


def embeds(host, GM, JM, K):
    """Is the Hermitian lattice (GM, JM) (spanned by a frame of minimal vectors
    of norm 2*hmin) an O_K-sublattice of the host?  Exhaustive.
    Returns (True, ys) / (False, None) / (None, reason) if the frame is not minimal."""
    N, t = _K(K)
    fr = frame(GM, JM, K, max_norm=4)
    if fr is None or any(h != 2 for h in fr[1]):
        return None, "frame not made of norm-4 vectors"
    vs = fr[0]; r = len(vs)
    G, J = _np(GM), _np(JM)
    tg = {}
    for i in range(r):
        for k in range(i + 1, r):
            a = int(vs[i] @ G @ vs[k]); c = int(vs[i] @ G @ (J @ vs[k]))
            tg[(i, k)] = _code(2 * N * a - t * c, 2 * c - t * a)
    Vm = np.array([w for v in vs for w in (v, J @ v)], dtype=float).T
    Vmi = np.linalg.inv(Vm)
    def check(ys):
        Ym = np.array([w for i in ys for w in (host.V[:, i], host.J @ host.V[:, i])], dtype=float).T
        F = Ym @ Vmi
        return np.abs(F - np.rint(F)).max() < 1e-6
    sol = embed_frame(host, tg, r, check)
    return (sol is not None), sol


def embeds_fast(host, GM, JM, K, freq, prepare_only=False):
    """As `embeds`, but the frame is chosen greedily among all minimal vectors of M
    so that the prescribed products are as rare as possible in P.
    freq: {code: average number of minimal vectors w of P with h(y, w) = code}."""
    N, t = _K(K)
    G, J = _np(GM), _np(JM)
    r = len(G) // 2
    V = short_vectors(GM, 4)
    if V.shape[1] == 0:
        return None, "no minimal vectors"
    V = np.concatenate([V, -V], axis=1)
    nr = np.einsum("ik,ij,jk->k", V, G, V)
    V = V[:, nr == 4]
    a = V.T @ G @ V; c = V.T @ G @ (J @ V)
    C = (2 * N * a - t * c) * 4096 + (2 * c - t * a)
    F = np.vectorize(lambda x: freq.get(int(x), 0.0))(C)
    tot = float(host.V.shape[1])
    def indep(idx):
        B = np.array([w for i in idx for w in (V[:, i], J @ V[:, i])], dtype=float)
        return np.linalg.matrix_rank(B, tol=1e-8) == 2 * len(idx)
    best = None
    m = V.shape[1]
    pairs = sorted(((F[i, j], i, j) for i in range(m) for j in range(m) if i != j), key=lambda x: x[0])
    tried = 0
    for _, i, j in pairs:
            if tried >= 40:
                break
            if not indep([i, j]):
                continue
            tried += 1
            chosen = [i, j]; cost = F[i, j]; level = F[i, j]
            while len(chosen) < r:
                bestk = None
                for k in range(m):
                    if k in chosen or not indep(chosen + [k]):
                        continue
                    lv = level * F[chosen[0], k] * np.prod([F[x, k] / tot for x in chosen[1:]])
                    if bestk is None or lv < bestk[0]:
                        bestk = (lv, k)
                if bestk is None:
                    break
                level = bestk[0]; cost += level; chosen.append(bestk[1])
            if len(chosen) == r:
                B = np.array([w for i in chosen for w in (V[:, i], J @ V[:, i])], dtype=float)
                unimod = abs(abs(np.linalg.det(B)) - 1) < 1e-6
                key = cost * (1.0 if unimod else 3.0)
                if best is None or key < best[0]:
                    best = (key, chosen)
    if best is None:
        return None, "minimal vectors do not span"
    idx = best[1]; vs = [V[:, i] for i in idx]
    tg = {(p, q): int(C[idx[p], idx[q]]) for p in range(r) for q in range(p + 1, r)}
    Vm = np.array([w for v in vs for w in (v, J @ v)], dtype=float).T
    Vmi = np.linalg.inv(Vm)
    def check(ys):
        Ym = np.array([w for i in ys for w in (host.V[:, i], host.J @ host.V[:, i])], dtype=float).T
        Fm = Ym @ Vmi
        return np.abs(Fm - np.rint(Fm)).max() < 1e-6
    if prepare_only == "map":
        sol = embed_frame(host, tg, r, check)
        if sol is None:
            return None
        Ym = np.array([w for i in sol for w in (host.V[:, i], host.J @ host.V[:, i])], dtype=float).T
        return np.rint(Ym @ Vmi).astype(np.int64)
    if prepare_only == "cost":
        return best[0], idx
    if prepare_only:
        return tg, r, check
    sol = embed_frame(host, tg, r, check)
    return (sol is not None), sol


def embeds_any(host, shapes, K, freq, chunk=2000):
    """First shape (index, ys) that embeds, scanning all shapes chunk by chunk of
    starting vectors; (None, None) if none embeds (exhaustive)."""
    preps = [embeds_fast(host, G, J, K, freq, prepare_only=True) for G, J in shapes]
    for s0 in range(0, len(host.reps), chunk):
        sub = host.reps[s0:s0 + chunk]
        for k, pr in enumerate(preps):
            if pr[0] is None:
                continue
            sol = embed_frame(host, pr[0], pr[1], pr[2], reps=sub)
            if sol is not None:
                return k, sol
    return None, None


def code_freq(host, sample=200, seed=0):
    """Number of minimal vectors with each product code (max over a sample of orbit
    representatives).  Only used to ORDER the search: exactness does not depend on it."""
    allidx = np.arange(host.V.shape[1])
    cnt = {}
    reps = host.reps
    if len(reps) > sample:
        reps = list(np.random.default_rng(seed).choice(reps, sample, replace=False))
    for i in reps:
        u, n = np.unique(host.code(np.array([i]), allidx)[0], return_counts=True)
        for a, b in zip(u.tolist(), n.tolist()):
            cnt[a] = max(cnt.get(a, 0), b)          # worst case over orbits
    return cnt
