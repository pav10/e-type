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
