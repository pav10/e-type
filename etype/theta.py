"""Exact extremal theta series for l-modular lattices, l in {1,2,3,5,7}.

Quebbemann: theta of an even l-modular lattice of rank n (weight k = n/2) lies in
C[theta_0, Delta_l], with
    theta_0 = theta series of the smallest l-modular lattice (rank r0, weight r0/2)
    Delta_l = (eta(z) eta(l z))^(24/(l+1)),  weight k1 = 24/(l+1),  = q + ...
The Fricke eigenvalue of a weight-k form in this ring is (-i)^k -- NOT +1; e.g.
Delta_3 = eta(z)^6 eta(3z)^6 has eigenvalue -1 = (-i)^6 and IS admissible (K12:
theta = theta_A2^6 - 36 Delta_3... is checked below).
All arithmetic is exact (PARI integer power series).  Exponent = norm/2.
"""
from .core import pari

_BASE = {  # l: (Gram of smallest l-modular lattice, weight of Delta_l)
    1: ("[2,-1,0,0,0,0,0,0;-1,2,-1,0,0,0,0,0;0,-1,2,-1,0,0,0,-1;0,0,-1,2,-1,0,0,0;"
        "0,0,0,-1,2,-1,0,0;0,0,0,0,-1,2,-1,0;0,0,0,0,0,-1,2,0;0,0,-1,0,0,0,0,2]", 12),   # E8
    2: ("[2,-1,0,0;-1,2,-1,-1;0,-1,2,0;0,-1,0,2]", 8),                                   # D4
    3: ("[2,1;1,2]", 6),                                                                 # A2
    5: ("[2,-1,-1,-1;-1,4,-1,-1;-1,-1,4,-1;-1,-1,-1,4]", 4),                                         # det 25
    7: ("[2,1;1,4]", 3),
}

def _theta0(l, T):
    G = pari(_BASE[l][0])
    v = pari.qfrep(G, 2 * T - 2, 1)                 # v[i] = #{+-x : q(x) = 2(i+1)} (pairs)
    return [1] + [2 * int(v[i]) for i in range(T - 1)]

def _mul(a, b, T):
    c = [0] * T
    for i, x in enumerate(a):
        if x:
            for j in range(T - i):
                c[i + j] += x * b[j]
    return c

def _pow(a, e, T):
    r = [1] + [0] * (T - 1)
    for _ in range(e):
        r = _mul(r, a, T)
    return r

def _delta(l, T):
    e = 24 // (l + 1)
    s = pari("(l,e,T)->my(q='q+O('q^T)); Vec(q*(eta(q)*eta(q^l))^e)")(l, e, T + 1)
    v = [0] + [int(x) for x in s]
    return (v + [0] * T)[:T]

def extremal_theta(n, l, T=None):
    """Coefficients a_0..a_{T-1} of the unique form 1 + O(q^(1+J)) in the weight
    n/2 part of C[theta_0, Delta_l]; J = floor(n / (2*r0*k1/ ... )) as in Quebbemann.
    Returns (coeffs, extremal_minimum)."""
    G = pari(_BASE[l][0]); r0 = len(G); k0 = r0 // 2; k1 = _BASE[l][1]
    k = n // 2
    assert n % r0 == 0 or l in (3, 7), "rank not allowed at this level"
    J = k // k1
    while (k - k1 * J) % k0:
        J -= 1
    T = T or J + 4
    th, de = _theta0(l, T), _delta(l, T)
    basis = [_mul(_pow(th, (k - k1 * j) // k0, T), _pow(de, j, T), T) for j in range(J + 1)]
    f = basis[0][:]
    for j in range(1, J + 1):                       # kill q^1 .. q^J (triangular)
        c = -f[j]
        f = [x + c * y for x, y in zip(f, basis[j])]
    assert f[0] == 1 and all(x == 0 for x in f[1:J + 1])
    return f, 2 * J + 2

def target(n, l):
    f, mn = extremal_theta(n, l)
    return {"n": n, "l": l, "min": mn, "kissing": f[mn // 2], "next": f[mn // 2 + 1],
            "nonneg_checked_to": len(f) - 1, "all_nonneg": all(x >= 0 for x in f)}

def first_negative(n, l, T=60):
    f, mn = extremal_theta(n, l, T)
    return next((i for i, x in enumerate(f) if x < 0), None)
