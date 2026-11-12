#!/usr/bin/env python3
"""Order of magnitude of the mass of the rank-4 Hermitian genus over Q(zeta_17)
behind 'extremal 5-modular, dim 64, automorphism of order 17' (brief, action 1).

Standard part (Shimura / Gan-Hanke-Yu, rank m = 4 over K/F, [F:Q] = 8):
     2^(1 - 4*8) * |zeta_F(-1) L_F(0,chi) ... |  with
     prod_{j=1..4} L_F(1-j, chi^j) = L_F(0,chi) zeta_F(-1) L_F(-2,chi) zeta_F(-3),
each a product of Dirichlet L-values of conductor 17 (even chars for zeta_F, odd
for L_F(.,chi)).  Local factor at the inert prime 5 for the Jordan type
(unimodular rank 2) + (5-modular rank 2): the parahoric with reductive quotient
U_2 x U_2 over F_q, q = 5^8, has index-type factor
     q^(-4) |U_4(q)| / |U_2(q)|^2  ~  q^4 = 5^32.
The factor at 17 (ramified) is O(1) and ignored: this is an ORDER OF MAGNITUDE.
"""
import cypari2
pari = cypari2.Pari()
pari("G = znstar(17, 1)")
vals = {}
for s in (0, -1, -2, -3):
    pr = pari("(s)->my(p=1.); for(m=1,16, if(zncharisodd(G,znconreylog(G,m)) == (s%2==0), p*=lfun(Mod(m,17),s))); p")(s)
    vals[s] = abs(complex(pr))
    print("s = %2d : |product over %s characters mod 17| = %.6g" % (s, "odd" if s % 2 == 0 else "even", vals[s]))
std = 2.0 ** (1 - 32) * vals[0] * vals[-1] * vals[-2] * vals[-3]
q = 5.0 ** 8
U = lambda n: q ** (n * (n - 1) / 2) * __import__("math").prod(q ** i - (-1) ** i for i in range(1, n + 1))
loc5 = q ** -4 * U(4) / U(2) ** 2
print("standard part  ~ %.3g" % std)
print("local factor 5 ~ %.3g  (5^32 = %.3g)" % (loc5, 5.0 ** 32))
print("mass           ~ %.3g   => class number >= 2 * mass / 34 is astronomically large" % (std * loc5))
