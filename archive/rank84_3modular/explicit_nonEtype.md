# An explicit non-E-type lattice, and how small one could be

## Part I. The tower construction (explicit, rank ≈ 7040)

**Setup.** For lattices $L, M$, every $z \in L\otimes M$ of rank $r$ (minimal number
of pure-tensor summands) with $r \ge 2$ is *non-split*. $L$ is of **E-type** if
$\min(L\otimes M) = \min(L)\min(M)$ for all $M$; equivalently no non-split
minimal vectors ever occur.

**Lemma (the witness).** For any lattice $L$ of rank $n$ with dual bases
$\{e_i\}$ of $L$ and $\{e_i^\#\}$ of $L^\#$, put
$$z = \sum_{i=1}^{n} e_i \otimes e_i^\# \ \in\ L\otimes L^\#.$$
Then $z$ has rank $n$ and
$$(z,z) = \sum_{i,j}(e_i,e_j)(e_i^\#,e_j^\#) = \operatorname{tr}(GG^{-1}) = n .$$
Hence if $\min(L)\cdot\min(L^\#) > n$, then $z$ is strictly shorter than every
split vector, so $L\otimes L^\#$ has non-split minimal vectors and **$L$ is not
of E-type.** (Mechanism: Coulangeon, *Tensor products of Hermitian lattices*,
Prop. 3.2 equality case; the criterion $\gamma'(L)^2>n$ is the Bergé–Martinet
form of Kitaoka's inequality.)

**Theorem.** Let $K$ be totally imaginary of degree $n$, $\mathcal{O}_K$ its ring
of integers, $\mathfrak{d}$ the different, and $\mathrm{rd}(K)=|d_K|^{1/n}$. Put
$$L = \big(\mathcal{O}_K,\ \operatorname{Tr}_{K/\mathbb{Q}}(x\bar y)\big),
\qquad L^\# = \big(\mathfrak{d}^{-1},\ \text{same form}\big).$$
If $\;n > \mathrm{rd}(K)^2\;$ then $L$ is not of E-type.

*Proof.* For $0\ne x\in\mathcal{O}_K$, $\operatorname{Tr}(x\bar x)=\sum_\sigma|\sigma(x)|^2
\ge n|N(x)|^{2/n}\ge n$ by AM–GM, with equality at $x=1$; so $\min(L)=n$ exactly.
For $0\ne y\in\mathfrak{d}^{-1}$ the ideal $(y)\subseteq\mathfrak{d}^{-1}$ gives
$|N(y)|\ge 1/|d_K|$, so $\operatorname{Tr}(y\bar y)\ge n|d_K|^{-2/n}=n/\mathrm{rd}^2$.
Therefore $\min(L)\min(L^\#)\ \ge\ n^2/\mathrm{rd}^2\ >\ n$ iff $n>\mathrm{rd}^2$. $\square$

Two remarks. The bound is **scale- and ideal-independent**: replacing
$\mathcal{O}_K$ by any fractional ideal $I$ gives the same product $n^2/\mathrm{rd}^2$.
And **no shortest-vector computation is required** — the inequality is automatic,
which is what makes the construction usable at huge rank.

Sharp form: the condition fails only if some $y\in\mathfrak{d}^{-1}$ has
$\sum_\sigma|\sigma(y)|^2\le 1$, which by AM–GM forces $|d_K|^2\ge n^n$, i.e.
$\mathrm{rd}^2\ge n$. So $n>\mathrm{rd}^2$ is exactly the right condition for this argument.

**Instantiation.** Root discriminant is *constant* along an unramified tower
while the degree tends to infinity, so any infinite tower supplies fields of
arbitrarily large degree at fixed $\mathrm{rd}$. Records for totally complex fields:
Martinet (1978) $\mathrm{rd}\approx 92.4$; **Hajir–Maire (2001), $\mathrm{rd} < 83.9$.**

$$\boxed{\ n > 83.9^2 \approx 7040\ }$$

Any field of degree $\ge 7040$ in the Hajir–Maire tower gives an explicit
non-E-type lattice.

**Relation to Steinberg.** Steinberg's unpublished theorem gives unimodular
$L,M$ with $\min(L\otimes M)<\min L\min M$ in every dimension $\ge 292$, non-constructively
— it is the Minkowski–Hlawka statement $\gamma'_n \ge n/2\pi e$ (Lagarias–Lenstra–Schnorr).
The above is its **explicit form**, and the cost of explicitness is exactly
$$\left(\frac{83.9}{2\pi e}\right)^2 = \left(\frac{83.9}{17.08}\right)^2 \approx 24\times$$
in rank: $292 \times 24 \approx 7040$. Any improvement to the small-root-discriminant
record lowers the rank **quadratically**.

**Prior art — must be checked.** Every ingredient is standard: algebraic lattices
from Hilbert class field towers with bounded root discriminant are routine in
coding theory (Litsyn–Tsfasman; Vehkalahti–Luzzi), $\gamma'_n\ge n/2\pi e$ is
Lagarias–Lenstra–Schnorr, and the criterion is Kitaoka/Coulangeon. The
*combination* may be folklore. Check Litsyn–Tsfasman's constructive
sphere-packing papers and Bergé–Martinet (1989) before claiming anything.

---

## Part II. How small can a witness be? (corrected)

**Correction.** I earlier said the witness must be *extremal*, giving exactly four
candidates. That was wrong. "Level $\ell$ dies at dimension $D$" means the
*extremal* theta series acquires a negative coefficient — it rules out
$\min = E$, not $\min = E-2$. At levels 7, 11, 23, where extremal lattices stop
existing early, **non-extremal** lattices can still be witnesses.

For an $\ell$-modular lattice of rank $n$: $\gamma'^2 = \min^2/\ell$, so the
required minimum is $R = $ least even integer with $R^2 > \ell n$, subject to
$R \le E = 2+2\lfloor n/m_\ell\rfloor$.

Since $R \in (\sqrt{\ell n},\, \sqrt{\ell n}+2]$, the required Hermite invariant is
$$\gamma_{\text{req}} = R/\sqrt{\ell} \ \in\ \big(\sqrt n,\ \sqrt n + 2/\sqrt\ell\,\big],$$
so **higher level = tighter requirement.** Smallest excesses over $\sqrt n$:

| excess | $n$ | $\ell$ | need min | extremal | $\gamma$ required | $\sqrt n$ |
|---|---|---|---|---|---|---|
| 0.023 | **56** | 23 | 36 | 58 | 7.507 | 7.483 |
| 0.025 | **52** | 11 | 24 | 28 | 7.236 | 7.211 |
| 0.050 | 64 | 5 | 18 | 18 | 8.050 | 8.000 |
| 0.072 | 84 | 3 | 16 | 16 | 9.238 | 9.165 |
| 0.076 | 56 | 7 | 20 | 20 | 7.559 | 7.483 |

**The structural point.** Because $\gamma_{\text{req}}$ is always *just barely*
above $\sqrt n$, every witness at every rank requires a lattice essentially **at
the linear-programming bound**. This is precisely why $\gamma_r < \sqrt r$ —
certified for $r\le 48$, borderline $49$–$51$ — closes everything below $\approx 52$.

**Consequence: the LP extension and the witness hunt are the same computation.**
Pushing the certification $\gamma_r<\sqrt r$ upward *kills* candidates; the rank
at which it fails is exactly the rank at which the smallest witness becomes
possible. Extending the certification from 48 to 52 would eliminate the
$(n{=}52,\ \ell{=}11,\ \min 24)$ candidate; if instead the LP permits
$\gamma_{52} > \sqrt{52}$, that candidate is the smallest conceivable witness.

So the double-digit frontier is $n = 52$–$56$, and the deciding computation is
the LP one already in progress.
