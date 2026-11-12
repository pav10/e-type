# The p = 17 structure for the rank-64 level-5 candidate

Setting: L a rank-64, 5-modular (det 5³²), minimum-18 lattice with a
**fixed-point-free automorphism σ of order 17**. Then K = Q(ζ₁₇) acts and L is
an O_K-module of rank 4 (4 × φ(17) = 4 × 16 = 64).

## Arithmetic of the level 5 in K = Q(ζ₁₇)

- **ord(5) mod 17 = 16**, so 5 is a *primitive root* ⟹ **5 is inert in K**,
  residue field O_K/(5) = **F₅¹⁶**.
- Frobenius at 5 has order 8 in Gal(K⁺/Q) ≅ Z/8 ⟹ **5 is inert in K⁺** too,
  with N_{K⁺/Q}(5) = 5⁸.
- disc(K) = 17¹⁵ = 2 862 423 051 509 815 793 (confirmed exactly; the Gram of
  (Z[ζ₁₇], Tr(x x̄)) is 17I − J with eigenvalues 1 and 17^{×15}).
- **h(Q(ζ₁₇)) = 1**, so O_K is a PID and **L is free: L ≅ O_K⁴**.

## The consistency check that could have failed — it passes

L is 5-modular ⟹ L*/L ≅ (Z/5)³². But L is an O_K-module, so L*/L is a module
over O_K/(5) = F₅¹⁶, forcing dim_{F₅}(L*/L) to be a **multiple of 16**:

        32 = 2 × 16  ✓

So L*/L is a 2-dimensional F₅¹⁶-vector space. This is a genuine constraint — had
the determinant exponent not been divisible by 16, p = 17 would have died here.
It doesn't.

## Determinant condition on the Hermitian Gram

With the codifferent-twisted trace form S(x,y) = Tr_{K/Q}(γ h(x,y)),
det(L_Z) = N_{K⁺/Q}(det H)². Requiring 5³² gives N_{K⁺/Q}(det H) = 5¹⁶, and
since 5 is inert in K⁺ with norm 5⁸:

> **det H = (totally positive unit) × 25.**

So the search space is: **4 × 4 Hermitian matrices over Z[ζ₁₇] with determinant
a unit times 25** — a very concrete, finite-flavoured target.

## Constraint from minimal vectors

For a minimal v (S(v,v) = 18), put a := h(v,v) ∈ O_{K⁺}, totally positive.
The rank-16 submodule O_K·v has determinant N_{K⁺/Q}(a)² (the |d_K| cancels
against N(γ) = 17⁻¹⁵) and minimum exactly 18. Two bounds:

- **AM–GM** on the 16 conjugates of γa: N(a)^{1/8} ≤ (18/16)·17^{15/16} = 16.02
  ⟹ N_{K⁺/Q}(a) ≤ 4.34 × 10⁹
- **Hermite** (γ₁₆ ≤ 3.03): 18 ≤ γ₁₆·N(a)^{1/8} ⟹ N_{K⁺/Q}(a) ≥ 1.55 × 10⁶

        1.55 × 10⁶  ≤  N_{K⁺/Q}(a)  ≤  4.34 × 10⁹

Non-empty (ratio ≈ 2800). **No contradiction — p = 17 survives.**

## The trivial form fails, by exactly 2

The untwisted case a = 1 gives the lattice (Z[ζ₁₇], Tr(x x̄)) = 17I − J, whose
minimum is **16**, attained by the 34 vectors ±ζ^i. Required: ≥ 18.

So the naive Hermitian form ⟨1⟩ misses by exactly 2, and the form **must** be
genuinely twisted by a totally positive a ∈ O_{K⁺} with norm in the window above.

## Full circle: the problem reduces to degree 8

        Tr_{K/Q}(γ a) = Tr_{K⁺/Q}( a · Tr_{K/K⁺}(γ) ),   K⁺ = Q(ζ₁₇)⁺, degree 8.

So the minimum-18 condition becomes a statement about **totally positive elements
of a codifferent-type ideal in a degree-8 totally real field with prescribed
trace** — exactly the species of question this whole thread opened with (there in
degree 96), but now in degree 8, where it is genuinely computable in PARI/Magma.

## Status

Everything checked so far is *consistent*: the splitting of 5, the F₅¹⁶-module
structure of L*/L, the class number, the determinant condition, and the
norm window. None of it produces a lattice, and none of it rules one out.

**Concrete next computation:** in K⁺ = Q(ζ₁₇)⁺ (degree 8), enumerate totally
positive a ∈ O_{K⁺} with 1.55×10⁶ ≤ N(a) ≤ 4.34×10⁹ and Tr_{K⁺/Q}(a·β) = 18
where β = Tr_{K/K⁺}(γ); then test whether any 4×4 Hermitian H over Z[ζ₁₇] with
det H = unit×25 and diagonal entries among those a has trace-form minimum 18.
