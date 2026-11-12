# Prior-art check: is the p = 17 structure already known impossible?

Verdict: **no — but it came within one step of it, and the methodology I was
reconstructing is already established. Credit is owed.**

## 1. Bayer-Fluckiger kills the rank-1 version outright

*Cyclotomic modular lattices*, J. Théor. Nombres Bordeaux 12 (2000), 273–280.

> **Theorem 1.** There exists a modular ideal lattice over Q(ζ_m) **iff** m is
> not a power of a prime p with p ≡ 1 (mod 4).

Our field is K = Q(ζ₁₇): m = 17 is a prime power and **17 ≡ 1 (mod 4)**.

> **⟹ There is NO modular ideal lattice over Q(ζ₁₇), of any level.**

An *ideal lattice* is rank 1 over O_K (Z-rank 16). So any attempt to realise the
p = 17 structure as an ideal lattice is dead on arrival. Had I proposed the rank-1
route, I would have been reproving a 25-year-old impossibility.

**Second warning from the same paper.** Every construction there (Props 1, 2, 6;
Theorems 2, 3) produces ℓ-modular lattices only for **ℓ = 1, 2, or ℓ ≡ 3 (mod 4)**.
Our level is **ℓ = 5 ≡ 1 (mod 4)** — outside the scope of the classification
entirely. Trace-type cyclotomic ideal lattices simply do not produce level 5.

## 2. The underlying obstruction is a parity condition — and rank 4 passes

Working out *why* Bayer-Fluckiger's condition bites here, in our situation:

Since ord(5) mod 17 = 16, **(5) is inert** in K, a prime ideal of norm 5¹⁶.
An O_K-linear similarity σ with multiplier 5 on a rank-m Hermitian lattice
satisfies σ̄ᵀHσ = 5H, so (det σ)(det σ̄) = (5)^m. As (det σ) and (det σ̄) are
conjugate ideals of equal norm, N(det σ) = 5^{8m}, and (det σ) = (5)^k forces
16k = 8m, i.e.

> **k = m/2 — so the O_K-rank m must be EVEN.**

| O_K-rank | needs (5)^{m/2} | status |
|---|---|---|
| 1 (ideal lattice) | (5)^{1/2} | **impossible** — matches Bayer-Fluckiger Thm 1 |
| 2 | (5)¹ | fine |
| 3 | (5)^{3/2} | impossible |
| **4 (our case)** | **(5)²** | **fine** |

So the obstruction is exactly a parity obstruction on the O_K-rank. Rank 4 is
even and survives. **We are not reproving a known impossibility.**

*Caveat:* this argument assumes the similarity is O_K-linear. A semi-linear
similarity (composed with a Galois automorphism) is not covered and could carry
a further obstruction; Bayer-Fluckiger's Theorem 1 presumably handles the general
rank-1 case. Worth checking before relying on the rank-4 conclusion.

## 3. The method is Nebe–Jürgens, not new

The "prime-order automorphism ⟹ constraints on the fixed and free parts" analysis
I was assembling is established machinery:

- **M. Jürgens, PhD thesis (TU Dortmund)** — introduces the *type* of an
  automorphism of prime order p of a lattice L: it determines the Z_p[σ]-module
  structure of Z_p ⊗ L, and is the standard tool for classifying extremal
  **modular** lattices with prescribed automorphisms.
- **G. Nebe, *Automorphisms of modular lattices*** (arXiv:1910.06587, J. Pure
  Appl. Algebra 225 (2021) 106606) — extends the extremal-unimodular methods to
  modular lattices and introduces the finer **det-type**, which "determines the
  det-type of all partial dual lattices" and thereby "excludes quite a few
  det-types left open in [Jürgens]".

My constraint (p−1 ≤ 64; fixed-point-free ⟹ (p−1) | 64 ⟹ p ∈ {2,3,5,17}) is a
crude special case of the type analysis. The proper tool is the **det-type**.

**What is NOT yet classified**, per Nebe 2021: the paper classifies extremal
ℓ-modular lattices of level 14, 15 in dim 12 and level 6 in dim 16. And Nebe's
Theorem 3.8 (dim72ideal paper) covers **even unimodular** dim 96/120/144 — where,
notably, p ∈ {13, 17} also survives for dim 96. Neither covers **level 5,
dimension 64**.

## 4. Required before any further work

1. **Check Jürgens' thesis tables for (ℓ = 5, dim 64).** He tabulates types of
   prime-order automorphisms across levels and dimensions; if dim 64 level 5
   appears, the p = 17 case may already be settled there.
2. **Redo the analysis with det-type**, not the crude type argument — Nebe 2021
   shows det-type kills cases the type analysis leaves open, so the p = 17 case
   might not survive it.
3. Only then is any construction attempt warranted.

## Summary

| claim | status |
|---|---|
| rank-1 / ideal-lattice route over Q(ζ₁₇) | **dead** (Bayer-Fluckiger Thm 1) |
| level 5 from trace-type cyclotomic ideal lattices | **not available** (ℓ ≡ 1 mod 4) |
| rank-4 Hermitian, p = 17, level 5, dim 64 | survives the parity obstruction |
| the automorphism-type method | **established** (Jürgens, Nebe) — not new |
| (ℓ=5, dim 64) already classified? | not in Nebe 2021; **check Jürgens' thesis** |
