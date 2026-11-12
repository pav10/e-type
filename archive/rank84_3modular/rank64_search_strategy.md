# Rank-64 level-5 candidate: search strategy (and a retracted lead)

## 1. CORRECTION — the Quebbemann rank-64 lead was wrong

I flagged arXiv:2106.04828 (*"A note on Quebbemann's extremal lattices of rank
64"*, Shimada) as directly relevant machinery. Having read it: **it is not.**

Those lattices are **even unimodular** of rank 64, where "extremal" means
min = 2 + 2⌊64/24⌋ = **6**, det 1, |Min| = 2 611 200. So

        γ = 6/1 = 6  <  √64 = 8.

They do **not** witness γ' ≥ √n and are irrelevant to the E-type question.
They are level 1 with minimum 6, not level 5 with minimum 18 — a different genus
entirely. Retract that lead.

What *is* transferable from the paper is methodology, not objects: enumerating
Min(L) by norm-types over a decomposition, distinguishing isomorphism classes by
the distribution of intersection patterns ⟨v,v'⟩, and the Γ-rigidifying-basis
trick for proving O(L) = Γ.

## 2. Genus enumeration is not viable — quantified

Minkowski–Siegel mass for even unimodular rank n
(formula validated exactly: n=8 gives 1/696729600 = 1/|Aut(E₈)|; n=32 gives
4.03×10⁷, the known value):

| rank | log₁₀(mass) |
|---|---|
| 32 | +7.6 |
| 48 | +121.2 |
| **64** | **+346.7** |
| 80 | +700.5 |

Shimada also quotes: the number of isomorphism classes of extremal rank-40
lattices already exceeds **8.45 × 10⁵¹**.

The rank-64 5-modular genus is comparably vast. **Kneser/Schiemann neighbour
enumeration of the genus is hopeless by hundreds of orders of magnitude.** The
suggestion I made two turns ago ("mass formula plus neighbours") is not a real
plan at this rank, and I should not have offered it without checking the size.

## 3. The only viable route: prescribed automorphisms

This is what Nebe did for dim 72, and what Jürgens' dissertation does for the
open 3-modular dim-36 problem: instead of enumerating the genus, fix a group
G ≤ Aut(L) and classify the (far fewer) G-invariant lattices.

**Constraint derived here.** Let σ ∈ Aut(L) have prime order p. The free part
has Z-rank k(p−1) and the fixed part f = 64 − k(p−1) ≥ 0, so p−1 ≤ 64:

        p ∈ {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61}

**Fixed-point-free case** (f = 0) requires (p−1) | 64, and 64's divisors are
{1,2,4,8,16,32,64}, so

        **p ∈ {2, 3, 5, 17}.**

Two of these are structurally interesting:
- **p = 5** is the level itself — σ interacts with the 5-modularity directly.
- **p = 17**: the free part is a single Z[ζ₁₇]-module of rank 4 (4 × 16 = 64).
  This is the exact analogue of the p = 13, 17 cases Nebe isolates for dim 96,
  and is the most promising place to start.

## 4. Suggested concrete programme

1. **p = 17, fixed-point-free.** Classify rank-4 Hermitian Z[ζ₁₇]-lattices whose
   trace lattice is 5-modular of rank 64; test min = 18. Small search space.
2. **p = 5.** Exploit the coincidence of automorphism order with the level;
   the σ-eigenspace decomposition is compatible with the modularity.
3. **p = 3, 2.** Larger spaces, use only if 1–2 fail.
4. For any candidate produced, verify against the certified extremal theta series
   (computed earlier): min 18, kissing number **21 281 280 000**.

## 5. Honest status

Nothing here produces a lattice. What this round establishes:
- the Quebbemann lead is retracted (wrong genus, wrong minimum);
- genus enumeration is quantitatively impossible at rank 64;
- the search must be automorphism-prescribed, and the fixed-point-free prime
  orders are pinned to {2, 3, 5, 17}, with p = 17 the natural first target.
