# (ℓ=5, n=64) theta positivity: **PASSED**

The candidate survives. Rank 64 at level 5 is not killed by the modular-form
obstruction, and is therefore a live target 20 dimensions below the rank-84 one.

## Generators — both validated, not guessed

**Weight 2:** the Eisenstein series for Γ₀(5)+,
E = 1 + 6∑σ₁(n)qⁿ − 30∑σ₁(n)q^{5n} = 1 + 6q + 18q² + 24q³ + 42q⁴ + 6q⁵ + …

*Validation:* I searched for the rank-4 even 5-modular lattice directly
(even, positive definite, det 25 = 5^{4/2}, min 2, 5G⁻¹ integral and even) and
found it with Gram

        [ 2 -1 -1 -1 ]
        [-1  2  0  1 ]
        [-1  0  4  2 ]
        [-1  1  2  4 ]

Its theta series (exponent (x,x)/2, box radius 7) is
`[1, 6, 18, 24, 42, 6, 72, 48, 90, 78, 18, 72, 168, 84]`
— **identical to E through q¹³**. So E is a genuine 5-modular theta series,
not an assumed generator.

**Weight 4:** g = η(z)⁴η(5z)⁴ = q − 4q² + 2q³ + 8q⁴ − 5q⁵ − 8q⁶ + 6q⁷ + 0q⁸ − 23q⁹ + …
the level-5 weight-4 newform, reproduced exactly by the eta product.

**Ring check:** basis Eᵃgᵇ with a + 2b = n/4 gives dim = ⌊n/8⌋+1 and hence
extremal minimum 2 + 2⌊n/8⌋ — reproducing the Quebbemann bound at *every* rank,
and matching Nebe's ExtLat table wherever extremal 5-modular lattices are known
(dims 4, 8, 12, 16, 20, 24, 28 → minima 2, 4, 4, 6, 6, 8, 8).

## Result at rank 64

Space of weight-32 forms has dimension 9; the extremal form kills q¹…q⁸, so
**minimum = 18**, exactly as predicted.

        Θ = 1 + 21 281 280 000 q⁹ + 117 273 692 160 q¹⁰
              + 6 701 253 120 000 q¹¹ + 71 579 948 400 000 q¹² + …

**All coefficients are non-negative integers out to q³¹. No negativity.**
Kissing number of a hypothetical extremal 5-modular rank-64 lattice:
**21 281 280 000** vectors of norm 18.

## Where level 5 *does* die

Scanning n = 4, 8, …, 132:

| rank | status |
|---|---|
| 4 – 84 | all coefficients ≥ 0 |
| **88** | **first negative coefficient** (at q¹³) — no extremal lattice exists |
| 96, 104, 112 | negative |
| 120 and beyond | negative throughout |

So ℓ=5 survives comfortably past 64 and dies at 88. (This also kills the
ℓ=5 entries at ranks 88 and 96 from the earlier candidate table.)

## Status of the candidate

| test | result |
|---|---|
| parity (4 \| n) | pass |
| LP exclusion (γ_r < √r certified to 48) | pass, 64 > 51 |
| level persistence | pass — ℓ=5 dies at 88, not before |
| **theta positivity** | **pass** |
| γ' > √n | 18/√5 = 8.0498 > 8.0000 |
| existence | **OPEN** |

## What this does and does not mean

- It does **not** produce the lattice. Positivity is necessary, not sufficient —
  exactly as for ranks 84 and 96.
- γ' ≥ √n is a **necessary** condition for a non-split minimal vector
  (Kitaoka's r ≤ γ_r²/ρ), not sufficient. A rank-64 extremal 5-modular lattice
  would remove the obstruction; exhibiting an actual non-E-type L ⊗ M would
  still require realising it as a rank-r section with M_r ≅ L_r^#.
- What it does establish: **the smallest known candidate rank for a non-E-type
  witness drops from 84 to 64**, and that candidate has now passed every cheap
  test that kills the others.

## Next steps

1. The rank-64 5-modular genus: mass formula + Kneser/Schiemann neighbours.
   Note arXiv:2106.04828, *"A note on Quebbemann's extremal lattices of rank 64"* —
   produces many rank-64 extremal lattices from codes, with small automorphism
   groups, at exactly this rank.
2. Same positivity test for (ℓ=2, n=96), the other new candidate — needs the
   Γ₀(2)+ ring (weight-4 and weight-8 generators), not done here.
3. If a rank-64 extremal 5-modular lattice is found, check whether it embeds as
   a self-dual-up-to-similarity section of a tensor product.
