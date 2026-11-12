# Positivity tests: all four candidates — and a correction

## Summary

| rank | level | extremal min | γ' vs √n | theta positivity | verdict |
|---|---|---|---|---|---|
| **64** | 5 | 18 | 8.0498 > 8.0000 | **PASSED** | live, smallest |
| 84 | 3 | 16 | 9.2376 > 9.1652 | **PASSED** | live |
| **96** | 2 | 14 | 9.8995 > 9.7980 | **PASSED** | live |
| 96 | 1 | 10 | 10.0000 > 9.7980 | (classical) | live |

All four survive. Existence remains open for every one of them.

Level death points (first negative extremal-form coefficient):
ℓ=23 → 6, ℓ=11 → 16, ℓ=7 → 36, **ℓ=5 → 88**, ℓ=3 → beyond 84, **ℓ=2 → none up to 160**.

---

## New this round: (ℓ=2, n=96) — PASSED

**Generators, both validated:**
- weight 2: θ(D₄) = 1 + 24∑σ₁^odd(n)qⁿ, confirmed against direct D₄ enumeration
  through q¹⁴: `[1,24,24,96,24,144,96,192,24,312,144,288,96,336,192]`. det(D₄)=4=2^{4/2} ✓
- weight 8: η(z)⁸η(2z)⁸, Fricke eigenvalue (−i)⁸ = **+1** ✓ (valid Γ₀(2)+ generator)

**Ring check:** reproduces min = 2+2⌊n/16⌋ at every rank, and two known kissing
numbers exactly — **D₄ → 24** and **Barnes-Wall Λ₁₆ → 4320**.

**Rank 96** (weight 48, dim 7, min 14):

        Θ = 1 + 256 628 736 000 q⁷ + 131 789 882 952 000 q⁸
              + 33 471 344 664 576 000 q⁹ + …

All coefficients non-negative out to q²⁵. Kissing number would be
**256 628 736 000** vectors of norm 14.

---

## CORRECTION to the rank-84 (ℓ=3) series I gave earlier

While cross-checking with a second generator pair, the two methods disagreed
from q¹⁰ onward. Tracing it: the earlier computation used θ(K12) values I had
**hardcoded from memory**, and one was wrong.

Direct Fincke–Pohst enumeration from the verified catalogue Gram (det 729) gives

        θ(K12) = [1, 0, 756, 4032, 20412, 60480, 139860, 326592, 652428, 1020096, **2000376**]

The last entry is **2 000 376**; I had used 1 693 512.

**Impact:**
- The linear solve for the extremal form uses only q⁰…q⁷ — all verified. The form
  itself was correctly determined.
- **q⁸ = 119 421 161 244 — CORRECT** (this is the kissing number, the key number).
- **q⁹ = 12 343 472 418 816 — CORRECT.**
- **q¹⁰ was WRONG.** Correct value **958 224 043 687 680** (I earlier said
  958 224 041 539 632). Everything beyond q¹⁰ in that earlier output is likewise
  unreliable.

Corrected series:

        Θ₈₄ = 1 + 119 421 161 244 q⁸ + 12 343 472 418 816 q⁹
                + 958 224 043 687 680 q¹⁰ + 47 464 428 506 148 864 q¹¹ + …

Two independent generator pairs now agree on these values.

**Lesson recorded:** the error came from hardcoded lattice data rather than from
the modular-forms machinery. Every generator should be enumerated or derived,
never recalled — the same failure mode as the first ℓ=3 attempt.

---

## Method note on generator validity

An eta quotient η(z)^k η(Nz)^k has Fricke eigenvalue (−i)^k, and only **+1**
eigenforms lie in M(Γ₀(N)+) where ℓ-modular theta series live:

| generator | k | eigenvalue | status |
|---|---|---|---|
| η(z)⁴η(5z)⁴ | 4 | +1 | valid (used for ℓ=5) |
| η(z)⁸η(2z)⁸ | 8 | +1 | valid (used for ℓ=2) |
| η(z)⁶η(3z)⁶ | 6 | −1 | **not** in M(Γ₀(3)+) |

So the ℓ=3 work must use genuine 3-modular theta series (θ(A₂)², θ(K12)) as
generators, which is what the corrected computation does. Curiously the invalid
η-basis still produced the right q¹⁰ value and is what exposed the K12 error —
but it should not be relied on.

---

## Where this leaves the programme

The smallest candidate rank for a non-E-type witness stands at **64** (level 5,
minimum 18, kissing number 21 281 280 000), 20 below the rank-84 target. It has
passed parity, LP exclusion, level persistence, and theta positivity.

Next: the rank-64 5-modular genus — mass formula plus Kneser/Schiemann
neighbours. arXiv:2106.04828 constructs many rank-64 extremal lattices from
codes at exactly this rank.
