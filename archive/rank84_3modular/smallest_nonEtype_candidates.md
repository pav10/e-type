# Smallest-rank candidates for a non-E-type witness

## The logic (why extremal modular lattices are the right search space)

Kitaoka: if L ⊗ M has a **non-split** minimal vector of rank r, then
r ≤ γ_r²/ρ with ρ ≥ 1, so a **necessary condition** is

        γ_r ≥ √r.

Coulangeon (Acta Arith. 92, Prop 3.3) equality case: such a witness needs
rank-r sections with M_r ≅ L_r^#, i.e. the section is **similar to its dual** —
modular or isodual. For an ℓ-modular lattice of rank n, det = ℓ^{n/2}, hence

        γ = γ' = min / √ℓ,

so the witness condition min/√ℓ ≥ √n becomes

        **min² ≥ ℓ · n.**

This is exactly why the rank-84 3-modular (min 16) and rank-96 unimodular
(min 10) targets are the objects of interest:
  ℓ=3, n=84:  16/√3 = 9.2376 > √84 = 9.1652
  ℓ=1, n=96:  10/1    = 10.0000 > √96 = 9.7980

## Method

Quebbemann extremal bound: min = 2 + 2⌊n/m_ℓ⌋ with jump m_ℓ = 48/(ℓ+1).

**Formula validated 5/5 against Nebe's ExtLat table:**

| ℓ | n | predicted min | known lattice |
|---|---|---|---|
| 2 | 16 | 4 | Barnes-Wall Λ16 |
| 3 | 12 | 4 | Coxeter-Todd K12 |
| 5 | 16 | 6 | N16 (unique extremal 5-modular) |
| 11 | 4 | 4 | det 11², min 4 |
| 23 | 2 | 4 | det 23, min 4 |

Filters applied:
- **(a) parity** — required divisor of n per level, read off the table's `-` pattern
  (ℓ=1: 8 | n; ℓ=2,5,11: 4 | n; ℓ=3,7,23: 2 | n).
- **(b) level persistence** — the table's `-` marks where the extremal modular form
  acquires a **negative coefficient**, so no extremal lattice exists. Levels
  **7 (from dim 36), 11 (from dim 16), 23 (from dim 6)** die well below rank 52.
  Only **ℓ = 1, 2, 3, 5** are still `?` (no negativity) at the top of the table.
- **(c) LP exclusion** — γ_r < √r is certified for r ≤ 48, borderline 49–51,
  and fails from ~56. So no candidate below ~52 can be real.
- **(d) tightness** — a *large* margin over √n means the Quebbemann bound has
  overshot the geometric (Cohn–Elkies) ceiling, so the lattice cannot exist.
  **Small margin = plausible.** This is the key discriminator.

## Result: the four tightest candidates

| rank n | level ℓ | extremal min | γ' | √n | margin | status |
|---|---|---|---|---|---|---|
| **64** | **5** | **18** | 8.0498 | 8.0000 | **+0.0498** | **smallest & tightest — NEW** |
| 84 | 3 | 16 | 9.2376 | 9.1652 | +0.0725 | your target |
| 96 | 2 | 14 | 9.8995 | 9.7980 | +0.1015 | alternative at rank 96 — NEW |
| 96 | 1 | 10 | 10.0000 | 9.7980 | +0.2020 | your target |

Everything else in ranks 52–100 has margin ≥ 0.22 and is progressively more
likely vacuous (Quebbemann bound above the geometric ceiling).

## The headline

**A rank-64 extremal 5-modular lattice of minimum 18 would be a witness
20 dimensions below the rank-84 target** — and it is simultaneously the
*tightest* candidate in the whole range, which is the signature of a bound
that has not yet overshot. Both facts point the same way.

Secondary finding: at rank 96 there is a **level-2** route (min 14) that is
tighter than the unimodular one (min 10), i.e. a-priori more plausible.

## What is NOT established (the immediate next checks)

1. **Theta positivity at (ℓ=5, n=64) has NOT been verified here.** Nebe's table
   stops at dim 48 and shows `?` (not `-`) for ℓ=5 at 32–48, so no negativity
   has set in by 48 — encouraging, not proof. Compute the extremal weight-32
   form for Γ₀(5)+ and check every coefficient ≥ 0. **If any is negative the
   candidate dies immediately**, exactly as levels 7/11/23 do.
   (Scharlau–Schulze-Pillot proved each level's extremal lattices are finite in
   number, so ℓ=5 *does* die somewhere; the only question is whether before 64.)
2. Existence, if positivity holds, remains open — as for 84 and 96.
3. **γ' ≥ √n is necessary, not sufficient.** Such a lattice removes the
   obstruction but does not by itself produce a non-E-type example: one would
   still have to exhibit L ⊗ M with an actual non-split minimal vector, using
   the lattice as the rank-r section.

## Suggested order of attack

1. Extremal theta for (5, 64) — cheap, decisive, kills or keeps the target.
2. If positive: the ℓ=5 genus at rank 64 is the search space. Note the lead
   *"A note on Quebbemann's extremal lattices of rank 64"* (arXiv:2106.04828),
   which produces many rank-64 extremal lattices from codes with small
   automorphism groups — directly relevant machinery at exactly this rank.
3. Same positivity check for (2, 96), the other new candidate.
