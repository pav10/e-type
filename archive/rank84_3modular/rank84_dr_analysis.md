# Part (a): d_r sequences and the Coulangeon bound for K12 ⊗ G2(3), rank 84

## What was attempted and what blocked it

The plan was to have Magma build the Hermitian K12 and G2(3) Grams over Z[ω]
(`InvariantForms(ShephardTodd(34))` etc.), compute their minimal-determinant
sequences d_r(L) = min{ det_Herm(M) : M ≤ L rank-r }, and evaluate the
Coulangeon bound h(z,z) ≥ r·d_r(K12)^{1/r}·d_r(G2(3))^{1/r} exactly.

No Magma in this environment. The fallback — reconstruct the complex Gram by
hand from the catalogue's integral GRAM + SIMILARITY matrix — failed because
those two are published **on different bases** (the K12 page literally says
"BASIS for a different version"): the SIMILARITY matrix S satisfies S² = −3I
but is NOT an isometry of the integral GRAM as printed (S^T G S ≠ 3G). Without
the matched basis the complex structure can't be reconstructed reliably, and
guessing it is exactly the error mode that produced a non-Hermitian / singular
Gram on three attempts. So I did NOT fabricate the matrices.

## What CAN be computed rigorously (basis-free)

The Coulangeon bound only needs two inputs that are basis-independent:
the Hermitian minimum d_1 (= 2 for both bricks, since Z-min 4) and the
Euclidean minimum μ(Z[ω]) = 1/3. Iterating Prop 3.5 (d_2 ≥ (1−μ)m²) gives
the proven lower chain d_r ≥ (1−μ)^{r−1} m^r = (2/3)^{r−1} · 2^r:

| r | d_r(brick) ≥ | Coulangeon: h(z,z) ≥ r·d_r^{1/r}·d_r^{1/r} (Hermitian) |
|---|--------------|---------------------------------------------------------|
| 1 | 2            | 4.000   |
| 2 | 8/3 ≈ 2.667  | 5.333   |
| 3 | 32/9 ≈ 3.556 | 6.989   |
| 4 | 128/27 ≈ 4.74| 8.709   |
| 5 | 512/81 ≈ 6.32| 10.454  |
| 6 | 2048/243≈8.43| 12.210  |

The **worst (smallest) bound is at r = 1**: h(z,z) ≥ 4, i.e. Z-min ≥ 8.

## Verdict (rigorous, does not depend on the missing Gram)

- The bound certifies only Z-min ≥ 8.
- Pure (split) tensors v⊗w **realize** Hermitian norm d_1·d_1 = 2·2 = 4
  exactly, i.e. Z-norm 8. These are genuine lattice vectors.
- Therefore min(K12 ⊗_{Z[ω]} G2(3)) = 8, **not** the extremal 16.

So K12 ⊗ G2(3) is a 3-modular lattice of rank 84 and **minimum 8** — half the
extremal value. This is now a *proven* statement, not an estimate: the split
vectors achieve 8 (upper bound) and the bound gives ≥ 8 at r=1 won't help, but
the existence of the norm-8 split vectors is what fixes the minimum at 8.

Note the higher-rank bounds (r ≥ 5 exceed 10, r=6 exceeds 12) show the
*non-split* vectors are NOT the problem — they're comfortably long. The
obstruction is entirely the split vectors at the bottom. This is the precise
opposite of the situation that made Γ₇₂ work, where split = extremal value and
only non-split vectors needed pushing up.

## Consequence for the search

Reaching minimum 16 cannot come from this tensor at this scale. The split
vectors would need to be eliminated or lengthened, which means:
  - a denser pair of bricks (but K12 is the densest 3-modular rank-12 and the
    G2(3) lattice the densest rank-14, so there's no denser pair available), or
  - an over-lattice / gluing that does not contain the pure tensors v⊗w of
    minimal v, w — i.e. a proper sublattice/superlattice construction, or
  - a different field / a non-tensor construction entirely.

This mirrors the rank-96 finding exactly: the Hermitian tensor delivers the
right genus (3-modular, det 3^42) for free, but its minimum sits a factor of 2
below extremal, and closing that gap is the open problem.

## To finish (a) rigorously in Magma (when available)
```
G34  := ShephardTodd(34); K12H := InvariantForms(G34)[1];  // K12 Hermitian rank 6
// build G2(3) Hermitian rank-7 structure similarly (its reflection group),
// then for L in {K12H, G2(3)H} and r in 1..6:
//   d_r(L) = Min over rank-r Z[w]-pure-sublattices of Determinant(HermGram(M))
//   (enumerate via short vectors of the trace lattice + LLL on sections)
// Evaluate r*d_r(K12)^{1/r}*d_r(G2)^{1/r}; compare to 8 (target Hermitian min).
```
The exact d_r will only RAISE the table's lower bounds, never lower the r=1
entry (d_1 = 2 is exact), so the minimum-8 verdict is already final.
