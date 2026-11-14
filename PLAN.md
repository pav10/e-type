# Non-E-type lattices at small rank — plan, code, first results

Targets: an **extremal even unimodular lattice of rank 96** (min 10) and/or an
**extremal 3-modular lattice of rank 84** (min 16). Either one is an explicit
non-E-type lattice, because `min(L)·min(L#) > rank` (100 > 96, 256/3 > 84) and the
identity tensor `Σ e_i ⊗ e_i#` has norm = rank.

Everything here is Python (cypari2 + fpylll) or PARI/GP; nothing needs Magma.
Run `./setup.sh`, then `python3 scripts/validate.py` (≈ 1 min) before anything else.

**Status in one paragraph.** The uploaded material was read critically (§1). A
tested pipeline now builds Hermitian tensor products and ideal lattices and
tries to kill them (§2, §3). It reproduces Coulangeon–Nebe's `P48n = P ⊗ T`
over ℤ[(1+√−11)/2] as a positive control. A first systematic run of the rank-96
family `A ⊗ Leech` over ℚ(√−d), d = 19, 20, 23, 24 killed all but a handful of
candidates; ten left over at d = 19 are **undecided** (§4): they contain vectors
of norm 10 and no vector of norm ≤ 8 was found, but the tests run so far do not
cover tensors of rank 3 and 4, so they cannot yet separate minimum 8 from 10. Deciding them is step 1 of the plan (§5).

---

## 1. Critical review of the brief and the uploaded files

| # | Claim in the brief / skill | Verdict |
|---|---|---|
| 1 | Extremal theta data: rank 84 ℓ=3 kissing 119 421 161 244; rank 64 ℓ=5 kissing 21 281 280 000; level deaths ℓ=5 → 88, ℓ=7 → 30 | **Confirmed** independently (`etype/theta.py`, exact integers). Rank 96: kissing **565 866 362 880**. |
| 2 | Error log #4: "η(z)⁶η(3z)⁶ has Fricke eigenvalue −1, so it is not admissible; only +1 forms qualify" | **Misstated.** θ of an ℓ-modular lattice of rank 2k has eigenvalue (−i)^k, and a weight-k generator must have eigenvalue (−i)^k. Δ₃ = η(z)⁶η(3z)⁶ has −1 = (−i)⁶ and *is* a generator (θ_K12 = θ_A2⁶ − 36·… uses it). My module uses Δ₃ and reproduces the brief's numbers, so the numbers are right; the rule as written is not. |
| 3 | Dead end #1: "Hermitian tensor products give half the extremal minimum" | **Over-claimed.** Shown only over ℤ[ω]. In general split vectors of `A ⊗_{O_K} B` (A Hermitian unimodular, B in the target genus) have norm `hmin(A)·min(B)`, and `hmin(A)` grows like √d. Over ℚ(√−19) there is a rank-4 A with hmin 3, so split vectors of `A ⊗ Leech` have norm 12 > 10. The obstruction is the *non-split* vectors, exactly as in Γ72. |
| 4 | Dead end #2: "cyclotomic ideal lattices are capped by root discriminant, γ ≥ n/rd" | γ ≥ n/rd is a **lower** bound; "capped" is an empirical observation on 17 cyclotomic fields. Non-cyclotomic abelian CM fields (composita like ℚ(ζ28, √−19)) were not covered, and they are where today's candidates live. |
| 5 | First action: `Mass(L)` for rank 4 over ℚ(ζ17) decides the p = 17 route | **Decided: dead for enumeration.** Standard part 4·10¹⁷ (exact L-values; check: ∏L(0,χ_odd) = 256/34), local factor at 5 ≈ 5³², mass ≈ 10⁴⁰ (`scripts/mass_p17.py`; order of magnitude, factor at 17 ignored). |
| 6 | Skill P2′: Kneser-neighbour random walk in the rank-84 genus | **Hopeless, not merely large.** A random lattice of the genus has ≈ 10¹⁰ vectors of norm < 16 (Siegel mean); an unstructured walk never sees an extremal lattice. Only structured families make sense. Neighbour steps are still useful *inside* small Hermitian genera (bricks). |
| 7 | "Use BKZ / `ShortVectors` as the gate" | BKZ kills only when short vectors are plentiful. With the flat profile of these lattices, a norm-8 vector in a rank-96 unimodular lattice needs block size ≈ 90 to be found. **A BKZ "survivor" is not evidence of extremality.** |
| 8 | Uploaded `.magma` / `.m` files | None was ever run; all bricks are placeholders. Superseded by this repository. |
| 9 | "Rank 96 is a pure packing problem" (skill P4) | The tensor route applies there too; it is where the first undecided candidates are. |
| 10 | LP table (γ₄₄…γ₇₂) and "43 → 48" | Not re-verified: cohn.mit.edu is not reachable from this environment. See §6. |

Prior art read in full: Coulangeon–Nebe, arXiv:1201.1832. Relevant facts: Prop. 3.3
(`h(z,z) ≥ r·(d_r(L)d_r(M))^{1/r}` for a rank-r tensor, equality iff dual minimal
sections), §4 (`P48n = P₃ ⊗ T` over ℤ[(1+√−11)/2]; the Hermitian structures of
the Leech lattice over that ring are **not classified**), Thm 5.3 (Γ72: minimum 8
iff the Leech structure does not contain the Barnes lattice — 1 of 9 structures).
The paper does not treat rank 96 or 84.

---

## 2. The constructions (what the code can build)

**Convention.** K = ℚ(√−d), O_K = ℤ[η], η² − tη + N = 0. A lattice with O_K-structure
is a pair `(G, J)`: Gram matrix and the integer matrix of η. `G = Tr∘h`.

**(T) Hermitian tensor products** `A ⊗_{O_K} B` (`tensor_OK`). If A is Hermitian
unimodular (trace lattice of det d^m) and B has an O_K-structure, then `A ⊗ B`
lies in the genus of `B^m`. Necessary condition from split vectors:
`hmin(A)·min(B) ≥ target minimum`.

Rank 96 (B even unimodular with O_K-structure):

| rank_{O_K} A | B | need hmin(A) ≥ | first d where A exists (walks, `scripts/brick_survey.py`) |
|---|---|---|---|
| 2 | extremal rank 48 (min 6) | 2 | d = 11, 15, 19, 20, … |
| 3 | rank 32, min 4 | 3 | d = 20, 23, 31 |
| 4 | Leech | 3 | d = 19, 20, 23, 24, 31 (not 3, 7, 8, 11, 15) |
| 12 | E8 | 5 | open (ℤ-rank 24, min 10, det d¹²) |

Rank 84 (B 3-modular via an element c ∈ O_K with c·c̄ = 3, so K ∈ {ℚ(√−2), ℚ(√−3), ℚ(√−11)}):
the only arithmetically possible brick pair is rank-2 A with hmin 2 (exists for
d = 8, 11) and **B = an extremal 3-modular lattice of rank 42 (min 8) with
O_K-structure**, itself not known to exist. So for rank 84 the tensor route
reduces to a rank-42 search; ideal lattices over degree-84 fields are the
direct alternative.

**(I) Ideal lattices** over a CM field L (`gp/ideal_lattices.gp`): pairs (I, α) with
`α·I·Ī·Diff_L = (c)`; c = 1 unimodular, c = √−d Hermitian unimodular over O_K,
c with c·c̄ = ℓ gives ℓ-modular. Full enumeration over the class group works up
to degree ≈ 32. For degree 84/96, `twist(nf, s, I, Q)` replaces I by `I·Q/Q̄`,
which preserves the equation for *every* ideal Q: the family can be sampled
with `nfinit` only, no class group.

**(N) Kneser neighbours** of Hermitian lattices at split primes (`neighbour`),
used to enumerate brick genera.

**Kill tests.** `bkz_shortest` (exact norm of the vector found; upper bound on
the minimum) and `rank2_kill` (all tensors `a₁⊗p₁ + a₂⊗p₂` on sections spanned
by short vectors; exact minima of the 8-dimensional sublattices).

---

## 3. Validation (all passing)

- Theta ring: E8 240, Leech 196 560, P48 52 416 000, Γ72 6 218 175 600, D4 24,
  BW16 4320, K12 756 and its full θ to q¹⁰, Barnes 42, ℓ=5 rank 8: 120.
- Tensor: `O⊗O`, `T⊗T` (min 4 = the identity tensor, as predicted).
- Ideal lattices over ℚ(ζ13, √−11): class group C35, even unimodular rank 24;
  root-free ones are Leech with a ℤ[η]-structure.
- **Control (C–N Cor. 4.2):** `P ⊗ T`, rank 48: for some structures P the kill
  test finds norm 4 at once, for the others nothing below 6 (the extremal P48n).

---

## 4. First results

### 4.1 Rank 96: `A ⊗_{O_K} Leech`, A of rank 4, hmin 3

Leech structures: ideal lattices over ℚ(ζ_m, √−d), φ(m) = 12 (m = 13, 21, 28, 36).
Counts are ideal-lattice classes, not isometry classes (many coincide).

| d | rank-4 bricks with hmin ≥ 3 (ℤ-min, kissing) | Leech structures used | candidates | killed by LLL/BKZ-50 (norm found) | not killed |
|---|---|---|---|---|---|
| 19 | (6,20), (6,24) | 420 | 840 | 827 (norm 4, 6 or 8) | **13**, all with brick (6,20) |
| 20 | (6,16) | 359 | 359 | 359 | 0 |
| 24 | (6,16) | 348 | 348 | 348 | 0 |
| 23 | (6,8), (6,12), (6,14), (6,16), (6,24), (8,66) | 926 | 5556 | 5556 (norm 4, 6 or 8) | 0 |

Follow-up on the 13 (d = 19):
- `scripts/rank2_test.py` (tensors `a₁⊗p₁ + a₂⊗p₂` on sections spanned by minimal
  vectors, exact 8-dimensional minima) **kills 3 more** with explicit norm-8
  vectors (Leech #93, #136, #318) — vectors BKZ-50 had missed, which is the
  practical proof of row 7 of §1.
- The remaining **10** (Leech #30, 81, 82, 109, 149, 160, 259, 260, 366 from
  m = 28 and #364 from m = 13) show the same 66 section types and are probably
  one or two isometry classes. For #81, `scripts/section_test.py` searched all
  50 sublattices `A₂ ⊗ P` (dimension 48) for sections A₂ spanned by vectors of
  norm ≤ 8: smallest norm found is **10**, attained on the four sections of
  determinant 4. BKZ up to block 70 on the full lattice also returns 10.
- The ten have the same Hermitian automorphism group order, **|Aut_{O_K}(P)| = 4368 =
  2·|SL₂(13)|** (PARI `qfauto` with the forms G and G(2J−1)), the order of
  Coulangeon–Nebe's P₃ (Aut SL₂(13).2), the structure that gives P48n over ℚ(√−11).
  They are very probably one lattice. #0 (killed) has |Aut| = 52.

#### 4.1.1 Exact decision procedure (running for #81)

For z ∈ A ⊗ P of rank r there is a primitive rank-r section A_r of A with
z = Σ bᵢ ⊗ xᵢ, (bᵢ) a basis of A_r, and h(z,z) = tr(W_r · X̄), where W_r = Gram(A_r)
and X = Gram(x₁,…,x_r) (positive definite, entries in 𝔇⁻¹, integer diagonal ≥ 2,
no vector of norm 1 in (O_K^r, X)).
- r = 1: h ≥ 3·2 = 6 > 4. Excluded.
- r ≥ 2: Coulangeon's Prop. 3.2 plus Hermite bounds for d_r(P) (Z-sections of
  Leech have minimum 4; γ₄ = √2, γ₆ = (64/3)^{1/6}) give det A₂ ≤ 9.50 and
  det A₃ ≤ 14.18. Enumerated completely: **19** classes of rank-2 sections
  (via a minimal vector w₁ and a second vector bounded with the covering radius
  1.316 of O_K) and **233** classes of rank-3 sections (w^⊥ for h(w,w) ≤ 14).
- For W = W_r, `etype/certify.py: enumerate_X` lists every X with tr(W X̄) ≤ 4
  (Fincke–Pohst in the Cholesky coordinates of W⁻¹; finite and exact).
  Result: 23 (r = 2), 88 (r = 3), 239 337 (r = 4) shapes.
- Necessary filter: the Hermitian products between minimal vectors of
  (O_K^r, X) must occur between minimal vectors of P. The complete table
  for P (all 196 560² pairs) has 67 values; the filter removes 94 182 shapes.
- O_K-isometry classes (`qfisom` on the pair (G, G(2J−1))): **7926** (2 / 15 / 7909).
- For each class: exhaustive search for an O_K-linear isometric embedding into P.
  The frame is chosen from minimal vectors so that the prescribed products are
  rare in P; the first frame vector runs over the **45 orbits of Aut(P)** on the
  minimal vectors; the final map is checked to be integral.

**z of norm ≤ 8 exists ⟺ some class embeds.** All classes embed nowhere ⟹ min(A ⊗ P) = 10,
i.e. an extremal even unimodular lattice of rank 96. Classes whose minimal
vectors do not span are reported separately and need norm-6 frame vectors.

Control (`scripts/control_rank48.py`, Coulangeon–Nebe Prop. 4.1): for P ⊗ T over
ℚ(√−11) the same pipeline enumerates exactly the three shapes T, L_K, L̄_K, finds an
embedding for structures with min 4, and finds none for the SL₂(13).2 structure
(P48n). **Passed** for the negative case (P48n: no shape embeds). Positive case on d = 19
itself (`scripts/control_d19.py`, production brick and shapes): for killed candidates
whose BKZ norm-8 vector has tensor rank 4, the vector's X occurs *exactly* in the
production list and the search finds its embedding (**passed**, Leech #56, #71).
Also checked numerically: h(z,z) = tr(W̄·X) with the conjugation used in the code
(an earlier version of the d = 19 control failed only because it used a rank-2 vector).

**Result: DECIDED, minimum 8.** The embedding search found 14 shape classes (all rank 4)
that embed into P. Each gives an explicit z ∈ A ⊗ P (`scripts/verify_kill.py`):
integral coordinates in the lattice basis, norm **exactly 8** in the 96-dimensional
Gram matrix (three verified and stored in `results/decide_d19_81/kill_vectors.pkl`).
The other nine undecided structures (#30, 82, 109, 149, 160, 259, 260, 364, 366) are
O_K-isometric to #81 (`qfisom`), so they give the same lattice.

**Conclusion: every A ⊗ Leech candidate over ℚ(√−19), ℚ(√−20), ℚ(√−23), ℚ(√−24)
has minimum ≤ 8.** None is extremal. Lessons:
- BKZ up to block 70 did not find these norm-8 vectors; the exact shape+embedding
  method did. For every future candidate, the exact method is the gate.
- The method is general: it decides min(A ⊗_{O_K} P) for any brick A of small rank
  and any P with a large automorphism group, at a cost of hours.

### 4.2 Bricks (rank 4, Hermitian unimodular) — complete genera, checked against Schiemann

`scripts/bricks.py` walks from O_K⁴ *and* from rank-4 ideal lattices over ℚ(ζ_k)·K
(these reach the even genera), with classes up to O_K-isometry (`qfisom`). Class numbers
agree exactly with Schiemann, J. Symb. Comp. 26 (1998), Tables 1–2:

| K | genus | classes found | Schiemann | best Hermitian min (classes) |
|---|---|---|---|---|
| ℚ(√−5), d = 20 | odd / even | 18 / 13 | 18 / 13 | 3 / **4** (3 classes, kissing 120) |
| ℚ(√−6), d = 24 | odd / even | 28 / 16 | 28 / 9 + 7 (two spaces) | 3 / **4** (4 classes) |
| ℚ(√−31) | I | 216 | 216 | 4 (3 classes) |
| ℚ(√−43) | I | 69 | 69 | **5** (1 class) |
| ℚ(√−67) | I | 301 | 301 | **6** (2 classes); 5 (30 classes) |
| ℚ(√−11) | I | — | 6 | 2 → `A ⊗ Leech` over ℚ(√−11) impossible |

### 4.3 Rank 96 with the new bricks (BKZ-50 pre-filter, then exact test)

BKZ-50 kills (rigorous: explicit vectors) and survivors, `A ⊗ Leech`:

| K | bricks used (hmin) | Leech structures | candidates | killed by BKZ | BKZ survivors |
|---|---|---|---|---|---|
| ℚ(√−5) | 3 even bricks (4) | 359 | 1077 | 1035 | 42 (brick #3: 2, #4: 24, #5: 16) |
| ℚ(√−6) | 4 even bricks (4) | 348 | 1392 | 1378 | 14 (all brick #6, not free) |
| ℚ(√−43) | 1 brick (5) | 357 (m = 28, 21) | 357 | 357 (mostly norm 6) | 0 |
| ℚ(√−43) | 12 bricks (4) | 357 | 4284 | running | 424 of the first 3132 |
| ℚ(√−67) | 2 bricks (6) | — | — | Leech structures not yet computed (PARI stack) | — |

A larger brick minimum does **not** help: the hmin-5 brick over ℚ(√−43) gives
vectors of norm 6, shorter than the hmin-3/4 bricks.

Exact follow-up (BKZ survivors are not evidence; §4.1.1):
- ℚ(√−5), brick #4 ⊗ Leech #111 (|Aut_h(P)| = 84): **killed**, explicit vector of
  norm 8 (`scripts/kill_stream.py`, 3397 shape classes searched).
- The killing shape transferred: Leech #145 **killed** (explicit norm-8 vector);
  for the other brick-4 survivors tested so far the shape does not occur — those
  stay undecided until their own kill test or full decision.
- Brick #3 survivors (|Aut_h(P)| = 2184 = |SL₂(13)|, the same symmetry as the dead
  d = 19 lattice): kill test running.
- BKZ with rerandomised bases is *weaker* than plain BKZ here (d = 19 control:
  best 12 over 12 tries vs. 10), so no cheap heuristic screen is available.

### 4.4 p = 17 route for (ℓ = 5, rank 64)

Mass ≈ 10⁴⁰: genus enumeration is out. A prescribed *group* or a tensor/ideal
construction is the only way in; the pipeline applies unchanged (5 = c·c̄ in
ℚ(√−1), ℚ(√−11), ℚ(√−19), ℚ(√−5)).

---

## 5. Plan

**Step 1 — decide the undecided rank-96 candidates.** **Done** (§4.1.1): minimum 8.
Stored in `results/SURVIVOR_rank96_d19_*.json` (regenerate with
`scripts/rank96_tensor.py 5 1`). In order of cost:
1. Identify the isometry classes among the 10, then extend `rank2_kill` /
   `section_test.py` to rank 3 (72-dimensional sublattices `A₃ ⊗ P`) and rank 4,
   and to sections spanned by longer vectors. By C–N Prop. 3.3
   a vector of norm 8 of rank r needs sections with `d_r(A)·d_r(P) ≤ (4/r)^r`;
   r = 1 is excluded, r = 2, 3, 4 are finite enumerations.
2. A full sieve (G6K) in dimension 96 on a larger machine. Heuristic, but it
   returns essentially all vectors up to √(4/3)·λ₁, so a norm-8 vector would show.
3. If still alive: proof of minimum 10 along C–N §5 (bounds on `d_r` of the
   Leech structure + representation numbers).
Expectation, stated plainly: an extremal rank-96 lattice would be a major
result, and the prior is strongly that these candidates have minimum 8.

**Step 2 — widen the rank-96 tensor family.** (Now the main line; hours each.)
- Bricks of Hermitian minimum 4 from Schiemann's *even* genera (ℚ(√−5), ℚ(√−6), rank 4)
  and from ℚ(√−23) / ℚ(√−31) (rank 4, min 4): split vectors of norm 16, so only the
  non-split ones matter. Generate them by neighbour walks from an even start lattice.
- Leech structures with a large Hermitian automorphism group: by Nebe (JNT 161, 2016,
  Thm 3.8), an extremal 96-dimensional lattice with an automorphism of prime order
  p > 7 has p ∈ {13, 17} acting fixed-point-freely, or types 11−(8,16)−8 / 13−(7,12)−7.
  Structures with SL₂(13)-type symmetry (as here) or order-17 symmetry are the natural
  ones; the dead candidate had exactly the order-13 symmetry.
- More fields: `scripts/rank96_tensor.py N t` for d = 31, 35, 39, 40, 43, 47.
- More Leech structures: m with φ(m) = 12 covers only structures with a cyclic
  symmetry; add neighbour walks in the rank-12 Hermitian genus.
- Rank-2 bricks × extremal rank-48 lattices with O_K-structure (ideal lattices
  over ℚ(ζ_m, √−d), φ(m) = 24: degree 48, `bnfinit` is slow — use `twist`).
- Rank-3 bricks × rank-32 lattices of minimum 4 (d = 20, 23, 31).

**Step 3 — rank 84.**
- Search extremal 3-modular lattices of rank 42 with ℤ[(1+√−11)/2]- or
  ℤ[√−2]-structure (ideal lattices over ℚ(√−11)·ℚ(ζ43)⁺ and ·ℚ(ζ49)⁺, degree 42),
  then tensor with the rank-2 brick.
- Ideal lattices over degree-84 CM fields containing c with c·c̄ = 3
  (ℚ(ζ129), ℚ(ζ147), ℚ(ζ43, √−11), ℚ(ζ49, √−11), ℚ(ζ43, √−2), …) via `twist`.

**Step 4 — rank 64, ℓ = 5** with the same code (bricks × 5-modular lattices of
rank 16/32 with O_K-structure).

**Step 5 — E-type bound.** See §6.

---

## 6. Extending Kitaoka's bound (43 → ?)

Kitaoka's argument gives E-type for rank ≤ R as soon as γ_r < √r for 2 ≤ r ≤ R.
With the linear-programming bounds this is claimed for R = 48, and the trend in
the brief's table crosses zero near r ≈ 51.4. To turn this into a statement:
1. obtain the certified LP bounds for r = 44…51 (Cohn's table; not reachable
   from here) or recompute them with rigorous arithmetic;
2. r = 49, 50, 51 have margins of roughly 1.4 %, 0.8 %, 0.2 % in γ — the last
   needs a near-optimal, rigorously verified LP certificate;
3. beyond 51 the LP bound itself exceeds √r, so new input is needed (the
   "defect estimate" of the skill notes, or any bound on γ′ₙ below γₙ).
This track is independent of the constructions and not started in code.

---

## 7. Literature that would help (please send if available)

- Cohn's sphere-packing bound table (cohn.mit.edu/sphere-packing) — for §6.
- Hentschel–Krieg–Nebe, *On the classification of even unimodular lattices with
  a complex structure* (Int. J. Number Theory 2012) — Leech structures, masses.
- Nebe, *On automorphisms of extremal even unimodular lattices* (arXiv:1008.2862
  and 1212.0865) and anything on dimension 96 — to prune automorphism types.
- Coulangeon, *Tensor products of Hermitian lattices*, Acta Arith. 92 (2000).
- Schiemann, *Classification of Hermitian forms with the neighbour method* —
  to certify completeness of the brick lists.

---

## 8. Files

| path | purpose |
|---|---|
| `setup.sh` | dependencies + fplll BKZ strategies |
| `gp/ideal_lattices.gp` | ideal lattices over CM fields, twisting, Kronecker product (pure GP) |
| `etype/core.py` | Gram tools, Hermitian ↔ ℤ, `tensor_OK`, `neighbour`, `bkz_shortest` |
| `etype/search.py` | `ideal_family`, `leech_structures`, `brick_walk`, `rank2_kill` |
| `etype/theta.py` | exact extremal theta series, ℓ ∈ {1, 2, 3, 5, 7} |
| `scripts/validate.py` | validation harness (run first) |
| `scripts/brick_survey.py` | Hermitian minima of small unimodular genera |
| `scripts/rank96_tensor.py` | the rank-96 `A ⊗ Leech` search for one field |
| `scripts/attack.py`, `scripts/rank2_test.py`, `scripts/section_test.py` | heuristic kill tests (superseded by the exact method) |
| `etype/certify.py` | exact enumeration of Gram shapes X with tr(W̄X) ≤ T |
| `etype/represent.py` | Hermitian automorphisms/isometries (qfauto/qfisom), exhaustive embedding search |
| `scripts/shapes.py`, `scripts/decide.py` | exact decision of min(A ⊗ P) ≤ 8: sections → shapes → classes → embeddings |
| `scripts/verify_kill.py` | turns an embedding into an explicit short vector and checks its norm exactly |
| `scripts/control_rank48.py`, `scripts/control_d19.py` | negative / positive controls of the exact method |
| `scripts/leech_types.py`, `scripts/leech_with_zeta.py` | inner-product table of minimal vectors; ζ_m-action of ideal-lattice structures |
| `scripts/mass_p17.py` | mass estimate for the p = 17 genus |
| `results/`, `logs/` | outputs of the runs reported in §4 |
| `archive/` | the uploaded brief and files, unchanged |
