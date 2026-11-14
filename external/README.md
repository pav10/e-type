# External candidate: SRFP311T1 "entangled cyclotomic Lagrangian gluing", rank 96

Source: srfp311t1.com/extremal_lattice_96 (text supplied by the user; the site does not resolve
from this environment) and the two scripts `reconstruct.py` / `reconstruct_and_audit.py`
(supplied by the user, not stored here).

Files
- `srfp_G96.txt` — Gram matrix of the *seed-free* version, rebuilt by `scripts/srfp_rebuild.py`
  (their construction stages 1–4 copied verbatim; their PARI calls / apt-get installs omitted).
- `srfp_G96_seed142_uploaded.txt` — the seed-142 matrix `G96_reduced.txt` as supplied by the user
  (SHA-256 7a6a3b8f…; their page states 74611ca2…, possibly a formatting difference).
- `srfp_G96_norm8_vector.txt` — a vector of norm **8** in `srfp_G96.txt` (coordinates in that
  basis), found by BKZ-60 and rechecked exactly in two coordinate systems
  (glued basis: 8; raw basis: (7v)ᵀG_raw(7v)/49 = 392/49 = 8).

Findings
1. Both matrices: even, unimodular, rank 96, positive definite (checked).
2. Both contain vectors of norm 10 (BKZ-50), so min ≤ 10 (their own bound: ≤ 16).
3. Seed-free version: **minimum 8** (explicit vector) — not extremal.
4. Seed-142 version (the user's `G96_reduced.txt`): BKZ-70 finds a vector of norm **8** (norm computed
   exactly in that Gram matrix) — **also not extremal**. Vector being re-extracted to `srfp_G96_seed142_norm8_vector.txt`.
5. Their own audit prints "Extremal Min 10 Established: NO"; the page's theorem is conditional.
6. Their "min(Λ₈₈) ≥ 12" step is not a proof: the Schur-complement search keeps only ONE
   candidate tail vector (`qfminim(..., 1, 2)`) and lifts only that direction.
7. Structure: the glued lattice is stable under ω = (1+√−7)/2 (a rank-48 Hermitian lattice), but it
   is not a tensor product E ⊗ P′ (its intersection with each E-block is J itself, index 7²⁴),
   so the exact tensor method of this repository does not apply.
