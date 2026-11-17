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
   exactly in that Gram matrix) — **also not extremal**. Vector saved in `srfp_G96_seed142_norm8_vector.txt` (coordinates in the basis of the supplied file; norm rechecked with plain Python integers: 8).
5. Their own audit prints "Extremal Min 10 Established: NO"; the page's theorem is conditional.
6. Their "min(Λ₈₈) ≥ 12" step is not a proof: the Schur-complement search keeps only ONE
   candidate tail vector (`qfminim(..., 1, 2)`) and lifts only that direction.
7. Structure: the glued lattice is stable under ω = (1+√−7)/2 (a rank-48 Hermitian lattice), but it
   is not a tensor product E ⊗ P′ (its intersection with each E-block is J itself, index 7²⁴),
   so the exact tensor method of this repository does not apply.


# External candidate 2: SRFP311T1 "Four-Leech graph-code gluing in dimension 96" (Nov 2026)

Supplied by the user: manuscript `rank96_four_leech.pdf` and `verify_publication_core.py`.
The manuscript itself does not claim extremality (shells 4, 6, 8 stated as open).

- `four_leech_G96.txt` — their Gram matrix G, rebuilt by `scripts/four_leech_rebuild.py`
  (data and construction copied verbatim from their script; their later checks omitted).
  Checked: even, det 1, positive definite; their rootless certificate for P passes.
- `four_leech_short_vector.txt` — a vector of norm **8** in G (coordinates in their basis),
  found by BKZ-50 in 29 s; norm rechecked with plain Python integers.
- `four_leech_norm8_tuple.txt` — the same vector in their four-Leech coordinates:
  component P-norms (|a1|², |a2|², |b1|², |b2|²) = (12, 22, 10, 12), total 56,
  b ≡ F a (mod 7), and (56)/7 = 8 by their formula (1.2).  This violates their
  criterion (8.1) (total ≥ 70), so **min(L_F) = 8: the candidate is not extremal.**
