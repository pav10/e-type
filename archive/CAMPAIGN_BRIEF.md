# Non-E-type lattice search — campaign brief

Handoff for a Claude Code session. Everything below was established or corrected
in a long chat session; the error log at the end is as important as the results.

---

## 1. The goal, and the criterion

Find an **explicit** lattice that is **not of E-type** (Kitaoka), at the smallest
possible rank. Kitaoka: rank ≤ 43 ⟹ E-type; Pavlo's LP work extends this to 48.
Steinberg: non-constructive existence from rank 292. **No explicit example is
known at any rank** (Coulangeon–Nebe 2012, Remark 2.4(1)).

**Sufficient criterion (verified).** For any lattice $L$ of rank $n$ with dual
bases $\{e_i\},\{e_i^\#\}$, the vector
$$z=\sum_{i=1}^n e_i\otimes e_i^\#\in L\otimes L^\#$$
has rank $n$ and $(z,z)=\operatorname{tr}(GG^{-1})=n$. So

$$\boxed{\ \min(L)\cdot\min(L^\#)>n\ \Longrightarrow\ L \text{ is not of E-type}\ }$$

witnessed explicitly by $z$. For an $\ell$-modular lattice this is $\min^2>\ell n$,
i.e. $\gamma=\gamma'=\min/\sqrt\ell>\sqrt n$. Mechanism = Coulangeon, *Tensor
products of Hermitian lattices*, Acta Arith. 92 (2000), Prop. 3.2 equality case.

**Structural fact.** The required minimum $R$ = least even integer with
$R^2>\ell n$ satisfies $R\in(\sqrt{\ell n},\sqrt{\ell n}+2]$, so
$\gamma_{\text{req}}\in(\sqrt n,\sqrt n+2/\sqrt\ell]$. Every witness, at every
rank and level, needs a lattice **essentially at the LP bound**. Higher $\ell$ =
tighter requirement. This is why the LP certification $\gamma_r<\sqrt r$ and the
witness hunt are the *same computation*.

---

## 2. PRIMARY TARGET

> ### Extremal 5-modular lattice, dimension 64
> - $\det = 5^{32}$, minimum **18**, kissing number **21 281 280 000**
> - $\gamma = 18/\sqrt5 = 8.0498 > \sqrt{64}=8$
> - LP density bound at $n=64$: $\gamma_{64}\le 8.5948$ — **not excluded**
> - best known $\gamma_{64} = 6.928$ (Nebe, cyclo-quaternionic) ⟹ **a ~16% sphere-packing density record**

Why this one and not the smallest surviving candidate: its theta test is
**binding** (unique form, 0 free parameters, every coefficient forced and all
non-negative even integers). At $\ell=23$ the theta test is **vacuous** —
Rains–Sloane: for $N=23$ extremal lattices "almost certainly do not exist in
dimensions above 4"; Nebe: "for $N=15$ and $N=23$ there is nothing to show since
$k_{\max}(N)=0$". $D_{23}=2$ makes the modular-form bound grow like $n+2$, far
above anything geometric. So $\ell=23$ candidates "surviving" means nothing.

**Secondary:** extremal 5-modular dim 72, min 20 ($\gamma=8.944$ vs best known
8.000 = $\Gamma_{72}$, +11.8% — the mildest record among survivors).

### Worked structure for the p = 17 prescribed automorphism

Fixed-point-free $\sigma$ of prime order $p$ needs $(p-1)\mid 64$ ⟹ $p\in\{2,3,5,17\}$.
$p=17$ gives the smallest search space: $L$ is a $\mathbb Z[\zeta_{17}]$-module of rank 4.

| fact | status |
|---|---|
| $\operatorname{ord}(5)\bmod 17=16$ ⟹ (5) **inert** in $K=\mathbb Q(\zeta_{17})$, residue field $\mathbb F_{5^{16}}$ | verified |
| 5 inert in $K^+$ too, $N_{K^+/\mathbb Q}(5)=5^8$ | verified |
| $h(\mathbb Q(\zeta_{17}))=1$ ⟹ $L$ free, $L=\mathcal O_K^4$ | verified |
| $L^*/L\cong(\mathbb Z/5)^{32}$ must be an $\mathbb F_{5^{16}}$-space: $32=2\times16$ | **consistency test that could have failed — passes** |
| $\det(L_{\mathbb Z})=N_{K^+/\mathbb Q}(\det H)^2=5^{32}$ ⟹ $\det H=(\text{tot. pos. unit})\times 25$ | verified |
| Bayer-Fluckiger parity: $\mathcal O_K$-rank must be **even**; 4 passes, 1 fails | verified |

The parity result explains BF's Theorem 1 (no modular ideal lattice over
$\mathbb Q(\zeta_{17})$, since $17\equiv1\bmod4$): rank 1 is odd.

---

## 3. Verified results (reuse, don't recompute)

**Extremal theta series, rank 84, 3-modular** (acceptance test):
$$\Theta=1+119\,421\,161\,244\,q^8+12\,343\,472\,418\,816\,q^9+958\,224\,043\,687\,680\,q^{10}+\cdots$$
(exponent = norm/2). Cross-validated by two independent generator pairs.
$\theta(K_{12})$ verified by Fincke–Pohst: `[1,0,756,4032,20412,60480,139860,326592,652428,1020096,2000376]`.

**Level death points** (first negative extremal-theta coefficient):
$\ell=23\to\dim 6$, $\ell=11\to16$, $\ell=7\to30$, $\ell=5\to88$, $\ell=3>84$, $\ell=2>160$.

**LP density bounds** from Cohn's table (https://cohn.mit.edu/sphere-packing/,
Afkhami-Jeddi–Cohn–Hartman–de Laat–Tajdini 2020), converted via
$\gamma_n=4(\Delta/V_n)^{2/n}$:

| n | 44 | 48 | 56 | 64 | 72 |
|---|---|---|---|---|---|
| $\gamma_n\le$ | 6.3347 | 6.7911 | 7.6969 | 8.5948 | 9.4864 |
| $\sqrt n$ | 6.6332 | 6.9282 | 7.4833 | 8.0000 | 8.4853 |
| certified? | yes | yes | **no** | no | no |

Margin closes at −0.0394, −0.0401, −0.0407, −0.0413 per dimension; zero at
$n\approx51.4$ — independently reproducing Pavlo's "49–51 borderline". **Cohn–Elkies
published values stop at dim 36; there is nothing to look up at 52 or 56.**

**Eliminated by exact computation:**
- $(n{=}52,\ \ell{=}11,\ \min 24)$ — exact rational 2D vertex enumeration, theta
  infeasible, obstruction pinned to $q^{12},q^{13},q^{14}$
- $(n{=}52,\ \ell{=}23,\ \min 36)$ — infeasible under slack maximisation
- $(n{=}56,\ \ell{=}7,\ \min 20)$ — exact solve, negative coefficient at $q^{11}$
- $(n{=}56,\ \ell{=}11,\ \min 26)$ — $\gamma=7.839>7.697$, killed by LP bound
- $(n{=}60,\ \ell{=}7)$, $(n{=}64,\ \ell{=}7)$ — exact, negative at $q^{12}$

**Integrality check**: theta coefficients must be non-negative **even integers**.
Common denominator = 1 in all multi-parameter families ⟹ no extra obstruction.
Not decisive, but ruled out.

---

## 4. PROVEN DEAD ENDS — do not revisit

1. **Hermitian tensor products.** $K_{12}\otimes_{\mathbb Z[\omega]}G_2(3)$ gives
   rank 84, 3-modular, **min 8 not 16**. Coulangeon bound at $r=1$ gives exactly 8;
   higher-rank bounds clear 10 and 12, so the obstruction is entirely the **split
   vectors**. Leech$\otimes$rank-4 gives rank 96, det 1, split vectors at norm 4.
   Both targets: correct genus for free, factor of 2 below extremal.
2. **Cyclotomic ideal lattices at rank 96.** Levels 1 and 2 both capped. Run data
   (level 2, $\varphi(m)=96$): min 2–8, best $\gamma=5.66$, need 9.90. Governed by
   root discriminant: $\gamma\ge n/\mathrm{rd}$ and $\mathrm{rd}\in[58.6,92.5]$ for all 17 such
   fields. **Smaller rd ⟹ larger min**, confirmed on 7 data points.
3. **Genus enumeration.** Minkowski–Siegel mass $\approx10^{347}$ at rank 64
   ($10^{7.6}$ at 32, $10^{121}$ at 48). Kneser/Schiemann neighbours are hopeless.
   *This kills the previously stated research plan item "Kneser-neighbour search
   in the rank-84 3-modular genus".*
4. **Construction A from self-dual codes.** Provably never works: $\min=\min(d,p^2)/p\le p$,
   and $\min^2>n$ forces both $p>\sqrt n$ and $p<\sqrt n$. Barnes–Wall sits at exactly
   $\gamma(BW_n)^2=n/2$ for every $n$.
5. **Rank-1 / ideal lattice over $\mathbb Q(\zeta_{17})$** — BF Thm 1, $17\equiv1\bmod4$.
6. **Kirschmer's `forms/` tables** — class number ≤ 2 only, and "definite hermitian
   lattices of class number at most 2 exist only up to rank 9". Does not contain the target.

**Works but huge:** class-field-tower construction. $K$ totally imaginary of degree
$n$, $L=(\mathcal O_K,\operatorname{Tr}(x\bar y))$; then $\min(L)=n$ exactly,
$\min(L^\#)\ge n/\mathrm{rd}^2$, so $n>\mathrm{rd}^2$ suffices — **no shortest-vector
computation needed**. With Hajir–Maire ($\mathrm{rd}<83.9$): rank ≈ **7040**. This is
the explicit form of Steinberg; cost of explicitness $(83.9/2\pi e)^2\approx24\times$.
Worth writing up. Prior-art check needed (Litsyn–Tsfasman, Bergé–Martinet 1989).

---

## 5. ERROR LOG — read before writing any code

Every one of these produced a *plausible wrong answer* that was caught only by an
independent check. Budget for validation harnesses.

| # | error | symptom | fix |
|---|---|---|---|
| 1 | **Float LP with huge dynamic range** | matrix entries spanning $3.3\times10^{27}$; HiGHS returned *spurious infeasibility for every candidate*, including one previously verified feasible | exact `Fraction` arithmetic; or row-normalise **and** re-verify exactly. One verdict flipped dead→alive after normalisation |
| 2 | **Monomial conversion of Laguerre polys** | degree 13, $\alpha=15$ ⟹ catastrophic cancellation, roots meaningless | stay in the Laguerre basis |
| 3 | **Grid verification with global tolerance** | $Q(t)\ge-\epsilon\max\lvert Q\rvert$ is blind near $t=0$ where $\lvert Q\rvert$ is $10^{20}$ smaller | tolerance scaled to *local* conditioning $S(t)=\sum_k\lvert c_k\rvert\lvert L_k(t)\rvert$ |
| 4 | **Wrong Fricke eigenvalue** | $\eta(z)^6\eta(3z)^6$ has eigenvalue $(-i)^6=-1$, so it is **not** in $M(\Gamma_0(3)^+)$ | eigenvalue is $(-i)^k$; only $+1$ forms qualify. $\eta^4\eta(5z)^4$ ✓, $\eta^8\eta(2z)^8$ ✓, $\eta^6\eta(3z)^6$ ✗ |
| 5 | **Hardcoded lattice data from memory** | $\theta(K_{12})$ at $q^{10}$: used 1 693 512, truth is **2 000 376** | enumerate (Fincke–Pohst) or derive; never recall |
| 6 | **Trace convention** | plain $\operatorname{Tr}$ gives $\det=2^{2m}$, not unimodular | use $\tfrac12\operatorname{Tr}$ (Coulangeon–Nebe §3.1) |
| 7 | **Wrong ideal equation** | 2-modular: splitting $\mathfrak d=P^a\mathfrak d_{\text{odd}}$ is wrong; the condition is $I\bar I\mathfrak d=(1+\zeta_4)$. Symptom was a non-integral Gram | assert integrality rather than dropping it |
| 8 | **`Minimum(L)` in dim 96** | exact SVP, runs forever | reduced-basis upper bound first; `ShortVectors(L,bound)` only for the positive direction |

### Mandatory validation harnesses

- **Theta rings**: validate every $\Gamma_0(\ell)^+$ ring against known data before
  use. $\theta(D_4)$ kissing **24**; Barnes–Wall $\Lambda_{16}$ kissing **4320**;
  extremal minima must equal $2+2\lfloor n/m_\ell\rfloor$, $m_\ell=48/(\ell+1)$;
  known level deaths must be reproduced.
- **LP/density bounds**: any upper bound on $\gamma_n$ **below** a certain lower
  bound is invalid. Checkpoints: $\gamma_8=2$ (sharp), $\gamma_{24}=4$ (sharp),
  $\gamma_{32}\ge4$, $\gamma_{48}\ge6$ ($P_{48}$), $\gamma_{72}\ge8$ ($\Gamma_{72}$).
  My first three LP implementations each failed at least one of these.
- **Any scan where everything dies** is a bug until proven otherwise.

---

## 6. Tooling

- **Magma package** (the right one): `lat.zip` from
  https://www.math.uni-bielefeld.de/~mkirschm/magma/lat.html — `AttachSpec("lat.spec")`.
  Intrinsics: `HermitianFormWithInvariants(E,Dim,Finite,N)`, `HermitianLattice`,
  `GenusRepresentatives`, `IsIsometric`, `AutomorphismGroup`, **`Mass`**.
  *Not* shipped with Magma.
- **`Mass(L)` must gate `GenusRepresentatives(L)`.** Cheap vs. non-terminating.
  Their documented example is $\mathbb Q(\zeta_5)/\mathbb Q(\sqrt5)$ (deg 4 over 2);
  ours is deg 16 over 8. No evidence it terminates.
- Nebe's catalogue http://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/ ;
  `ExtLat.html` stops at dim 48. Hentschel ϑ-lattice tables (RWTH) pulled already.
- Speed-up: $\sigma$ of order 17 acts freely off 0 ⟹ minimal vectors in orbits of
  17; enumerate orbit representatives (Stehlé–Watkins trick, dim 80).

## 7. Attribution discipline

Coulangeon is already aware of the rank-84 and rank-96 targets — any contribution
must be genuinely new. **Prior-art check precedes every novelty claim.** Key prior
art: Kitaoka (E-type, rank ≤ 43); Coulangeon Acta Arith. 92 (2000); Coulangeon–Nebe
arXiv:1201.1832; Bayer-Fluckiger 1984 + *Cyclotomic modular lattices* (2000);
Jürgens thesis + Nebe arXiv:1910.06587 (**type** and **det-type** of prime-order
automorphisms — the established framework; the det-type refinement kills cases
the crude type analysis leaves open, and has **not** been applied to $(\ell=5,\dim 64)$);
Kirschmer JNT 197 (2019) *Automorphisms of even unimodular lattices over number fields*.

## 8. First actions for the Claude Code session

1. `Mass(L)` for the rank-4 Hermitian genus over $\mathbb Q(\zeta_{17})$ with
   $\det H=\text{unit}\times25$. **Decides whether the whole $p=17$ route is viable.**
2. If small: `GenusRepresentatives` → trace form → cheap kill (any vector of norm
   ≤ 16) → kissing number vs 21 281 280 000.
3. If large: apply **det-type** (Nebe 2021) to $(\ell=5,\dim 64)$ before anything
   else — it may eliminate $p=17$ outright, and that is a cheap result either way.
4. Independently: check Jürgens' thesis tables for $(\ell=5,\dim 64)$.
5. Parallel track: write up the rank-7040 tower construction (Section 4) after a
   prior-art check — it is a complete, correct, explicit result as it stands.

### Files produced in the chat session
`smallest_witness_final.md`, `explicit_nonEtype.md`, `prior_art_check.md`,
`positivity_all_candidates.md`, `rank64_search_strategy.md`, `p17_structure.md`,
`rank84_dr_analysis.md`, `scan_exact.py`, `ell5_extremal.py`, `ell2_extremal.py`,
`k12_theta.py`, `ce_lp4.py` (sound but 32% loose), `ext5mod64_v2.m`,
`rank96_2modular_v3.m`, Hentschel lattice tables.
