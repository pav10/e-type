# The smallest possible non-E-type witness: current state

## Criterion

$L$ is not of E-type if $\min(L)\min(L^\#) > n$, witnessed explicitly by
$z=\sum_i e_i\otimes e_i^\#\in L\otimes L^\#$ with $(z,z)=n$. For an
$\ell$-modular lattice of rank $n$ this reads $\min^2>\ell n$, so the required
minimum is $R=$ least even integer with $R^2>\ell n$.

Three filters, all necessary:
1. **Modular forms:** $R\le E=2+2\lfloor n/m_\ell\rfloor$, $m_\ell=48/(\ell+1)$.
2. **LP density bound** (Cohn's table, Afkhami-Jeddi–Cohn–Hartman–de Laat–Tajdini):
   $\gamma=R/\sqrt\ell$ must not exceed the tabulated bound on $\gamma_n$.
3. **Theta non-negativity:** some form of weight $n/2$ for $\Gamma_0(\ell)+$ with
   constant term 1, vanishing $q^1\ldots q^{R/2-1}$, and **all** remaining
   coefficients $\ge 0$.

## Everything below rank 56 is eliminated

| n | why |
|---|---|
| ≤ 48 | LP density bound certifies $\gamma_n<\sqrt n$ (margin +0.137 at 48, closing ~0.041/dim) |
| 49–51 | same, by extrapolation; crossover at $n\approx 51.4$ |
| **52, ℓ=11, min 24** | **exact rational proof of theta-infeasibility**; obstruction pinned to $q^{12},q^{13},q^{14}$ |
| 52, ℓ=23, min 36 | theta-infeasible (normalised LP, 9 free params) — *not yet verified exactly* |
| 56, ℓ=7, min 20 | extremal case; exact solve gives a **negative coefficient at $q^{11}$** |
| 56, ℓ=11, min 26 | $\gamma=7.839$ exceeds the LP bound $7.697$ |

## The smallest survivor

> ### A 23-modular lattice of rank 56 with minimum 36
> $\det=23^{28}$, $\gamma=\gamma'=36/\sqrt{23}=7.5065$
> against $\sqrt{56}=7.4833$ (below) and the LP bound $7.6969$ (above).

The theta space has dimension 29 with 11 free parameters after imposing
vanishing through $q^{17}$, and an **exact rational certificate** (denominator
$\le 10^3$) exhibits a form with all coefficients $\ge 0$ through $q^{39}$.

Other survivors, all with exact certificates: (60, ℓ=11, min 26), (60, ℓ=23, min 38),
(64, ℓ=5, min 18) [exact unique solve, all coefficients positive, kissing number
21 281 280 000], (64, ℓ=23, min 40), (68, ℓ=7/11/23), (72, ℓ=5/11/23).

Also proven dead: (60, ℓ=7, min 22) and (64, ℓ=7, min 22), both by exact solve
with a negative coefficient at $q^{12}$; (64, ℓ=11, min 28) and (72, ℓ=7, min 24)
by LP.

## Methodological warning, recorded because it nearly cost a false result

The first pass of this scan returned **"DEAD" for every candidate**, including
(64, ℓ=5, min 18) — a case verified earlier in the same session as feasible with
all-positive coefficients. The cause: the LP matrix entries span a ratio of
$3.3\times10^{27}$, far outside double precision, and HiGHS returned spurious
infeasibility throughout. After row normalisation the verdict for
(56, ℓ=23, min 36) **flipped from dead to alive**.

Consequences for how to read the table: verdicts marked *exact* (unique-form
solves, and the 2-parameter vertex enumeration at $n=52$) are rigorous. Verdicts
resting on the normalised float LP have already been wrong once. In particular
**(52, ℓ=23, min 36) should be redone in exact arithmetic before claiming 56 is
minimal** — it is the single load-bearing unverified cell.

## Where this leaves the search

The target has moved from rank 84/96 (where the session began) to **rank 56**,
and the object is concrete: a 23-modular lattice of rank 56, determinant
$23^{28}$, minimum 36. All three filters pass. Nothing about existence follows —
these are necessary conditions only — but no smaller candidate survives on
current knowledge.
