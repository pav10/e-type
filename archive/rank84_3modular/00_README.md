# Hentschel theta-lattice (vartheta-lattice) data — consolidated

Source: M. Hentschel, RWTH Aachen
  http://www.matha.rwth-aachen.de/de/mitarbeiter/hentschel/
A "vartheta-lattice" (a.k.a. Eisenstein lattice for d=3) is a positive
definite Hermitian O_K-lattice whose trace lattice (x.y = Tr_{K/Q} h(x,y))
is an EVEN UNIMODULAR Z-lattice. K = Q(sqrt(-d)).

## Retrieval status (network-limited; large files truncated by fetch)
| file                | K          | ranks      | #lattices | status here            |
|---------------------|------------|------------|-----------|------------------------|
| dim12Eisenstein     | Q(sqrt-3)  | 12 (Z=24)  | 5         | COMPLETE               |
| d-3 rank16          | Q(sqrt-3)  | 16 (Z=32)  | many      | partial (1 of N)       |
| d-7 rank4+8         | Q(sqrt-7)  | 4, 8       | 1 + 3     | COMPLETE               |
| d-7 rank12          | Q(sqrt-7)  | 12 (Z=24)  | ~9 + more | partial (first ~9)     |
| d-11 rank4+8        | Q(sqrt-11) | 4, 8       | 1 + 6     | COMPLETE               |
| d-19 rank4+8        | Q(sqrt-19) | 4, 8       | 2 + 83    | partial (first ~60)    |
| d-43 rank4          | Q(sqrt-43) | 4          | 4         | COMPLETE               |
| d-67 rank4          | Q(sqrt-67) | 4          | 6         | COMPLETE               |
| d-163 rank4         | Q(sqrt-163)| 4          | 16        | COMPLETE               |

For the truncated three, re-download the originals directly (they are plain
Magma source). The COMPLETE files below are ready to `load` in Magma.

## Relevance to the targets
- dim12 d=3 (Eisenstein, Z-rank 24): the FIVE even unimodular Eisenstein
  lattices — these realise 3E8, 4E6, 6D4, 12A2, Leech as Hermitian Z[w]-rank
  12. The Leech one (Gram24) is the key 3-modular tensor brick.
- d-7 rank12: the Hermitian Z[(1+sqrt-7)/2]-structures whose trace lattice is
  the Leech lattice — these are the P_i used in Nebe's Gamma_72 = P_b (x) P_i.

## VERIFIED inner-product convention (important!)
Checked here on the d=3 3E8 Gram (Gram02): with the HALF-TRACE form
    x . y = (1/2) Tr_{K/Q}( h(x,y) )
the trace lattice determinant is exactly 1 (even unimodular), diagonal
norm 2, minimum 2.  This is the Coulangeon-Nebe convention (1201.1832 sec 3.1:
"x.y = (1/2)Tr_{K/Q}(x ybar)").  If you instead use the plain trace Tr (no 1/2)
you get det = 2^{2m} and norms doubled -- a similar lattice, not unimodular.
So: ALWAYS use (1/2)Tr when forming the Z-lattice from these Hermitian Grams.
The Hermitian matrices themselves were verified Hermitian (A = conj-transpose A).

## Files
- d3_dim12_eisenstein.magma   COMPLETE (5 lattices: 3E8,4E6,6D4,12A2,Leech)
- d7_rank4and8.magma          COMPLETE (E8-form rank4 ; 3 rank-8)
- d11_rank4and8.magma         COMPLETE (1 rank4 ; 6 rank8)
- d43_67_163_rank4.magma      COMPLETE (4 + 6 + 16 rank-4 lattices)
- d7_rank12_PARTIAL.magma     PARTIAL  (Leech Hermitian structures; re-download for all)
- d7_rank12_more_raw.txt      raw reference for the extra rank-12 matrices fetched
