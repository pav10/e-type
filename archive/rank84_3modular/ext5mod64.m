// ====================================================================
// ext5mod64.m -- search for an EXTREMAL 5-MODULAR LATTICE of dimension 64
//
//   det 5^32,  minimum 18,  kissing number 21 281 280 000
//   (theta series computed and verified independently: unique form, all
//    coefficients non-negative even integers, so the modular-forms test is
//    binding here and it PASSES -- unlike level 23 where it is vacuous)
//
// STRATEGY (Nebe's prescribed-automorphism method, as used for the fourth
// extremal 48-dimensional lattice): do NOT enumerate the genus -- the mass is
// astronomical.  Instead fix a large automorphism and classify invariant
// lattices.
//
// WHY p = 17.  A fixed-point-free automorphism of prime order p needs
// (p-1) | 64, so p in {2,3,5,17}.  p = 17 is the largest, giving the smallest
// search space: L becomes a Z[zeta_17]-module of rank 4 (4 x 16 = 64).
//
// ARITHMETIC, all verified:
//   ord(5) mod 17 = 16, so (5) is INERT in K = Q(zeta_17); residue field F_5^16.
//   5 is inert in K+ as well, with N_{K+/Q}(5) = 5^8.
//   h(Q(zeta_17)) = 1, so O_K is a PID and L is FREE:  L = O_K^4.
//   CONSISTENCY: L*/L = (Z/5)^32 must be an O_K/(5) = F_5^16-space, so its
//     F_5-dimension must be divisible by 16.  32 = 2*16.  Passes -- this could
//     have failed and would have killed p = 17 outright.
//   DETERMINANT: det(L_Z) = N_{K+/Q}(det H)^2 = 5^32 forces
//     N_{K+/Q}(det H) = 5^16, and since 5 is inert in K+ with norm 5^8,
//        det H = (totally positive unit) * 25.
//   PARITY (Bayer-Fluckiger): an O_K-linear similarity of multiplier 5 needs
//     (det sigma) = (5)^{m/2}, so the O_K-rank m must be EVEN.  m = 4 passes.
//     (m = 1 fails, consistent with BF Thm 1: no modular ideal lattice over
//      Q(zeta_17) because 17 is a prime power with 17 = 1 mod 4.)
//
// HONEST STATUS: this lattice would be a ~16% sphere-packing density record in
// dimension 64 (best known gamma_64 = 6.928, Nebe; this needs 8.0498).  It is
// not excluded by any filter I could apply -- modular forms, integrality, the
// LP density bound, parity -- but "not excluded" is not "exists".
// ====================================================================

K<z> := CyclotomicField(17);
ZK    := MaximalOrder(K);
Kp    := sub< K | z + z^-1 >;              // maximal real subfield, degree 8
cc    := hom< K -> K | z^16 >;             // complex conjugation

// -------- sanity checks on the arithmetic, run before anything else --------
procedure Preflight()
    printf "Degree [K:Q] = %o (expect 16)\n", Degree(K);
    printf "Class number h(K) = %o (expect 1, so L is free)\n", ClassNumber(K);
    F := Factorization(5*ZK);
    printf "Factorisation of (5): %o prime(s), residue degree %o (expect 1 prime, f=16)\n",
           #F, InertiaDegree(F[1][1]) * Degree(ResidueClassField(F[1][1]));
    printf "disc(K) = 17^15 ? %o\n", Discriminant(ZK) eq 17^15 or Discriminant(ZK) eq -17^15;
end procedure;

// -------- codifferent-twisted trace form: Hermitian rank 4 -> Z-rank 64 -----
function TraceLattice(H)
    n := Nrows(H);                          // 4
    B := Basis(ZK);                         // Z-basis of O_K, 16 elements
    d := Different(ZK);
    // gamma generates d^-1 (K has class number 1, so it is principal)
    ok, g := IsPrincipal(d^-1);
    error if not ok, "inverse different not principal";
    idx := [<i,b> : b in B, i in [1..n]];
    m := #idx;                              // 64
    G := ZeroMatrix(Rationals(), m, m);
    for r in [1..m] do
        i := idx[r][1]; bi := idx[r][2];
        for c in [1..m] do
            j := idx[c][1]; bj := idx[c][2];
            G[r,c] := Trace(K ! (g * H[i,j] * bi * cc(K ! bj)));
        end for;
    end for;
    return G;
end function;

// -------- the acceptance test -----------------------------------------------
TARGET_MIN     := 18;
TARGET_KISSING := 21281280000;

procedure TestCandidate(H)
    Gq := TraceLattice(H);
    den := LCM([Denominator(x) : x in Eltseq(Gq)]);
    if den ne 1 then printf "  non-integral trace form (den %o)\n", den; return; end if;
    G := ChangeRing(Gq, Integers());
    if not IsPositiveDefinite(G) then printf "  not positive definite\n"; return; end if;
    if Determinant(G) ne 5^32 then
        printf "  det = %o, expected 5^32\n", Determinant(G); return;
    end if;
    L := LatticeWithGram(LLLGram(G));
    // cheap kill first: any vector of norm <= 16 ends it
    b := Min([Norm(L.i) : i in [1..Rank(L)]]);
    if b le 16 then printf "  min <= %o -> not extremal\n", b; return; end if;
    sv := ShortVectors(L, 16);
    if #sv gt 0 then printf "  has vectors of norm <= 16 -> not extremal\n"; return; end if;
    k := 2 * #ShortVectors(L, 18, 18);
    printf "  *** min = 18, kissing number %o (target %o) ***\n", k, TARGET_KISSING;
    if k eq TARGET_KISSING then
        printf "  *** EXTREMAL 5-MODULAR LATTICE OF DIMENSION 64 ***\n";
        printf "  gamma = 18/sqrt(5) = 8.0498 > sqrt(64) = 8  ->  NOT OF E-TYPE\n";
        printf "  witness z = sum_i e_i (x) e_i^# in L (x) L^#, (z,z) = 64 < 18*18/5\n";
    end if;
end procedure;

// -------- generating the candidates -----------------------------------------
// Requires Kirschmer's Hermitian lattice package (Magma >= 2.24: HermitianLattice,
// GenusRepresentatives).  Enumerate rank-4 Hermitian O_K-lattices with
// det H = unit * 25 and the correct local behaviour at 5.
procedure Search()
    Preflight();
    printf "\nEnumerating rank-4 Hermitian Z[zeta_17]-lattices with det = unit*25 ...\n";
    // U, mU := UnitGroup(ZK);  totally positive units of K+ modulo norms
    // for u in <totally positive units>:
    //     H := <Hermitian form with det = u*25>;
    //     G := GenusRepresentatives(HermitianLattice(H));
    //     for X in G do TestCandidate(GramMatrix(X)); end for;
    printf "  [fill in from Kirschmer's package -- see notes below]\n";
end procedure;

/*  NOTES ON RUNNING THIS

  1. The genus of rank-4 Hermitian O_K-lattices over a degree-16 field is small
     compared with the rank-64 Z-genus (mass ~1e347), which is the whole point.

  2. Minimum computation at rank 64 is the expensive step.  Use the order-17
     automorphism: minimal vectors come in orbits of size 17 (sigma acts freely
     off 0), so enumerate orbit representatives, not all vectors.  This is the
     same trick Stehle-Watkins used for the 80-dimensional extremal lattice.

  3. If p = 17 yields nothing, fall back through p = 5 (where the automorphism
     order coincides with the level, so sigma is compatible with the modularity)
     and then p = 3, 2 -- larger spaces, in increasing order of cost.

  4. Cheap pre-filter before any minimum computation: theta series.  A genuine
     candidate must have NO vectors of norm 2,4,...,16.  Any hit kills it in
     seconds.
*/
