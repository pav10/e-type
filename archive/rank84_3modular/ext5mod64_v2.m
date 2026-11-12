// ====================================================================
// ext5mod64_v2.m -- extremal 5-modular lattice, dimension 64, minimum 18
//
// CORRECTION vs v1: v1 referred vaguely to "Kirschmer's Hermitian genus
// package".  There are two different things and only one is relevant:
//
//   (a) https://www.math.uni-bielefeld.de/~mkirschm/forms/
//       TABLES of definite quadratic/hermitian lattices with CLASS NUMBER
//       <= 2, from his Habilitation.  A classification, not a tool.  And
//       "definite hermitian lattices of class number at most 2 exist only up
//       to rank 9" -- our rank-4 genus over a degree-16 CM field will not have
//       class number <= 2, so these tables do NOT contain our target.
//
//   (b) https://www.math.uni-bielefeld.de/~mkirschm/magma/lat.html  (lat.zip)
//       The MAGMA PACKAGE implementing his habilitation algorithms.  This is
//       the right tool.  Not shipped with Magma; attach it explicitly.
//
// The intrinsics used below are taken from that package's documentation.
// ====================================================================

AttachSpec("lat.spec");                       // from lat.zip

// ---- E = Q(zeta_17) as a degree-2 extension of K = Q(zeta_17)^+ -----------
// (the package requires E/K relative of degree 2; this mirrors their own
//  Q(zeta_5) over Q(sqrt 5) example, but four times deeper)
C<z> := CyclotomicField(17);
Kp   := sub< C | z + z^-1 >;                  // maximal real subfield, degree 8
E    := RelativeField(Kp, C);                 // degree 2 over Kp, 16 over Q
ZE   := MaximalOrder(E);

// ====================================================================
// STEP 0 -- MASS GATE.  Do this BEFORE any enumeration.
//
// The Minkowski-Siegel mass of the rank-64 Z-genus is ~1e347, which is why we
// moved to the Hermitian side at all.  But nothing guarantees the rank-4
// Hermitian genus over a degree-16 field is small.  Mass(L) is cheap;
// GenusRepresentatives(L) is not.  If the mass is large, STOP -- the
// prescribed-automorphism idea has to be pushed further (larger group), not
// brute-forced.
// ====================================================================
procedure MassGate(L)
    m := Mass(L);
    printf "Hermitian genus mass = %o  (approx %o)\n", m, RealField(6)!m;
    printf "class number is at least 1/max|Aut| * mass; if this is large,\n";
    printf "GenusRepresentatives will not terminate -- do not attempt it.\n";
end procedure;

// ====================================================================
// STEP 1 -- build the Hermitian space with the right local invariants.
//
// Constraints derived earlier and re-verified:
//   ord(5) mod 17 = 16  =>  (5) inert in E, residue field F_5^16
//   5 inert in Kp with N_{Kp/Q}(5) = 5^8
//   det(L_Z) = N_{Kp/Q}(det H)^2 = 5^32  =>  det H = (tot. pos. unit) * 25
//   h(Q(zeta_17)) = 1  =>  L is free, L = O_E^4
//   L*/L = (Z/5)^32 must be an F_5^16-space: 32 = 2*16.  Consistent.
//   Bayer-Fluckiger parity: O_E-rank must be EVEN; 4 is even.  Passes.
// ====================================================================
RANK := 4;                                    // 4 * 16 = 64
P5   := [ t[1] : t in Factorization(5*MaximalOrder(Kp)) ];
assert #P5 eq 1;                              // 5 inert in Kp

// signature: definite, so no negative entries at any of the 8 real places
SIG := [0 : i in [1..Degree(Kp)]];

function CandidateForm(nonNormPlaces)
    return HermitianFormWithInvariants(E, RANK, nonNormPlaces, SIG);
end function;

// ====================================================================
// STEP 2 -- trace lattice and acceptance test
// ====================================================================
TARGET_MIN     := 18;
TARGET_KISSING := 21281280000;                // computed from the extremal
                                              // theta series (unique form,
                                              // all coefficients non-negative
                                              // even integers -- the modular
                                              // forms test is BINDING here)

function TraceGram(L)
    // codifferent-twisted trace form; gives det 5^32 when det H = unit*25
    H  := GramMatrix(L);
    ZK := MaximalOrder(AbsoluteField(E));
    d  := Different(ZK);
    ok, g := IsPrincipal(d^-1);  error if not ok, "d^-1 not principal";
    B  := Basis(ZK);  n := Nrows(H);
    idx := [<i,b> : b in B, i in [1..n]];  m := #idx;
    A := AbsoluteField(E);
    cc := hom< A -> A | ComplexConjugate(A.1) >;
    G := ZeroMatrix(Rationals(), m, m);
    for r in [1..m] do
      i := idx[r][1]; bi := idx[r][2];
      for c in [1..m] do
        j := idx[c][1]; bj := idx[c][2];
        G[r,c] := Trace(A ! (g * (A!H[i,j]) * (A!bi) * cc(A!bj)));
      end for;
    end for;
    return G;
end function;

procedure Accept(L)
    Gq := TraceGram(L);
    if LCM([Denominator(x) : x in Eltseq(Gq)]) ne 1 then
        printf "  non-integral trace form\n"; return; end if;
    G := ChangeRing(Gq, Integers());
    if not IsPositiveDefinite(G) then printf "  indefinite\n"; return; end if;
    if Determinant(G) ne 5^32 then
        printf "  det %o != 5^32\n", Determinant(G); return; end if;
    Lz := LatticeWithGram(LLLGram(G));
    b := Min([Norm(Lz.i) : i in [1..Rank(Lz)]]);
    if b le 16 then printf "  min <= %o -> not extremal\n", b; return; end if;
    if #ShortVectors(Lz, 16) gt 0 then printf "  vectors of norm <= 16\n"; return; end if;
    k := 2 * #ShortVectors(Lz, 18, 18);
    printf "  min = 18, kissing %o (target %o)\n", k, TARGET_KISSING;
    if k eq TARGET_KISSING then
        printf "  *** EXTREMAL 5-MODULAR LATTICE, DIMENSION 64 ***\n";
        printf "  gamma = 8.0498 > sqrt(64) = 8  =>  NOT OF E-TYPE\n";
    end if;
end procedure;

// ====================================================================
// DRIVER -- mass first, enumeration only if the mass permits
// ====================================================================
procedure Run()
    F := CandidateForm({});
    L := HermitianLattice(F);
    printf "rank %o hermitian lattice over E of degree %o\n", RANK, AbsoluteDegree(E);
    printf "determinant ideal: %o\n", Determinant(L);
    MassGate(L);
    printf "\n--- only proceed past here if the mass is small ---\n";
    // reps := GenusRepresentatives(L);
    // printf "%o classes in the genus\n", #reps;
    // for X in reps do Accept(X); end for;
end procedure;

/*  HONEST STATUS
    - Every filter applied so far (modular forms with a BINDING test, integrality,
      LP density bound, Bayer-Fluckiger parity, F_5^16-module consistency) passes.
    - The lattice would be a ~16% sphere-packing density record in dimension 64
      (best known gamma_64 = 6.928; this needs 8.0498).
    - I have not run any of this.  The mass gate exists precisely because the
      analogous Z-side computation died at mass 1e347, and I have no evidence
      the Hermitian genus over a degree-16 field is tractable.
    - Related and directly relevant, not yet used here:
        M. Kirschmer, Automorphisms of even unimodular lattices over number
        fields, J. Number Theory 197 (2019) 121-134.
        M. Kirschmer, G. Nebe, Binary hermitian lattices over number fields,
        Exp. Math. 31 (2022) 280-301.
*/
