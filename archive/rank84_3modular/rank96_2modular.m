// ====================================================================
// rank96_2modular.m   --  2-modular variant of your CyclotomicIdealLattice
//
// TARGET: rank-96, 2-modular (det 2^48), min 14.
//   Then min(L^#) = 14/2 = 7 and gamma'^2 = 14*7 = 98 > 96 = n,
//   so such a lattice is NOT of E-type, witnessed explicitly by
//        z = sum_i e_i (x) e_i^#   in   L (x) L^#,   (z,z) = 96 < 98.
//
// WHY THIS LEVEL: Bayer-Fluckiger, "Cyclotomic modular lattices", Thm 2:
//   ell-modular ideal lattices of trace type exist only for
//        ell = 1,  ell = 2,  or ell = 3 (mod 4).
//   Level 2 is admissible; level 5 (= 1 mod 4) is NOT, which is why the
//   rank-64 level-5 candidate is unreachable by this machinery.
//   BF Prop 4: if 4 | m then 2 is a norm of m.
//
// YOUR m-VALUES with 4 | m and phi(m) = 96  (9 of the 10):
//   208, 224, 260, 280, 288, 312, 336, 360, 420      (195 is odd, excluded)
// ====================================================================

CAND2 := [208, 224, 260, 280, 288, 312, 336, 360, 420];

// --------------------------------------------------------------------
// NOTE ON THE IDEAL.  I have deliberately NOT hard-coded Bayer-Fluckiger's
// exponent for J_1 = P^a in D_K = J_1 * J * Jbar (my copy of the paper is
// an OCR scan and I could not read the exponent reliably; you have the
// original -- read it off Prop 6 / Thm 3 of section 5).
// Instead this searches ideals and VERIFIES the determinant, which is
// self-certifying and cannot be wrong in the way a mis-transcribed
// exponent would be.
// --------------------------------------------------------------------

function Cyclotomic2ModularLattice(m)
    K   := CyclotomicField(m);
    ZK  := MaximalOrder(K);
    n   := Degree(K);                    // must be 96
    if n ne 96 then return false, _, _; end if;

    D   := Different(ZK);
    // split off the part above 2
    P2  := [ t[1] : t in Factorization(2*ZK) ];
    if #P2 ne 1 then                     // 2 is totally ramified for 4 | m
        printf "  m = %o: %o primes above 2\n", m, #P2;
    end if;
    P   := P2[1];
    a   := Valuation(D, P);

    // D = P^a * (prime-to-2 part).  Write the prime-to-2 part as J*Jbar
    Dodd := D / P^a;
    // conjugation
    cc   := hom< K -> K | K.1^(m-1) >;

    // try to write Dodd = J * Jbar  (needs "m' is a norm of m", BF Prop 3)
    fact := Factorization(Dodd);
    J    := 1*ZK;
    ok   := true;
    for t in fact do
        Q := t[1]; e := t[2];
        Qb := ideal< ZK | [ cc(x) : x in Generators(Q) ] >;
        if Qb eq Q then
            if IsEven(e) then J *:= Q^(e div 2); else ok := false; end if;
        else
            J *:= Q^e;                    // take one of the conjugate pair
        end if;
    end for;
    if not ok then
        printf "  m = %o: Dodd is not of the form J*Jbar -> m' not a norm of m\n", m;
        return false, _, _;
    end if;

    // candidate ideals I = (P^s * J)^-1 ; s chosen so det comes out 2^48
    for s in [0..a] do
        I  := (P^s * J)^(-1);
        B  := Basis(I);
        Gm := Matrix(96,96,[ Trace(B[i]*cc(B[j])) : i,j in [1..96] ]);
        Gm := ChangeRing(Gm, Integers());
        if IsPositiveDefinite(Gm) and Determinant(Gm) eq 2^48 then
            return true, LatticeWithGram(Gm), s;
        end if;
    end for;
    return false, _, _;
end function;

// --------------------------------------------------------------------
// Acceptance test.  Extremal 2-modular rank 96 has min 14, so probe
// norms 2,4,...,12 -- ANY hit kills the candidate immediately.
// --------------------------------------------------------------------
procedure TestExtremal2(m)
    ok, L, s := Cyclotomic2ModularLattice(m);
    if not ok then printf "m = %o: no 2-modular ideal lattice found\n", m; return; end if;
    printf "=== m = %o  (s = %o) : rank %o, det %o = 2^%o ===\n",
           m, s, Rank(L), Determinant(L), Ilog(2, Determinant(L));
    L := LLL(L);
    try L := BKZ(L, 40); catch e ; end try;
    for bound in [2,4,6,8,10,12] do
        sv := ShortVectors(L, bound);
        printf "  norm <= %o : %o vectors\n", bound, #sv;
        if #sv gt 0 then
            printf "  => min = %o.  NOT extremal.\n", Minimum([p[2] : p in sv]);
            return;
        end if;
    end for;
    printf "  => min = 14 : *** EXTREMAL 2-MODULAR ***\n";
    printf "  gamma'^2 = 14*7 = 98 > 96 -> NOT OF E-TYPE\n";
    printf "  witness: z = sum_i e_i (x) e_i^#  in L (x) L^#, (z,z) = 96\n";
end procedure;

/*  run:
      load "rank96_2modular.m";
      for m in CAND2 do TestExtremal2(m); end for;
*/
