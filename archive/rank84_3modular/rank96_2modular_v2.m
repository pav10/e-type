// ====================================================================
// rank96_2modular_v2.m   --  CORRECTED
//
// Fixes two bugs in v1:
//   (1) WRONG IDEAL EQUATION.  2-modularity = "mult by lam carries I* to I"
//       with lam = 1 + zeta_4  (lam * lambar = 2).  Since I* = (Ibar*D)^-1,
//       the condition is exactly
//                       I * Ibar * D  =  (lam).
//       v1 instead split D = P^a * Dodd and asked Dodd = J*Jbar, which does
//       NOT give I*Ibar*D integral -> non-integral Gram -> ChangeRing error.
//       With the correct equation, integrality is automatic and
//                 det = N_{K/Q}(1+i) = 2^(96/2) = 2^48
//       falls out with no search over exponents.
//   (2) v1 did P2[1], picking one prime above 2 arbitrarily when g > 1.
//
// SOLVABILITY: D/(lam) = A*Abar is possible iff every SELF-CONJUGATE prime
// has even valuation.  At the prime above 2, v_P(D/(lam)) = 2^(a-2)(2a-3)
// where 2^a || m -- ODD exactly when a = 2.  Hence m = 4*odd is excluded.
//   fails:        260, 420     (m = 4 * odd)
//   should work:  208, 224, 280, 288, 312, 336, 360
// ====================================================================

CAND2 := [208, 224, 280, 288, 312, 336, 360];   // 260, 420 excluded a priori

function ConjIdeal(ZK, cc, I)
    return ideal< ZK | [ cc(x) : x in ZBasis(I) ] >;
end function;

function Cyclotomic2Modular(m)
    K<z> := CyclotomicField(m);
    ZK   := MaximalOrder(K);
    n    := Degree(K);
    error if n ne 96, "phi(m) must be 96";
    error if m mod 4 ne 0, "need 4 | m";

    cc  := hom< K -> K | z^(m-1) >;         // zeta -> zeta^-1
    i4  := z^(m div 4);                     // primitive 4th root of unity
    lam := 1 + i4;                          // lam * conj(lam) = 2
    assert lam*cc(lam) eq 2;

    D   := Different(ZK);
    Lam := ideal< ZK | lam >;
    T   := D / Lam;                         // want T = A * Abar

    // ---- try to write T = A * Abar --------------------------------
    A    := 1*ZK;
    done := {};
    for t in Factorization(T) do
        Q, e := Explode(t);
        if Q in done then continue; end if;
        Qb := ConjIdeal(ZK, cc, Q);
        if Qb eq Q then
            if IsOdd(e) then
                printf "  m=%o: SELF-CONJUGATE prime of norm %o has ODD valuation %o -> no 2-modular ideal lattice\n",
                       m, Norm(Q), e;
                return false, _;
            end if;
            A := A * Q^(e div 2);
        else
            A := A * Q^e;                   // take one of each conjugate pair
            Include(~done, Qb);
        end if;
        Include(~done, Q);
    end for;

    I := A^(-1);
    // ---- verify the defining equation before building anything ----
    if I * ConjIdeal(ZK, cc, I) * D ne Lam then
        printf "  m=%o: ideal equation I*Ibar*D = (lam) FAILED\n", m;
        return false, _;
    end if;

    // ---- trace form -----------------------------------------------
    B  := ZBasis(I);
    Gq := Matrix(Rationals(), n, n,
          [ Trace(K!B[i] * cc(K!B[j])) : i, j in [1..n] ]);
    // integrality is now guaranteed; assert it rather than silently dropping
    assert &and[ Denominator(Gq[i,j]) eq 1 : i,j in [1..n] ];
    G := ChangeRing(Gq, Integers());
    if not IsPositiveDefinite(G) then
        // the trace form of I may be negative/indefinite for a bad choice of
        // generator; twist by a totally positive unit if so
        printf "  m=%o: trace form not positive definite (needs a tot. pos. twist)\n", m;
        return false, _;
    end if;
    return true, LatticeWithGram(G);
end function;

// --------------------------------------------------------------------
// Extremal 2-modular rank 96 has min 14.  Probe 2,4,...,12: any hit kills it.
// --------------------------------------------------------------------
procedure TestExtremal2(m)
    ok, L := Cyclotomic2Modular(m);
    if not ok then return; end if;
    printf "=== m = %o : rank %o, det = 2^%o, even = %o ===\n",
           m, Rank(L), Ilog(2, Determinant(L)), IsEven(L);
    assert Determinant(L) eq 2^48;
    L := LLL(L);
    try L := BKZ(L, 40); catch e ; end try;
    for bound in [2,4,6,8,10,12] do
        sv := ShortVectors(L, bound);
        if #sv gt 0 then
            printf "  min = %o  ->  not extremal (need 14)\n", Minimum(L);
            return;
        end if;
        printf "  no vectors of norm <= %o\n", bound;
    end for;
    printf "  *** min = 14 : EXTREMAL 2-MODULAR, rank 96 ***\n";
    printf "  gamma'^2 = 14*7 = 98 > 96  ->  NOT OF E-TYPE\n";
    printf "  witness z = sum_i e_i (x) e_i^# in L (x) L^#, (z,z) = 96 < 98\n";
end procedure;

/*  load "rank96_2modular_v2.m";
    for m in CAND2 do TestExtremal2(m); end for;                        */
