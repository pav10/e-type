// ====================================================================
// rank96_2modular_v3.m
//
// Fix vs v2:  ZBasis is for orders in associative algebras, not RngOrdIdl.
//             For ideals of a number-field order the Z-basis is Basis(I).
//             Also: conjugation is now only ever applied to INTEGRAL ideals.
//
// Why that is safe: T := D/(lam) is integral, because at each prime above 2
//   v_P(D) - v_P(lam) = (a-1)2^(a-1) - 2^(a-2) = 2^(a-2)(2a-3) >= 0   (2^a || m),
// and lam is a unit at every odd prime.  Hence A is integral too, and the
// defining equation can be checked as  T = A * Abar  entirely over integral
// ideals.  The fractional ideal I = A^-1 is formed only to extract a basis.
//
// Solvability: fails exactly when a = 2 (m = 4*odd), i.e. m = 260, 420.
// ====================================================================

CAND2 := [208, 224, 280, 288, 312, 336, 360];

function Cyclotomic2Modular(m)
    K<z> := CyclotomicField(m);
    ZK   := MaximalOrder(K);
    n    := Degree(K);
    error if n ne 96,        "phi(m) must be 96";
    error if m mod 4 ne 0,   "need 4 | m";

    cc  := hom< K -> K | z^(m-1) >;          // zeta -> zeta^-1
    i4  := z^(m div 4);
    lam := 1 + i4;
    assert lam*cc(lam) eq 2;                 // lam * lambar = 2

    // conjugate of an INTEGRAL ideal (Basis, not ZBasis)
    conj := function(J)
        return ideal< ZK | [ ZK ! cc(K ! b) : b in Basis(J) ] >;
    end function;

    D   := Different(ZK);
    Lam := ideal< ZK | ZK ! lam >;
    T   := D / Lam;
    assert IsIntegral(T);                    // guaranteed by the valuation count

    // ---- write T = A * Abar ---------------------------------------
    A    := 1*ZK;
    seen := {};
    for t in Factorization(T) do
        Q, e := Explode(t);
        if Q in seen then continue; end if;
        Qb := conj(Q);
        if Qb eq Q then
            if IsOdd(e) then
                printf "  m=%o: self-conjugate prime (norm %o) has odd valuation %o -> none\n",
                       m, Norm(Q), e;
                return false, _;
            end if;
            A := A * Q^(e div 2);
        else
            A := A * Q^e;
            Include(~seen, Qb);
        end if;
        Include(~seen, Q);
    end for;

    // ---- verify entirely over integral ideals ----------------------
    if T ne A * conj(A) then
        printf "  m=%o: T = A*Abar verification FAILED\n", m;
        return false, _;
    end if;

    // ---- build the trace lattice on I = A^-1 ------------------------
    I  := A^(-1);                            // fractional
    B  := Basis(I);                          // Z-basis, elements of K
    assert #B eq n;
    Gq := Matrix(Rationals(), n, n,
          [ Trace(K!B[i] * cc(K!B[j])) : i, j in [1..n] ]);
    assert &and[ Denominator(Gq[i,j]) eq 1 : i,j in [1..n] ];   // now guaranteed
    G := ChangeRing(Gq, Integers());
    if not IsPositiveDefinite(G) then
        printf "  m=%o: trace form indefinite -> needs a totally positive unit twist\n", m;
        return false, _;
    end if;
    return true, LatticeWithGram(G);
end function;

procedure TestExtremal2(m)
    ok, L := Cyclotomic2Modular(m);
    if not ok then return; end if;
    d := Determinant(L);
    printf "=== m = %o : rank %o, det = %o (2^%o), even = %o ===\n",
           m, Rank(L), d, Ilog(2,d), IsEven(L);
    assert d eq 2^48;
    L := LLL(L);
    try L := BKZ(L, 40); catch e ; end try;
    printf "  minimum = %o   (extremal would be 14)\n", Minimum(L);
    if Minimum(L) eq 14 then
        printf "  *** EXTREMAL 2-MODULAR rank 96 ***\n";
        printf "  gamma'^2 = 14*7 = 98 > 96 -> NOT OF E-TYPE\n";
        printf "  witness z = sum_i e_i (x) e_i^#, (z,z) = 96 < 98\n";
    end if;
end procedure;

/*  load "rank96_2modular_v3.m";
    for m in CAND2 do TestExtremal2(m); end for;                        */
