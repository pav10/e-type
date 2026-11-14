#!/usr/bin/env python3
"""Rebuild the SRFP311T1 rank-96 Gram matrix.  Construction stages 1-4 are copied verbatim
from their reconstruct_and_audit.py (helpers + modular linear algebra + stages 1-4);
the PARI subprocess calls, apt-get installs and later "certification" stages are omitted.
Writes external/srfp_G96.txt."""
from fractions import Fraction
import numpy as np, sympy as sp, sys, os

class CertificationError(Exception):
    """Base exception for mathematical and certification failures."""
    pass


class ArithmeticIntegrityError(CertificationError):
    """Raised when exact divisibility, symmetry, or integrality fails."""
    pass


class PARIExecutionError(CertificationError):
    """Raised when PARI/GP fails, times out, or produces malformed output."""
    pass


def strict_int(val, label=""):
    """Strictly validates and converts an input to integer without silent float truncation."""
    if isinstance(val, (int, np.integer)) and not isinstance(val, bool):
        return int(val)
    if isinstance(val, Fraction):
        if val.denominator != 1:
            raise ArithmeticIntegrityError(
                f"{label}: Expected integer, got Fraction with denominator {val.denominator}: {val}"
            )
        return int(val.numerator)
    if isinstance(val, sp.Integer):
        return int(val)
    if isinstance(val, sp.Rational):
        if val.q != 1:
            raise ArithmeticIntegrityError(
                f"{label}: Expected integer, got SymPy Rational with denominator {val.q}: {val}"
            )
        return int(val.p)
    raise ArithmeticIntegrityError(
        f"{label}: Strict integer expected, got invalid type {type(val)}: {val}"
    )


def exact_div(a, d, label=""):
    """Enforces exact integer divisibility; rejects floor division truncation."""
    a_int = strict_int(a, f"{label}.numerator")
    d_int = strict_int(d, f"{label}.denominator")
    if d_int == 0:
        raise ArithmeticIntegrityError(f"Division by zero in {label}")
    if a_int % d_int != 0:
        raise ArithmeticIntegrityError(
            f"Non-exact division in {label}: {a_int} % {d_int} = {a_int % d_int} != 0"
        )
    return a_int // d_int


def exact_div_vec(vec, d, label=""):
    return [exact_div(v, d, f"{label}[{i}]") for i, v in enumerate(vec)]


def exact_div_mat(mat, d, label=""):
    m, n = len(mat), len(mat[0])
    return [
        [exact_div(mat[r][c], d, f"{label}[{r},{c}]") for c in range(n)]
        for r in range(m)
    ]



# STEP 1: MODULAR LINEAR ALGEBRA IN PURE PYTHON INTEGERS
# ------------------------------------------------------------------------------
def rref_mod_p(A_list, p=7):
    M = [[strict_int(x) % p for x in row] for row in A_list]
    m, n = len(M), len(M[0])
    row, pivots = 0, []
    for col in range(n):
        piv = next((r for r in range(row, m) if M[r][col] % p != 0), None)
        if piv is None:
            continue
        M[row], M[piv] = M[piv], M[row]
        inv = pow(M[row][col], p - 2, p)
        M[row] = [(x * inv) % p for x in M[row]]
        for r in range(m):
            if r != row and M[r][col] % p != 0:
                factor = M[r][col]
                M[r] = [(M[r][c] - factor * M[row][c]) % p for c in range(n)]
        pivots.append(col)
        row += 1
        if row == m:
            break
    return M, pivots


def nullspace_mod_p(M_list, p=7):
    R, pivots = rref_mod_p(M_list, p)
    n = len(M_list[0])
    free_vars = [c for c in range(n) if c not in pivots]
    ns = [[0] * n for _ in range(len(free_vars))]
    for i, free in enumerate(free_vars):
        ns[i][free] = 1
        for r, piv in enumerate(pivots):
            ns[i][piv] = (-R[r][free]) % p
    return ns


def mat_mul_mod_p(A, B, p=7):
    m, k, n = len(A), len(A[0]), len(B[0])
    C = [[0] * n for _ in range(m)]
    for i in range(m):
        for j in range(n):
            s = 0
            for l in range(k):
                s = (s + strict_int(A[i][l]) * strict_int(B[l][j])) % p
            C[i][j] = s
    return C


def transpose(A):
    return [[A[r][c] for r in range(len(A))] for c in range(len(A[0]))]


def construct_deterministic_isotropic_subspace(B_K_list, p=7):
    B12 = [row[24:] for row in B_K_list[:24]]
    B12_sp = sp.Matrix(B12)
    B12_inv = B12_sp.inv_mod(p).tolist()
    J12 = [[0] * 24 for _ in range(24)]
    for i in range(12):
        J12[i][12 + i] = 1
        J12[12 + i][i] = -1
    M_symp = mat_mul_mod_p(J12, transpose(B12_inv), p=p)
    return [
        [(1 if i == j else 0) for j in range(24)] + M_symp[i] for i in range(24)
    ]


def exact_integer_quadratic_form(G_mat, v_vec):
    n = len(v_vec)
    total = 0
    for i in range(n):
        v_i = strict_int(v_vec[i])
        if v_i == 0:
            continue
        row_sum = 0
        for j in range(n):
            v_j = strict_int(v_vec[j])
            if v_j != 0:
                row_sum += strict_int(G_mat[i][j]) * v_j
        total += v_i * row_sum
    return total


# ------------------------------------------------------------------------------


def validate_integer_matrix(matrix, rows, cols, label):
    """Enforces rigorous structure, dimensions, non-emptiness, and integer type purity."""
    if not isinstance(matrix, list):
        raise ArithmeticIntegrityError(f"{label}: Expected list of lists, got {type(matrix)}")
    if len(matrix) != rows:
        raise ArithmeticIntegrityError(f"{label}: Expected {rows} rows, got {len(matrix)}")
    if rows == 0:
        raise ArithmeticIntegrityError(f"{label}: Matrix has zero rows")

    for i, row in enumerate(matrix):
        if not isinstance(row, list):
            raise ArithmeticIntegrityError(f"{label}: Row {i} is not a list ({type(row)})")
        if len(row) != cols:
            raise ArithmeticIntegrityError(f"{label}: Row {i} must contain {cols} entries, found {len(row)}")
        for j, value in enumerate(row):
            if type(value) is not int:
                raise ArithmeticIntegrityError(
                    f"{label}[{i},{j}] is not a pure Python integer: {type(value)}"
                )

    return matrix



def build():
    print("\n[*] Stage 1: Cyclotomic relative trace on Q(zeta_35)...", flush=True)
    x = sp.Symbol("x")
    phi_35 = sp.cyclotomic_poly(35, x)
    zeta7 = x**5
    g = zeta7 + zeta7**2 - zeta7**3 + zeta7**4 - zeta7**5 - zeta7**6
    g_inv = sp.rem(-g / 7, phi_35)
    denom = sp.rem(x - x**34, phi_35)
    denom_inv = sp.invert(denom, phi_35)

    def get_unit(a):
        return sp.rem(sp.rem(x**a - x ** (35 - a), phi_35) * denom_inv, phi_35)

    u = sp.rem(-1 * get_unit(4) * get_unit(8) * get_unit(16), phi_35)
    galois_F = [1, 2, 4, 8, 9, 11, 16, 18, 22, 23, 29, 32]
    y = sp.Symbol("y")
    P_y = 1
    for k in galois_F:
        P_y = sp.rem(P_y * (y - x**k), phi_35)

    poly_P = sp.Poly(P_y, y)
    coeffs_P = poly_P.all_coeffs()
    P_prime_x = sum(
        deg * c * (x ** (deg - 1))
        for deg, c in enumerate(reversed(coeffs_P))
        if deg > 0
    )
    delta = sp.rem(P_prime_x, phi_35)
    beta = sp.rem((x**30) * delta, phi_35)
    alpha0 = sp.rem(beta * g_inv, phi_35)
    alpha_pos = sp.rem(u * alpha0, phi_35)
    lam = sp.invert(alpha_pos, phi_35)

    coeffs_phi = [strict_int(c) for c in sp.Poly(phi_35, x).all_coeffs()]
    c_top = -np.array(coeffs_phi[1:], dtype=object)[::-1]
    C_comp = np.zeros((24, 24), dtype=object)
    C_comp[1:, :-1] = np.eye(23, dtype=object)
    C_comp[:, -1] = c_top

    C_pow = [np.eye(24, dtype=object)]
    for _ in range(1, 71):
        C_pow.append(C_pow[-1] @ C_comp)

    Tr_KF_mat = sum(
        np.column_stack([C_pow[(k * j) % 35][:, 0] for j in range(24)])
        for k in galois_F
    )
    conj_mat = np.column_stack([C_pow[(34 * j) % 35][:, 0] for j in range(24)])

    lam_poly = sp.Poly(lam, x)
    M_5lam = np.zeros((24, 24), dtype=object)
    for p_tuple, val in lam_poly.as_dict().items():
        if val.q not in (1, 5):
            raise ArithmeticIntegrityError(f"Unexpected denominator in lambda: {val.q}")
        M_5lam += strict_int(val * 5) * C_pow[p_tuple[0]]

    B = np.column_stack(
        [C_pow[(7 * a + 5 * b) % 35][:, 0] for a in range(4) for b in range(3)]
    )
    B_bar = conj_mat @ B

    g_vec = (
        C_pow[5][:, 0]
        + C_pow[10][:, 0]
        - C_pow[15][:, 0]
        + C_pow[20][:, 0]
        - C_pow[25][:, 0]
        - C_pow[30][:, 0]
    )
    idx_nz = int(np.argmax(np.abs(g_vec)))
    g_nz = strict_int(g_vec[idx_nz])

    H_J_c0 = np.zeros((12, 12), dtype=object)
    H_J_c1 = np.zeros((12, 12), dtype=object)
    for i in range(12):
        b_i = B[:, i]
        M_bi_5lam = sum(strict_int(b_i[k]) * C_pow[k] for k in range(24)) @ M_5lam
        for j in range(12):
            raw_tr = Tr_KF_mat @ (M_bi_5lam @ B_bar[:, j])
            tr_vec = exact_div_vec(raw_tr, 5, f"tr_vec({i},{j})")
            tr_fq = strict_int(tr_vec[0] + (conj_mat @ tr_vec)[0])
            diff_val = strict_int((tr_vec - conj_mat @ tr_vec)[idx_nz])
            c1 = exact_div(diff_val, g_nz, f"c1({i},{j})")
            c0 = exact_div(tr_fq - c1, 2, f"c0({i},{j})")
            H_J_c0[i, j] = c0
            H_J_c1[i, j] = c1

    H_J_div_theta_c0 = exact_div_mat(H_J_c0 + 4 * H_J_c1, 7, "H_J_div_theta_c0")
    H_J_div_theta_c1 = exact_div_mat(-(2 * H_J_c0 + H_J_c1), 7, "H_J_div_theta_c1")

    print("[*] Stage 2: Constructing exact E_8 factor and base lattice G_raw...", flush=True)
    A = np.array([[2, 3], [3, 5]], dtype=object)
    D = np.eye(2, dtype=object)
    R = np.array([[3, 2], [5, 3]], dtype=object)

    if not np.array_equal(7 * A - R @ R.T, D):
        raise ArithmeticIntegrityError("Schur inversion identity 7A - RR^T = I_2 failed!")

    H_48_c0 = np.zeros((48, 48), dtype=object)
    H_48_c1 = np.zeros((48, 48), dtype=object)
    H_48_c0[:24, :24] = np.kron(A, H_J_c0)
    H_48_c1[:24, :24] = np.kron(A, H_J_c1)
    H_48_c0[24:, 24:] = np.kron(D, H_J_c0)
    H_48_c1[24:, 24:] = np.kron(D, H_J_c1)

    UR_c0 = np.kron(R, H_J_div_theta_c0)
    UR_c1 = np.kron(R, H_J_div_theta_c1)
    LL_c0 = UR_c0.T + UR_c1.T
    LL_c1 = -UR_c1.T

    H_48_c0[:24, 24:] = UR_c0
    H_48_c1[:24, 24:] = UR_c1
    H_48_c0[24:, :24] = LL_c0
    H_48_c1[24:, :24] = LL_c1

    G_96_raw = np.zeros((96, 96), dtype=object)
    G_96_raw[0::2, 0::2] = 2 * H_48_c0 + H_48_c1
    G_96_raw[1::2, 0::2] = H_48_c0 - 3 * H_48_c1
    G_96_raw[0::2, 1::2] = H_48_c0 + 4 * H_48_c1
    G_96_raw[1::2, 1::2] = 4 * H_48_c0 + 2 * H_48_c1
    G_raw_list = [[strict_int(G_96_raw[r, c]) for c in range(96)] for r in range(96)]

    validate_integer_matrix(G_raw_list, 96, 96, "G_raw")
    if not all(G_raw_list[i][j] == G_raw_list[j][i] for i in range(96) for j in range(96)):
        raise ArithmeticIntegrityError("G_raw is asymmetric!")
    print("[+] Certified G_raw: Symmetric matrix in M_96(Z)", flush=True)

    print("[*] Stage 3: Computing 48D discriminant kernel over F_7...", flush=True)
    K_basis = nullspace_mod_p(G_raw_list, p=7)
    if len(K_basis) != 48:
        raise ArithmeticIntegrityError(f"Discriminant kernel dimension != 48: {len(K_basis)}")

    Z = []
    for i in range(48):
        z_row = []
        for j in range(96):
            s_val = sum(K_basis[i][k] * G_raw_list[k][j] for k in range(96))
            z_row.append(exact_div(s_val, 7, f"Z[{i},{j}]"))
        Z.append(z_row)

    B_K = mat_mul_mod_p(Z, transpose(K_basis), p=7)

    print("[*] Stage 4: Constructing deterministic symplectic isotropic subspace C_iso...", flush=True)
    C_iso = construct_deterministic_isotropic_subspace(B_K, p=7)
    iso_check = mat_mul_mod_p(mat_mul_mod_p(C_iso, B_K, p=7), transpose(C_iso), p=7)
    if not all(iso_check[i][j] == 0 for i in range(24) for j in range(24)):
        raise ArithmeticIntegrityError("C_iso is not Lagrangian isotropic mod 7!")

    V_glued = mat_mul_mod_p(C_iso, K_basis, p=7)
    R_glue, pivots = rref_mod_p(V_glued, p=7)
    if len(pivots) != 24:
        raise ArithmeticIntegrityError(f"Glued subspace rank != 24: {len(pivots)}")
    free_cols = [c for c in range(96) if c not in pivots]

    B_int = [[0] * 96 for _ in range(96)]
    for r in range(24):
        for c in range(96):
            B_int[r][c] = strict_int(R_glue[r][c])
    for j, f in enumerate(free_cols):
        B_int[24 + j][f] = 7

    validate_integer_matrix(B_int, 96, 96, "B_int")

    print("[*] Auditing index of gluing extension [Lambda_96 : Lambda_raw]...", flush=True)
    B_int_sp = sp.Matrix(B_int)
    det_B = abs(strict_int(B_int_sp.det()))
    expected_det_B = 7**72
    if det_B != expected_det_B:
        raise ArithmeticIntegrityError(
            f"Gluing index mismatch: det(B_int) = {det_B} != 7^72 ({expected_det_B})"
        )
    print(f"[+] Certified Gluing Index: det(B_int) = 7^72", flush=True)

    print("[*] Verifying exact gluing integrality modulo 49...", flush=True)
    G_scaled = [[0] * 96 for _ in range(96)]
    for i in range(96):
        for j in range(96):
            s = sum(
                B_int[i][k] * sum(G_raw_list[k][l] * B_int[j][l] for l in range(96))
                for k in range(96)
            )
            G_scaled[i][j] = exact_div(s, 49, f"G_scaled[{i},{j}]")

    G_96_list = G_scaled
    validate_integer_matrix(G_96_list, 96, 96, "G_96")
    if not all(G_96_list[i][i] % 2 == 0 for i in range(96)):
        raise ArithmeticIntegrityError("G_96 diagonal contains odd values!")
    if not all(G_96_list[i][j] == G_96_list[j][i] for i in range(96) for j in range(96)):
        raise ArithmeticIntegrityError("G_96 is asymmetric!")

    build.B_int = B_int; build.G_raw = G_raw_list
    return G_96_list


if __name__ == "__main__":
    G = build()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "external", "srfp_G96.txt")
    with open(out, "w") as f:
        f.write("[" + "; ".join(", ".join(str(v) for v in r) for r in G) + "]\n")
    print("written", out)


def omega_stability():
    """Is the glued lattice stable under omega = (1+sqrt-7)/2 (acting on the raw Z[omega]-coordinates)?"""
    import cypari2
    pari = cypari2.Pari()
    src = open(__file__).read()
    ns = {}
    G = build()                      # rebuild; we need B_int, so recompute it here
    return G
