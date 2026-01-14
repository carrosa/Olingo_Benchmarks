import sys
import time

sys.path.append("..")  # path to LaZer module
from lazer import *  # import LaZer Python module
import hashlib  # for SHAKE128
from lazer import lin_prover_state_t, lin_verifier_state_t  # Explicit import
from _ctxr_params_cffi import lib


def decompose_vector(m: polyvec_t, B: int, ring: polyring_t):
    m0 = []
    for pol in m.to_pol_list():
        tmp = []
        for coeff in pol.to_list():
            tmp.append(coeff % B)
        m0.append(poly_t(ring, tmp))
    m1 = []
    for pol in m.to_pol_list():
        tmp = []
        for coeff in pol.to_list():
            tmp.append(coeff // B)
        m1.append(poly_t(ring, tmp))
    return m0, m1


def check_large_coeffs(label, vec, qhat):
    print(f"Checking {label}...")
    for i, pol in enumerate(vec.to_pol_list()):
        coeffs = pol.to_list()
        for j, c in enumerate(coeffs):
            if not (0 <= c < qhat * 2):
                print(f"[{label}][{i}][{j}] = {c} is too large")
                return (i, j, c)
    print(f"{label}: all coefficients < 2·qhat")
    return None


def main():
    # === Setup ===
    from ctxr_params import deg, q_hat, dim, cmt_n, cmt_k, cmt_l, secdim_hat, m

    seed = b"\0" * 32  # public randomness (proof seed)
    prover = lin_prover_state_t(seed, lib.get_params("param"))
    verifier = lin_verifier_state_t(seed, lib.get_params("param"))

    seed_1 = b"\01"
    seed_2 = b"\02"
    seed_3 = b"\03"
    seed_4 = b"\04"

    Rqhat = polyring_t(deg, q_hat)

    one_pol = poly_t(
        Rqhat, [1] + [0] * (deg - 1)
    )  # Polynomial with one in the first coefficient

    A1_prime = polymat_t.urandom_static(Rqhat, cmt_n, cmt_k - cmt_n, q_hat, seed_1, 0)
    A1 = polymat_t(Rqhat, cmt_n, cmt_k, [one_pol, A1_prime])
    A_top_right = polymat_t(Rqhat, 1, cmt_l)
    A_top = polymat_t(Rqhat, 1, cmt_k + cmt_l, [A1, A_top_right])

    A2_prime = polymat_t.urandom_static(
        Rqhat, cmt_l, cmt_k - cmt_n - cmt_l, q_hat, seed_2, 0
    )
    A2_mid = polymat_t.identity(Rqhat, cmt_l)
    A2_left = polymat_t(Rqhat, cmt_l, cmt_n)
    A2 = polymat_t(Rqhat, cmt_l, cmt_k, [A2_left, A2_mid, A2_prime])
    I22 = polymat_t.identity(Rqhat, cmt_l)
    A_bot = polymat_t(Rqhat, cmt_l, cmt_k + cmt_l, [A2, I22])

    A = polymat_t(Rqhat, cmt_n + cmt_l, cmt_k + cmt_l)
    for row in range(1, cmt_n + cmt_l):
        A.set_row(row, A_bot.get_row(row - 1))
    A.set_row(0, A_top.get_row(0))

    # Binary
    r = polyvec_t.urandom_bnd_static(Rqhat, cmt_k, 0, 1, seed_3, 0)
    S = polyvec_t.urandom_bnd_static(Rqhat, cmt_l, 0, 1, seed_4, 0)
    # E = polyvec_t.grandom_static(Rqhat, 22, int(log2o), seed_4, 0)
    w = polyvec_t(Rqhat, cmt_k + cmt_l, [r, S])

    t = A * w

    # === Prover setup ===
    print("set_statement...")
    prover.set_statement(A, -t)

    print("set_witness...")
    print(type(w))

    prover.set_witness(w=w)

    print("generate proof ...")
    start_time = time.time()
    proof = prover.prove()
    print("proof done")
    end_time = time.time()

    # Report proof time and size
    elapsed_ms = (end_time - start_time) * 1000
    elapsed_ms_xm = (end_time - start_time) * 1000 * m
    proof_size_bits = len(proof) * 8
    proof_size_bits_xm = len(proof) * 8 * m
    proof_size_kb = proof_size_bits / 8000  # 1 KB = 8000 bits
    proof_size_kb_xm = proof_size_bits_xm / 8000  # 1 KB = 8000 bits
    print(f"Prover time: {elapsed_ms:.2f} ms per column")
    print(f"Prover time time: {elapsed_ms_xm:.2f} ms for {secdim_hat}x{m}")
    print(f"Proof size: {proof_size_bits} bits ({proof_size_kb:.2f} KB) per column")
    print(f"Proof size: {proof_size_bits_xm} bits ({proof_size_kb_xm:.2f} KB) for {secdim_hat}x{m}")

    # === Verifier ===
    verifier.set_statement(A, -t)

    print("verify proof ... ")
    try:
        start_time = time.time()
        verifier.verify(proof)
        end_time = time.time()
        elapsed_ms = (end_time - start_time) * 1000
        elapsed_ms_xm = (end_time - start_time) * 1000 * m
        print(f"Verify time: {elapsed_ms:.2f} ms per column")
        print(f"Verify time: {elapsed_ms_xm:.2f} ms for {secdim_hat}x{m}")
    except VerificationError as e:
        print("reject")
        print(e)
    else:
        print("accept")


if __name__ == "__main__":
    main()
