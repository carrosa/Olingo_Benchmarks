#ifndef GHKSS256_H
#define GHKSS256_H

#include "params.h"
#include "poly256.h"
#include "randombytes.h"
#include "reduce256.h"
#include "symmetric.h"
#include <gmp.h>

typedef struct {
  poly *u; // size LHAT
  poly *v; // size M
} ctx_t;

typedef struct {
  poly c;
  poly *z; // size L
  poly *h; // size K
} sig_t;

typedef struct {
  poly (*A)[L]; // [K][L]
  poly *yprime; // [K]
} pks_t;

typedef struct {
  poly (*A)[LHAT]; // [KHAT][LHAT]
  poly (*B)[M];    // [KHAT][M]
  poly (*S)[M];    // [LHAT][M]
} pke_t;

typedef struct {
  poly *share; // [L]
} decryption_share_t;

void gen_keys(poly A[M][M], poly B[M][L], poly S[M][L]);

void encrypt(poly U[M][K], poly V[L][K], poly A[M][M], poly B[M][L],
             poly Msg[L][K]);

void decrypt(poly Msg[L][K], poly S[M][L], poly U[M][K], poly V[L][K]);

void gen_sharing_poly(poly (*P_coeff)[LHAT][M], poly (*S)[M]);

// void compute_share(poly Share[LHAT][M], poly P_coeff[][LHAT][M], int32_t
// user);
void compute_share(poly (*Share)[M], poly (*P_coeff)[LHAT][M], int32_t user);

void compute_share_horner(poly (*Share)[M],         // [LHAT][M]
                          poly (*P_coeff)[LHAT][M], // [THRESHOLD][LHAT][M]
                          int32_t user);            // users[USERS]

void reconstruct_secret(poly (*S)[M], int32_t *users, poly (*Shares)[LHAT][M]);
// void reconstruct_secret(poly S[LHAT][M], int32_t *users,
//                         poly Shares[][LHAT][M]);

void tdecrypt(poly ds_i[L][K], poly S_i[M][L], poly U[M][K], int32_t user,
              int32_t *users, int t);

void combine(poly Msg[L][K], poly V[L][K], poly (*ds)[L][K], int usize);

void gen_randomness(uint8_t seed[SEEDBYTES]);
void lagrange_coeff(mpz_t lambda, int32_t user, int32_t *users, int t);

// Signing
void sigkeygen(poly A[K][L], poly yprime[L], poly s[L]);

void sign1p(poly s[L], poly A[K][L], poly yprime[K], uint8_t *mu, size_t mu_len,
            poly c, poly z[L], poly h[K]);

int verify(poly A[K][L], poly yprime[K], poly c, poly z[L], poly h[K],
           uint8_t *mu, size_t mu_len);

// void compute_challenge(poly *c, poly A[K][L], poly yprime[K], poly wprime[K],
//                        poly mu[L]);
void compute_challenge(poly *c,
                       poly (*A)[L],   // [K][L]
                       poly yprime[K], // [K]
                       poly wprime[K], // [K]
                       poly mu[L]);    // [L]

void keygen_e_1d(poly (*A)[LHAT], // size: [KHAT][LHAT]
                 poly (*B)[M],    // size: [KHAT][M]
                 poly (*S)[M]);   // size: [LHAT][M]

void encrypt_1d(poly *U,         // [LHAT]
                poly *V,         // [M]
                poly (*A)[LHAT], // [KHAT][LHAT]
                poly (*B)[M],    // [KHAT][M]
                poly *Msg);      // [M]

void decrypt_1d(poly Msg[M], poly S[LHAT][M], poly U[LHAT], poly V[M]);

// void tdecrypt_1d(poly ds_i[M], poly S_i[LHAT][M], poly U[LHAT], int32_t user,
//                  int32_t *users, int t);
void tdecrypt_1d(poly *ds_i, poly (*S_i)[M], poly *U, int32_t user,
                 int32_t *users, int t);

void combine_1d(poly *Msg,     // Output: poly[M]  — accumulator of the result
                poly *V,       // Input: poly[M]   — original ciphertext vector
                poly (*ds)[M], // Input: [usize][M] — decryption shares
                int usize);    // usize = number of users/shares to combine

void as_keygen_preenc(poly (*As)[L],     // [K][L], NTT domain
                      poly *si,          // [L], initialized
                      poly *yi,          // [K], initialized
                      uint8_t h_yi[32]); // output hash
void as_keygen_encrypt_only(
    ctx_t *ctx,       // ctx->u: [LHAT], ctx->v: [M], initialized
    poly (*Ae)[LHAT], // [KHAT][LHAT]
    poly (*Be)[M],    // [KHAT][M]
    poly *si);        // [L], NTT domain (as produced by as_keygen_preenc)
void as_sign_r1(poly (*As)[L],   // [K][L], NTT domain
                uint8_t h_w[32], // output hash of w
                poly **r_out);   // output: allocated r[M] (NTT domain)
void as_sign_r2(
    poly (*Ae)[LHAT], // [KHAT][LHAT]
    poly (*Be)[M],    // [KHAT][M]
    poly *u,          // [LHAT], initialized polys
    poly *v,          // [M],    initialized polys
    poly *r,          // [M], NTT domain (from as_sign_r1)
    int want_free_r); // if non-zero, this clears+frees r after encrypt
void as_sign_round2(poly (*Ae)[LHAT], // [KHAT][LHAT]
                    poly (*Be)[M],    // [KHAT][M]
                    poly *u,          // [LHAT]
                    poly *v,          // [M]
                    poly *r           // M
);
// void as_keygen_1(poly As[K][L], poly Ae[KHAT][LHAT], poly Be[KHAT][M],
//                  poly *si, // L
//                  poly *yi, // K
//                  uint8_t h_yi[32], ctx_t *ctx);

void as_keygen_1(poly (*As)[L],    // [K][L]
                 poly (*Ae)[LHAT], // [KHAT][LHAT]
                 poly (*Be)[M],    // [KHAT][M]
                 poly *si,         // [L]
                 poly *yi,         // [K]
                 uint8_t h_yi[32],
                 ctx_t *ctx); // ctx->u: [LHAT], ctx->v: [M]

void as_keygen_2(ctx_t *ctx_s,        // ctx_s->u [LHAT], ctx_s->v [M]
                 poly (*yi)[K],       // yi[USERS][K]
                 uint8_t (*h_yi)[32], // h_yi[USERS][32]
                 ctx_t *ctx,          // ctx[USERS]
                 poly *yprime);       // yprime[K]

void as_sign_1(poly (*As)[L],    // [K][L]
            //    poly (*Ae)[LHAT], // [KHAT][LHAT]
            //    poly (*Be)[M],    // [KHAT][M]
            //    poly *u,          // [LHAT]
            //    poly *v,          // [M]
               poly *r           // M
);

void as_sign_1_old(poly (*As)[L],    // [K][L]
                   poly (*Ae)[LHAT], // [KHAT][LHAT]
                   poly (*Be)[M],    // [KHAT][M]
                   ctx_t ctx_si);    // ctx_si.u [LHAT], ctx_si.v [M]

void as_sign_2(poly *c,           // just one poly
               poly (*A_pks)[L],  // KxL
               poly *yprime,      // K
               poly (*u_r)[LHAT], // THRESHOLD x LHAT
               poly (*v_r)[M],    // THRESHOLD x M
               poly *u_s,         // LHAT
               poly *v_s,         // M
               poly *u_z,         // LHAT
               poly *v_z,         // M
               poly (*ski)[M],    // LHAT x M
               poly *mu,          // L, but padded to M
               poly *dsi,         // L, but padded to M
               poly (*w_i)[K],    // THRESHOLD x K
               int user, int32_t *users);

void as_sign_round3(sig_t *sig, pks_t pks,
                    ctx_t ctx_r[THRESHOLD], // Should be THRESHOLD?
                    ctx_t ctx_s, ctx_t ctx_z, poly ski[LHAT][M], poly mu[M],
                    poly dsi[M], poly w_i[THRESHOLD][K],
                    int user, int32_t users[USERS]);

// void as_sign_3(
//     ctx_t *ctx_z,
//     const pks_t *pks,
//     sig_t *sig,
//     const poly *yprime, // K
//     poly *wprime,       // K
//     poly dsi[THRESHOLD][M],
//     const int32_t *users); // USERS

void as_sign_comb(ctx_t *ctx_z,       // ctx_z->v [M]
                  const pks_t *pks,   // pks->A [K][L]
                  sig_t *sig,         // sig->c, sig->z [L], sig->h [K]
                  const poly *yprime, // [K]
                  poly *wprime,       // [K]
                  poly (*dsi)[M]      // [THRESHOLD][M]
);

void as_sign_3(poly *u_z,        // LHAT
               poly *v_z,        // M
               poly (*A_pks)[L], // KxL
               poly *c_sig,      // just one poly
               poly *z_sig,      // L
               poly *h_sig,      // K
               poly *yprime,     // K
               poly *wprime,     // K
               poly *dsi,        // L, but padded to M
               int32_t *users);

int as_verify_sig(sig_t *sig, pks_t *pks, poly *mu);

void dk_gen_1(poly (*Ae)[LHAT], // Ae: KHAT x LHAT
              poly (*Si)[M],    // Si: LHAT x M
              poly (*Ei)[M],    // Ei: KHAT x M
              poly (*Bi)[M],    // Bi: KHAT x M
              uint8_t *h_Bi);   // 32-byte output hash

void dk_gen_2(poly (*Si)[M],    // [LHAT][M]
              poly (*Ei)[M],    // [LHAT][M]
              poly (*Ae)[LHAT], // [KHAT][LHAT]
              int32_t *users);  // users[USERS]

// void as_keygen_2(ctx_t *ctx_s,        // ctx_s->u [LHAT], ctx_s->v [M]
//                  poly (*yi)[K],       // yi[USERS][K]
//                  uint8_t (*h_yi)[32], // h_yi[USERS][32]
//                  ctx_t *ctx,          // ctx[USERS]
//                  poly *yprime);       // yprime[K]

void dk_gen_3(poly (*Bj)[KHAT][M],  // [USERS][KHAT][M]
              poly (*Sij)[LHAT][M], // [USERS][LHAT][M]
              poly (*Be)[M],        // [KHAT][M]
              poly (*Sei)[M]);      // [LHAT][M]

#endif