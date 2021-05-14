//------------------------------------------------------------------------------
// SPEX_Update/SPEX_Update_Chol_Rank1: perform Cholesky rank-1 update
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is used to perform Cholesky rank-1 update for
// A' = A + sigma*w*w^T, where A' and A are n-by-n matrix, w is n-by-1 vector
// and sigma is a scalor (positive for update and negative for downdate). L
// is the given Cholesky factorization for A such that A = L*D^(-1)*L^T. sd is
// the diagonal of L, P is the given row permutation for the Cholesky
// factorization such that L(P,:) is a lower triangular matrix, and P_inv is
// the inverse of P. The result of of this function will have L, w updated as
// L_out*D^(-1)*L_out^T = A' and L*D^(-1)*w_out = w.

#define SPEX_FREE_ALL               \
    SPEX_FREE(h);                   \
    SPEX_MPQ_CLEAR(sd_ratio);       \
    SPEX_MPQ_CLEAR(pending_scale);  \
    SPEX_MPQ_CLEAR(tmpq);           \
    SPEX_MPZ_CLEAR(sd0_old);        \
    SPEX_MPZ_CLEAR(sd1_old);        \
    SPEX_MPZ_CLEAR(tmpz);

#include "spex_update_internal.h"

SPEX_info SPEX_Update_Chol_Rank1
(
    SPEX_matrix *L,
    SPEX_matrix *rhos,// array of scaled pivots
    int64_t *P,
    int64_t *P_inv,
    SPEX_vector *w,
    int64_t sigma
)
{
    SPEX_REQUIRE(L, SPEX_DYNAMIC_CSC, SPEX_MPZ);
    SPEX_REQUIRE(rhos, SPEX_DENSE, SPEX_MPZ);
    // initialize workspace
    SPEX_info info;
    int sgn;
    int64_t i, j, p, n = L->n, w_top = 0;
    int64_t *h = NULL;
    mpz_t *sd = rhos->x.mpz;
    mpq_t sd_ratio;// = sd_new/sd_old
    mpq_t tmpq, pending_scale;
    mpz_t tmpz;
    // sd0_old and sd1_old are the original values of sd[j-1]
    // and sd[j], respectively.
    mpz_t sd0_old, sd1_old;
    SPEX_MPQ_SET_NULL(sd_ratio);
    SPEX_MPQ_SET_NULL(pending_scale);
    SPEX_MPQ_SET_NULL(tmpq);
    SPEX_MPZ_SET_NULL(tmpz);
    SPEX_MPZ_SET_NULL(sd0_old);
    SPEX_MPZ_SET_NULL(sd1_old);
    SPEX_CHECK(SPEX_mpq_init(sd_ratio));
    SPEX_CHECK(SPEX_mpq_init(pending_scale));
    SPEX_CHECK(SPEX_mpq_init(tmpq));
    SPEX_CHECK(SPEX_mpz_init(tmpz));
    SPEX_CHECK(SPEX_mpz_init(sd0_old));
    SPEX_CHECK(SPEX_mpz_init(sd1_old));

    // allocate and initialize h as all -1s
    h = (int64_t*) SPEX_malloc(n * sizeof(int64_t));
    if (!h)
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }
    for (i = 0; i < n; i++)
    {
        h[i] = -1;
    }

    // initialize for the loop
    // sd1_old = 1
    SPEX_CHECK(SPEX_mpz_set_ui(sd1_old, 1));
    // sd_ratio = 1
    SPEX_CHECK(SPEX_mpq_set_ui(sd_ratio, 1, 1));
    // update column j of L
    for (j = 0; j < n; j++)
    {
        int64_t Pj = P[j], hj = h[Pj], hi;

        // sd0_old = sd1_old, sd1_old = sd[j]
        SPEX_CHECK(SPEX_mpz_swap(sd0_old, sd1_old));
        SPEX_CHECK(SPEX_mpz_set(sd1_old, sd[j]));

        // pending_scale = sd_ratio *S(j)
        SPEX_CHECK(SPEX_mpq_mul(pending_scale, sd_ratio, L->v[j]->scale));

        SPEX_CHECK(SPEX_mpz_sgn(&sgn, w->x[Pj]));
        if (sgn != 0)// w[P[j]] != 0
        {
            // check if history update is needed for w[P[j]]
            if (hj < j-1)
            {
                // perform history update
                // w[P[j]] = w[P[j]] * sd_new[j-1]/sd_new[h[P[j]]]
                SPEX_CHECK(SPEX_mpz_mul(w->x[Pj], w->x[Pj], sd[j-1]));
                if (hj > -1)
                {
                    SPEX_CHECK(SPEX_mpz_divexact(w->x[Pj], w->x[Pj], sd[hj]));
                }
            }

            // tmpq = sigma*w[P[j]]*w[P[j]]/(sd_old[j]*sd_old[j-1])
            SPEX_CHECK(SPEX_mpz_mul   (SPEX_MPQ_NUM(tmpq), w->x[Pj], w->x[Pj]));
            SPEX_CHECK(SPEX_mpz_mul_si(SPEX_MPQ_NUM(tmpq),
                                       SPEX_MPQ_NUM(tmpq), sigma));
            SPEX_CHECK(SPEX_mpz_mul   (SPEX_MPQ_DEN(tmpq), sd0_old, sd1_old));
            SPEX_CHECK(SPEX_mpq_canonicalize(tmpq));

            // sd_ratio += tmpq, which gives sd_new[j]/sd_old[j]
            SPEX_CHECK(SPEX_mpq_add(sd_ratio, sd_ratio, tmpq));

            // update sd_new[j] as sd_old[j]*sd_ratio
            SPEX_CHECK(SPEX_mpz_divexact(sd[j], sd[j], SPEX_MPQ_DEN(sd_ratio)));
            SPEX_CHECK(SPEX_mpz_mul     (sd[j], sd[j], SPEX_MPQ_NUM(sd_ratio)));

            // set L(j,j) = sd_new[j]
            ASSERT(L->v[j]->i[0] == Pj);
            SPEX_CHECK(SPEX_mpz_set(L->v[j]->x[0], sd[j]));

            // reset S(j) = 1
            SPEX_CHECK(SPEX_mpq_set_ui(L->v[j]->scale, 1, 1));

            // iterate across nnz in L->v[j] but skip the pivot entry L(j,j)
            for (p = 1; p < L->v[j]->nz; p++)
            {
                i = L->v[j]->i[p];
                hi = h[i];
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, w->x[i]));
                if (sgn == 0) // w[i] == 0
                {
                    SPEX_CHECK(SPEX_mpz_sgn(&sgn, L->v[j]->x[p]));
                    if (sgn == 0) {continue;}

                    // update L(i,j)
                    SPEX_CHECK(SPEX_mpz_mul     (L->v[j]->x[p],
                                                 L->v[j]->x[p],
                                                 SPEX_MPQ_NUM(pending_scale)));
                    SPEX_CHECK(SPEX_mpz_divexact(L->v[j]->x[p],
                                                 L->v[j]->x[p],
                                                 SPEX_MPQ_DEN(pending_scale)));

                    // perform IPGE to update w using updated L(i,j). Since w[i]
                    // is zero, this could be a fill-in.
                    // w[i] = -L(i,j)*w[P[j]]/sd_new[j-1]
                    SPEX_CHECK(SPEX_mpz_submul(w->x[i],
                                               w->x[Pj], L->v[j]->x[p]));
                    if (j != 0)
                    {
                        SPEX_CHECK(SPEX_mpz_divexact(w->x[i], w->x[i],sd[j-1]));
                    }
                    // add this entry to nnz pattern of w
                    w->i[w->nz] = i;
                    w->nz++;
                }
                else // w[i] != 0
                {
                    // check if history update is needed
                    if (hi < j-1)
                    {
                        // perform history update
                        // w[i] = w[i] * sd_new[j-1]/sd_new[h[i]]
                        SPEX_CHECK(SPEX_mpz_mul(w->x[i], w->x[i], sd[j-1]));
                        if (hi > -1)
                        {
                            SPEX_CHECK(SPEX_mpz_divexact(w->x[i],
                                                         w->x[i], sd[hi]));
                        }
                    }

                    // tmpz = sigma*w[i]*w[P[j]]
                    SPEX_CHECK(SPEX_mpz_mul(tmpz, w->x[i], w->x[Pj]));
                    SPEX_CHECK(SPEX_mpz_mul_si(tmpz, tmpz, sigma));
#ifdef SPEX_DEBUG
                    // tmpz /= sd_old[j-1]
                    mpq_t r1, r2; mpq_init(r1); mpq_init(r2);
                    mpz_fdiv_qr(tmpz, SPEX_MPQ_NUM(r1),
                                tmpz, sd0_old);
                    mpq_set_den(r1, sd0_old);
                    mpq_canonicalize(r1);

                    SPEX_CHECK(SPEX_mpz_mul(L->v[j]->x[i],
                                            L->v[j]->x[i],
                                            SPEX_MPQ_NUM(pending_scale)));
                    mpz_cdiv_qr(L->v[j]->x[i], SPEX_MPQ_NUM(r2),
                                L->v[j]->x[i], SPEX_MPQ_DEN(pending_scale));
                    mpq_set_den(r2, SPEX_MPQ_DEN(pending_scale));
                    mpq_canonicalize(r2);
                    mpq_neg(r2, r2);
                    if (mpq_cmp(r1, r2) != 0)
                    {
                        mpq_clear(r1);
                        mpq_clear(r2);
                        SPEX_CHECK(SPEX_PANIC);
                    }
                    mpq_clear(r1);
                    mpq_clear(r2);
#else
                    // tmpz /= sd_old[j-1]
                    SPEX_CHECK(SPEX_mpz_fdiv_q(tmpz, tmpz, sd0_old));

                    // L(i,j) *= sd_ratio
                    SPEX_CHECK(SPEX_mpz_mul(L->v[j]->x[i],
                                            L->v[j]->x[i],
                                            SPEX_MPQ_NUM(pending_scale)));
                    SPEX_CHECK(SPEX_mpz_cdiv_q(L->v[j]->x[i],
                                            L->v[j]->x[i],
                                            SPEX_MPQ_DEN(pending_scale)));
#endif
                    SPEX_CHECK(SPEX_mpz_add(L->v[j]->x[i],
                                            L->v[j]->x[i], tmpz));

                    // ---------------------------------------------------------
                    // perform IPGE to update w using updated L(i,j).
                    // ---------------------------------------------------------
                    SPEX_CHECK(SPEX_mpz_sgn(&sgn, L->v[j]->x[p]));
                    if (sgn == 0)
                    {
                        // h[i] is supposed to be j-1 at this moment, but it is
                        // instead set as -(j-1)-3 = -j-2 to mark that L(i,j)
                        // is already in the nnz pattern of L(:,j)
                        h[i] = -j-2;
                        continue;
                    }

                    // w[i] = w[i]*sd_new[j]
                    SPEX_CHECK(SPEX_mpz_mul(w->x[i], w->x[i], sd[j]));
                    // w[i] = (w[i]-L(i,j)*w[P[j]])/sd_new[j-1]
                    SPEX_CHECK(SPEX_mpz_submul(w->x[i],
                                            w->x[Pj], L->v[j]->x[p]));
                    if (j != 0)
                    {
                        SPEX_CHECK(SPEX_mpz_divexact(w->x[i], w->x[i],sd[j-1]));
                    }
                }
                h[i] = j;
            }

            // iterate across all nnz in w to check if fillin should be added
            // to L->v[j]
            int64_t Lj_nz = L->v[j]->nz;
            for (p = w_top; p < w->nz; p++)
            {
                i = w->i[p];
                hi = h[i];
                if (P_inv[i] <= j)
                {
                    w->i[p] = w->i[w_top];
                    w->i[w_top] = i;
                    w_top++;
                    continue;
                }
                if (hi < -1) // updated L(i,j) becomes 0 
                {
                    // restore value of h[i]
                    h[i] = -hi-3;
                    continue;
                }
                else if (h[i] != j) // this is a fillin to L->v[j]
                {
                    SPEX_CHECK(SPEX_mpz_sgn(&sgn, w->x[i]));
                    if (sgn == 0) {continue;}

                    // perform history update to w[i]
                    // w[i] = w[i]*sd_new[j-1]/sd_new[h[i]]
                    SPEX_CHECK(SPEX_mpz_mul(w->x[i], w->x[i], sd[j-1]));
                    if (hi > -1)
                    {
                        SPEX_CHECK(SPEX_mpz_divexact(w->x[i], w->x[i], sd[hi]));
                    }

                    // allocate extra space if needed
                    if (Lj_nz == L->v[j]->nzmax)
                    {
                        SPEX_CHECK(SPEX_vector_realloc(L->v[j],
                            SPEX_MIN(n, 2*(L->v[j]->nzmax)), NULL));
                    }

                    // L(i,j) = sigma*w[P[j]]*w[i]/sd_old[j-1]
                    SPEX_CHECK(SPEX_mpz_mul   (L->v[j]->x[Lj_nz],
                                               w->x[Pj], w->x[i]));
                    SPEX_CHECK(SPEX_mpz_mul_si(L->v[j]->x[Lj_nz],
                                               L->v[j]->x[Lj_nz], sigma));
                    SPEX_CHECK(SPEX_mpz_divexact(L->v[j]->x[Lj_nz],
                                               L->v[j]->x[Lj_nz], sd0_old));
                    // insert new entry to L
                    L->v[j]->i[Lj_nz] = i;
                    Lj_nz++;

                    // w[i] can be efficiently updated as
                    // w[i] = w[i]*sd_old[j]/sd_old[j-1]
                    SPEX_CHECK(SPEX_mpz_mul(w->x[i], w->x[i], sd1_old));
                    SPEX_CHECK(SPEX_mpz_divexact(w->x[i], w->x[i], sd0_old));
                    h[i] = j;
                }
            }
            L->v[j]->nz = Lj_nz;
        }
        else // w[P[j]] == 0
        {
            SPEX_CHECK(SPEX_mpq_set(L->v[j]->scale, pending_scale));
            SPEX_CHECK(SPEX_mpz_divexact(sd[j], sd[j], SPEX_MPQ_DEN(sd_ratio)));
            SPEX_CHECK(SPEX_mpz_mul     (sd[j], sd[j], SPEX_MPQ_NUM(sd_ratio)));
        }
    }

    SPEX_FREE_ALL;
    return SPEX_OK;
}
