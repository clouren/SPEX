//------------------------------------------------------------------------------
// SPEX_Update/SPEX_update_cholesky_rank1: perform Cholesky rank-1 update
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2023, Chris Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// Purpose: This function performs rank-1 Cholesky update/downdate. The input
// factorization needs to be updatable with L being SPEX_DYNAMIC_CSC MPZ
// matrix. Otherwise (if F is non-updatable upon input), this function calls
// SPEX_factorization_convert to make F updatable, which requires that L in the
// factorization must be non-shallow SPEX_CSC MPZ matrix. The output
// factorization will always be updatable.  Since the factorization is modified
// during the update process, the returned F should be considered as undefined
// if this function fails for any reason.

// The matrix w is modified during the update. If the updated A is needed,
// user can compute A = A + sigma*w*w' *BEFORE* using this function (since w
// will be modified).

#define SPEX_FREE_ALL               \
{                                   \
    spex_scattered_vector_free(&w_dense, option); \
    SPEX_FREE(h);                   \
    SPEX_MPQ_CLEAR(sd_ratio);       \
    SPEX_MPQ_CLEAR(pending_scale);  \
    SPEX_MPQ_CLEAR(tmpq);           \
    SPEX_MPZ_CLEAR(sd0_old);        \
    SPEX_MPZ_CLEAR(sd1_old);        \
    SPEX_MPZ_CLEAR(tmpz);           \
}

#include "spex_update_internal.h"

// TODO allow w->v[0]->scale != 1

SPEX_info SPEX_update_cholesky_rank1
(
    SPEX_factorization F,   // The SPEX Cholesky factorization of A, including
                            // L, rhos, P and Pinv. This factorization will be
                            // modified during the update process. Therefore,
                            // if this function fails for any reason, the
                            // returned F should be considered as undefined.
    SPEX_matrix w,          // a n-by-1 dynamic_CSC matrix that contains the
                            // vector to modify the original matrix A, the
                            // resulting A is A+sigma*w*w^T. A->scale = w->scale
                            // and w->v[0]->scale = 1. In output, w is
                            // updated as the solution to L*D^(-1)*w_out = w
    const int64_t sigma,    // a nonzero scalar that determines whether
                            // this is an update (sigma > 0) or downdate
                            // (sigma < 0).
    const SPEX_options option // Command options
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------
    SPEX_info info;
    if (!spex_initialized()) {return SPEX_PANIC;}

    SPEX_REQUIRE(w , SPEX_DYNAMIC_CSC, SPEX_MPZ);

    if (!F || F->kind != SPEX_CHOLESKY_FACTORIZATION || sigma == 0 ||
        w->m != F->L->m || w->n != 1 || w->v[0]->nz < 0 ||
        (w->v[0]->nz > 0 && (!(w->v[0]->x) || !(w->v[0]->i))))
    {
        return SPEX_INCORRECT_INPUT;
    }

    // if w is a zero vector, no need to perform any update
    if (w->v[0]->nz == 0) {return SPEX_OK;}

    // make sure F is updatable
    info = SPEX_factorization_convert(F, true, option);
    if (info != SPEX_OK) return info;

    //--------------------------------------------------------------------------
    // initialize workspace
    //--------------------------------------------------------------------------
    SPEX_matrix L = F->L;
    int64_t *P = F->P_perm, *P_inv = F->Pinv_perm;
    int sgn;
    int64_t i, j, p, n = L->n, w_top = 0;
    int64_t *h = NULL;
    spex_scattered_vector w_dense = NULL;
    mpz_t *sd = F->rhos->x.mpz;
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
    SPEX_MPQ_INIT(sd_ratio);
    SPEX_MPQ_INIT(pending_scale);
    SPEX_MPQ_INIT(tmpq);
    SPEX_MPZ_INIT(tmpz);
    SPEX_MPZ_INIT(sd0_old);
    SPEX_MPZ_INIT(sd1_old);

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

    //--------------------------------------------------------------------------
    // initialize for the loop
    //--------------------------------------------------------------------------
    // get the scattered form of w
    SPEX_CHECK(spex_update_get_scattered_v(&w_dense, NULL, w->v[0], n, n, NULL,
        false, option));
    // sd1_old = 1
    SPEX_MPZ_SET_UI(sd1_old, 1);
    // sd_ratio = 1
    SPEX_MPQ_SET_UI(sd_ratio, 1, 1);

    //--------------------------------------------------------------------------
    // update column j of L
    //--------------------------------------------------------------------------
    for (j = 0; j < n; j++)
    {
        int64_t Pj = P[j], hj = h[Pj], hi;

        // sd0_old = sd1_old, sd1_old = sd[j]
        SPEX_MPZ_SWAP(sd0_old, sd1_old);
        SPEX_MPZ_SET(sd1_old, sd[j]);

        // pending_scale = sd_ratio *S(j)
        SPEX_MPQ_MUL(pending_scale, sd_ratio, L->v[j]->scale);

        SPEX_MPZ_SGN(&sgn, w_dense->x[Pj]);
        if (sgn != 0)// w[P[j]] != 0
        {
            // check if history update is needed for w[P[j]]
            if (hj < j-1)
            {
                // perform history update
                // w[P[j]] = w[P[j]] * sd_new[j-1]/sd_new[h[P[j]]]
                SPEX_MPZ_MUL(w_dense->x[Pj],
                                        w_dense->x[Pj], sd[j-1]);
                if (hj > -1)
                {
                    SPEX_MPZ_DIVEXACT(w_dense->x[Pj],
                                                 w_dense->x[Pj], sd[hj]);
                }
            }

            // tmpq = sigma*w[P[j]]*w[P[j]]/(sd_old[j]*sd_old[j-1])
            SPEX_MPZ_MUL(SPEX_MPQ_NUM(tmpq), w_dense->x[Pj], w_dense->x[Pj]);
            SPEX_MPZ_MUL_SI(SPEX_MPQ_NUM(tmpq),
                                       SPEX_MPQ_NUM(tmpq), sigma);
            SPEX_MPZ_MUL(SPEX_MPQ_DEN(tmpq), sd0_old, sd1_old);
            SPEX_MPQ_CANONICALIZE(tmpq);

            // sd_ratio += tmpq, which gives sd_new[j]/sd_old[j]
            SPEX_MPQ_ADD(sd_ratio, sd_ratio, tmpq);
            SPEX_MPQ_SGN(&sgn, sd_ratio);
            if (sgn == 0)
            {
                SPEX_FREE_ALL;
                return SPEX_SINGULAR;
            }

            // update sd_new[j] as sd_old[j]*sd_ratio
            SPEX_MPZ_DIVEXACT(sd[j], sd[j], SPEX_MPQ_DEN(sd_ratio));
            SPEX_MPZ_MUL     (sd[j], sd[j], SPEX_MPQ_NUM(sd_ratio));

            // set L(j,j) = sd_new[j]
            ASSERT(L->v[j]->i[0] == Pj);
            SPEX_MPZ_SET(L->v[j]->x[0], sd[j]);

            // reset S(j) = 1
            SPEX_MPQ_SET_UI(L->v[j]->scale, 1, 1);

            // TODO combine the above into following iteration as p=0
            // iterate across nnz in L->v[j] but skip the pivot entry L(j,j)
            for (p = 1; p < L->v[j]->nz; p++)
            {
                i = L->v[j]->i[p];
                hi = h[i];
                SPEX_MPZ_SGN(&sgn, w_dense->x[i]);
                if (sgn == 0) // w[i] == 0
                {
                    SPEX_MPZ_SGN(&sgn, L->v[j]->x[p]);
                    if (sgn == 0) {continue;}

                    // update L(i,j)
                    SPEX_MPZ_MUL     (L->v[j]->x[p],
                                                 L->v[j]->x[p],
                                                 SPEX_MPQ_NUM(pending_scale));
                    SPEX_MPZ_DIVEXACT(L->v[j]->x[p],
                                                 L->v[j]->x[p],
                                                 SPEX_MPQ_DEN(pending_scale));

                    // perform IPGE to update w using updated L(i,j). Since w[i]
                    // is zero, this could be a fill-in.
                    // w[i] = -L(i,j)*w[P[j]]/sd_new[j-1]
                    SPEX_MPZ_SUBMUL(w_dense->x[i],
                                               w_dense->x[Pj], L->v[j]->x[p]);
                    if (j != 0)
                    {
                        SPEX_MPZ_DIVEXACT(w_dense->x[i],
                                               w_dense->x[i], sd[j-1]);
                    }
                    // add this entry to nnz pattern of w if this was not in
                    // the nnz pattern. w_dense initially has no explicit zero.
                    // Therefore, any explicit zero found in the nnz
                    // pattern of w_dense must have hi > -1 (i.e., at least
                    // updated once).
                    if (hi == -1)
                    {
                        w_dense->i[w_dense->nz] = i;
                        w_dense->nz++;
                    }
                }
                else // w[i] != 0
                {
                    // check if history update is needed
                    if (hi < j-1)
                    {
                        // perform history update
                        // w[i] = w[i] * sd_new[j-1]/sd_new[h[i]]
                        SPEX_MPZ_MUL(w_dense->x[i],
                                                w_dense->x[i], sd[j-1]);
                        if (hi > -1)
                        {
                            SPEX_MPZ_DIVEXACT(w_dense->x[i],
                                                w_dense->x[i], sd[hi]);
                        }
                    }

                    // tmpz = sigma*w[i]*w[P[j]]
                    SPEX_MPZ_MUL(tmpz, w_dense->x[i],
                                            w_dense->x[Pj]);
                    SPEX_MPZ_MUL_SI(tmpz, tmpz, sigma);
//#if 0
#ifdef SPEX_DEBUG
                    // tmpz /= sd_old[j-1]
                    mpq_t r1, r2; mpq_init(r1); mpq_init(r2);
                    mpz_fdiv_qr(tmpz, SPEX_MPQ_NUM(r1),
                                tmpz, sd0_old);
                    mpq_set_den(r1, sd0_old);
                    mpq_canonicalize(r1);

                    SPEX_MPZ_MUL(L->v[j]->x[p],
                                            L->v[j]->x[p],
                                            SPEX_MPQ_NUM(pending_scale));
                    mpz_cdiv_qr(L->v[j]->x[p], SPEX_MPQ_NUM(r2),
                                L->v[j]->x[p], SPEX_MPQ_DEN(pending_scale));
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
                    SPEX_MPZ_FDIV_Q(tmpz, tmpz, sd0_old);

                    // L(i,j) *= sd_ratio
                    SPEX_MPZ_MUL(L->v[j]->x[p],
                                            L->v[j]->x[p],
                                            SPEX_MPQ_NUM(pending_scale));
                    SPEX_MPZ_CDIV_Q(L->v[j]->x[p],
                                            L->v[j]->x[p],
                                            SPEX_MPQ_DEN(pending_scale));
#endif
                    SPEX_MPZ_ADD(L->v[j]->x[p],
                                            L->v[j]->x[p], tmpz);

                    // ---------------------------------------------------------
                    // perform IPGE to update w using updated L(i,j).
                    // ---------------------------------------------------------
                    SPEX_MPZ_SGN(&sgn, L->v[j]->x[p]);
                    if (sgn == 0)
                    {
                        // h[i] is supposed to be j-1 at this moment, but it is
                        // instead set as -(j-1)-3 = -j-2 to mark that L(i,j)
                        // is already in the nnz pattern of L(:,j)
                        h[i] = -j-2;
                        continue;
                    }

                    // w[i] = w[i]*sd_new[j]
                    SPEX_MPZ_MUL(w_dense->x[i],
                                            w_dense->x[i], sd[j]);
                    // w[i] = (w[i]-L(i,j)*w[P[j]])/sd_new[j-1]
                    SPEX_MPZ_SUBMUL(w_dense->x[i],
                                            w_dense->x[Pj], L->v[j]->x[p]);
                    if (j != 0)
                    {
                        SPEX_MPZ_DIVEXACT(w_dense->x[i],
                                            w_dense->x[i],sd[j-1]);
                    }
                }
                h[i] = j;
            }

            // iterate across all nnz in w to check if fillin should be added
            // to L->v[j]
            int64_t Lj_nz = L->v[j]->nz;
            for (p = w_top; p < w_dense->nz; p++)
            {
                i = w_dense->i[p];
                hi = h[i];
                if (P_inv[i] <= j)
                {
                    w_dense->i[p] = w_dense->i[w_top];
                    w_dense->i[w_top] = i;
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
                    SPEX_MPZ_SGN(&sgn, w_dense->x[i]);
                    if (sgn == 0) {continue;}

                    // perform history update to w[i]
                    if (hi < j-1)
                    {
                        // w[i] = w[i]*sd_new[j-1]/sd_new[h[i]]
                        SPEX_MPZ_MUL(w_dense->x[i],
                                                w_dense->x[i], sd[j-1]);
                        if (hi > -1)
                        {
                            SPEX_MPZ_DIVEXACT(w_dense->x[i],
                                                w_dense->x[i], sd[hi]);
                        }
                    }

                    // allocate extra space if needed
                    if (Lj_nz == L->v[j]->nzmax)
                    {
                        SPEX_CHECK(SPEX_vector_realloc(L->v[j],
                            SPEX_MIN(n, 2*(L->v[j]->nzmax)), option));
                    }

                    // L(i,j) = sigma*w[P[j]]*w[i]/sd_old[j-1]
                    SPEX_MPZ_MUL   (L->v[j]->x[Lj_nz],
                                               w_dense->x[Pj], w_dense->x[i]);
                    SPEX_MPZ_MUL_SI(L->v[j]->x[Lj_nz],
                                               L->v[j]->x[Lj_nz], sigma);
                    SPEX_MPZ_DIVEXACT(L->v[j]->x[Lj_nz],
                                               L->v[j]->x[Lj_nz], sd0_old);
                    // insert new entry to L
                    L->v[j]->i[Lj_nz] = i;
                    Lj_nz++;

                    // w[i] can be efficiently updated as
                    // w[i] = w[i]*sd_old[j]/sd_old[j-1]
                    SPEX_MPZ_MUL(w_dense->x[i],
                                            w_dense->x[i], sd1_old);
                    SPEX_MPZ_DIVEXACT(w_dense->x[i],
                                            w_dense->x[i], sd0_old);
                    h[i] = j;
                }
            }
            L->v[j]->nz = Lj_nz;
        }
        else // w[P[j]] == 0
        {
            SPEX_MPQ_SET(L->v[j]->scale, pending_scale);
            SPEX_MPZ_DIVEXACT(sd[j], sd[j], SPEX_MPQ_DEN(sd_ratio));
            SPEX_MPZ_MUL     (sd[j], sd[j], SPEX_MPQ_NUM(sd_ratio));
        }
    }

    //--------------------------------------------------------------------------
    // construct w from w_dense
    //--------------------------------------------------------------------------
    // reallocate w if needed
    if (w_dense->nz > w->v[0]->nzmax)
    {
        SPEX_CHECK(SPEX_vector_realloc(w->v[0], w_dense->nz, option));
    }
    int64_t w_nz = 0;
    for (p = 0; p < w_dense->nz; p++)
    {
        i = w_dense->i[p];
        SPEX_MPZ_SWAP(w->v[0]->x[w_nz], w_dense->x[i]);
        w->v[0]->i[w_nz] = i;
        w_nz++;
    }
    w->v[0]->nz = w_nz;

    SPEX_FREE_ALL;
    return SPEX_OK;
}
