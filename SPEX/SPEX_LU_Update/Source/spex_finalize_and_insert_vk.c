//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_finalize_and_insert_vk.c: perform history update for
// entries that would be in L and insert entries that would in U to
// corresponding row of U.
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is called to perform history update for entries that
// would be in L and insert entries that would in U to corresponding row of U.
// Entries that would appear in U are removed from the nnz pattern, while
// entries that would be in L are updated and moved from vk_dense. The pivot
// entry is kept as the first entry in L->v[k] but not inserted to U.

#include "spex_lu_update_internal.h"

#define SU(k) S->x.mpq[1+2*(k)]

SPEX_info spex_finalize_and_insert_vk
(
    spex_scattered_vector *vk_dense, //scattered version of the solution for
                      // LDx=v using the first k-1 columns of L
    int64_t *h,       // history vector for vk_dense
    SPEX_mat *U,      // matrix U
    SPEX_mat *L,      // matrix L
    SPEX_matrix *S,   // a 2*n dense mpq matrix that stores pending scales
    const mpz_t *sd,  // array of scaled pivots
    const int64_t *Q, // the column permutation
    const int64_t *P_inv,// inverse of row permutation
    const int64_t k,  // the column index in L that vk_dense will be inserted
    const int64_t diag,// the index of entry in vk_dense that will be diagonal
    const mpq_t one   // a const mpq number, just to avoid constantly alloc
)
{
    SPEX_info info;
    int64_t i, p = 0, real_i, vk_nz = vk_dense->nz, Lk_nz;
    int sgn;

    // move entries to U
    while(p < vk_nz)
    {
        i = vk_dense->i[p];
        real_i = P_inv[i];
        if (real_i < k)
        {
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, vk_dense->x[i]));
            if (sgn != 0)
            {
                // insert vk_dense->x[i] to U(P_inv[i],Q[k]) by swapping
                SPEX_CHECK(spex_insert_new_entry(vk_dense->x[i], U->v[real_i],
                    SU(real_i), Q[k], one));
                // vk_dense->x[i] = 0
                SPEX_CHECK(SPEX_mpz_set_ui(vk_dense->x[i], 0));
            }
            vk_nz--;
            vk_dense->i[p] = vk_dense->i[vk_nz];
        }
        else
        {
            p++;
        }
    }

    // check if L->v[k] needs more space for all remaining entries
    if (vk_nz > L->v[k]->nzmax)
    {
        SPEX_CHECK(SPEX_vector_realloc(L->v[k], vk_nz));
    }

    // move the remaining nonzero entries to L->v[k]
    Lk_nz = 1; // reserve the first entry for the pivot
    for (p = 0; p < vk_nz; p++)
    {
        i = vk_dense->i[p];
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, vk_dense->x[i]));
        if (sgn == 0)        {continue;}

        h[i] = SPEX_FLIP(h[i]);
        if (h[i] < k-1)
        {
            // perform history update if needed
            SPEX_CHECK(SPEX_mpz_mul(vk_dense->x[i], vk_dense->x[i], sd[k-1]));
            if (h[i] > -1)
            {
                SPEX_CHECK(SPEX_mpz_divexact(vk_dense->x[i],
                                             vk_dense->x[i], sd[h[i]]));
            }
        }
        if (P_inv[i] == diag) // put pivot as the first entry
        {
            SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[0], vk_dense->x[i]));
            L->v[k]->i[0] = i;
        }
        else
        {
            SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[Lk_nz], vk_dense->x[i]));
            L->v[k]->i[Lk_nz] = i;
            Lk_nz++;
        }
        // vk_dense->x[i] = 0
        SPEX_CHECK(SPEX_mpz_set_ui(vk_dense->x[i], 0));
    }
    L->v[k]->nz = Lk_nz;

    return SPEX_OK;
}
