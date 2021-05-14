//------------------------------------------------------------------------------
// SPEX_Update/spex_update_get_scattered_v.c: build scattered vector for given
// sparse vector
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is called to build scattered mpz vector for column or
// row k of A, or the inserted column. This function eliminates explicit
// 0. If keep_v is false, this function swap the mpz values, and
// thus the original vector will become all zeros. Otherwise, mpz_set will be
// used to keep the original mpz values.

#define SPEX_FREE_ALL   \
    spex_scattered_vector_free(&sv, NULL);

#include "spex_update_internal.h"

SPEX_info spex_update_get_scattered_v
(
    spex_scattered_vector **sv_handle,// output vector in scattered form
    SPEX_vector *v,              // the vector in compressed form, whose
                                 // max index is n
    const int64_t n,             // number of entries in v
    const bool keep_v            // indicate if the mpz values should be kept
)
{
    if (!sv_handle || !v)
    {
        return SPEX_INCORRECT_INPUT;
    }
    *sv_handle = NULL;
    SPEX_info info;
    int64_t p, i;
    spex_scattered_vector *sv = NULL;

    // make sv a sparse vector so that we can have sv->i for nnz pattern
    SPEX_CHECK(spex_scattered_vector_alloc(&sv, n, true, NULL));

    int sgn;
    p = 0;
    while (p < v->nz)
    {
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, v->x[p]));
        if (sgn != 0)
        {
            i = v->i[p];
            /* TODO combine spex_find_next_nz with this function?
            if (k != perm_inv[i] && perm_inv[i] < inext)
            {
                inext = perm_inv[i];
            }
            */
            if (keep_v)
            {
                SPEX_CHECK(SPEX_mpz_set (sv->x[i], v->x[p]));
            }
            else
            {
                SPEX_CHECK(SPEX_mpz_swap(sv->x[i], v->x[p]));
            }
            sv->i[p] = i;
            p++;
        }
        else
        {
            v->nz--;
            SPEX_CHECK(SPEX_mpz_swap(v->x[p], v->x[v->nz]));
            v->i[p] = v->i[v->nz];
        }
    }
    sv->nz = v->nz;

    *sv_handle = sv;
    return SPEX_OK;
}

