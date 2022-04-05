//------------------------------------------------------------------------------
// SPEX_Update/spex_update_get_scattered_v.c: build scattered vector for given
// sparse vector
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is called to build scattered mpz vector for column or
// row k of A, or the inserted column. This function eliminates explicit 0
// regardless of keep_v. If keep_v is false, this function swap the mpz values,
// and thus the original vector will become all zeros. Otherwise, mpz_set will
// be used to make a copy of the original mpz values.
//
// This function also searches for the smallest column/row index of the nonzero
// entry with pivot excluded, if requested.

#define SPEX_FREE_ALL   \
    spex_scattered_vector_free(&sv, option);

#include "spex_update_internal.h"

SPEX_info spex_update_get_scattered_v
(
    // output
    spex_scattered_vector **sv_handle,// output vector in scattered form
    int64_t *next,               // the smallest col/row index of non-pivot nz.
                                 // If next == NULL, searching is not performed.
    // input
    SPEX_vector *v,              // the vector in compressed form, whose
                                 // max index is n
    const int64_t n,             // number of entries in v
    const int64_t k,             // column index of v in L, or row index in U.
                                 // Ignored if next == NULL.
    const int64_t *perm_inv,     // inverse of permutation applied on v.
                                 // This can be NULL if next == NULL.
    const bool keep_v,           // indicate if the mpz values should be kept
    const SPEX_options *option
)
{
    *sv_handle = NULL;
    SPEX_info info;
    int64_t p, i;
    spex_scattered_vector *sv = NULL;

    // make sv a sparse vector so that we can have sv->i for nnz pattern
    SPEX_CHECK(spex_scattered_vector_alloc(&sv, n, option));

    int sgn;
    if (next != NULL) {(*next) = n;}
    p = 0;
    while (p < v->nz)
    {
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, v->x[p]));
        if (sgn != 0)
        {
            i = v->i[p];
            if (next != NULL)
            {
                // search the smallest col/row index of the non-pivot nonzero
                int64_t real_i = perm_inv[i];
                if (k != real_i && real_i < (*next))
                {
                    (*next) = real_i;
                }
            }

            if (keep_v)
            {
                // make a copy of the mpz value
                SPEX_CHECK(SPEX_mpz_set (sv->x[i], v->x[p]));
            }
            else
            {
                // swapping mpz pointer, which is more efficient
                SPEX_CHECK(SPEX_mpz_swap(sv->x[i], v->x[p]));
            }
            sv->i[p] = i;
            p++;
        }
        else
        {
            // remove explicit zeros
            v->nz--;
            SPEX_CHECK(SPEX_mpz_swap(v->x[p], v->x[v->nz]));
            v->i[p] = v->i[v->nz];
        }
    }
    sv->nz = v->nz;

    *sv_handle = sv;
    return SPEX_OK;
}

