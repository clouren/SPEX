//------------------------------------------------------------------------------
// SPEX_Update/spex_update_find_next_nz.c: find the next index of next nz in
// given vector
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is to find the row index of first off diagonal of
//          L(P,k) or the column index of the first off diagonal of U(k,Q)


#include "spex_update_internal.h"

SPEX_info spex_update_find_next_nz
(
    int64_t *next,                  // the col/row index of next nnz
    spex_scattered_vector *Ak_dense,// the scattered vector
    int64_t *perm_inv,              // inverse of permutation
    int64_t k
)
{
    SPEX_info info;
    int sgn;
    int64_t p, i;
    *next = Ak_dense->nzmax;

    for (p = 0; p < Ak_dense->nz; p++)
    {
        i = Ak_dense->i[p];
        int64_t real_i = perm_inv[i];
        if (k == real_i) { continue; }

        SPEX_CHECK(SPEX_mpz_sgn(&sgn, Ak_dense->x[i]));
        if (sgn != 0 && real_i < *next)
        {
            *next = real_i;
        }
    }
    return SPEX_OK;
}
