//------------------------------------------------------------------------------
// SPEX_LU_Update/SPEX_mat_canonicalize.c: canonicalize a SPEX_mat matrix such
// that the pivot of the vector is found as the first entry of the nnz list.
//------------------------------------------------------------------------------

// SPEX_LU_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_LU_Update/License for the license.

//------------------------------------------------------------------------------

// SPEX_mat_canonicalize is called to canonicalize a SPEX_mat matrix such
// that the pivot of the vector is found as the first entry of the nnz list.

#include "spex_lu_update_internal.h"

SPEX_info SPEX_mat_canonicalize
(
    SPEX_mat *A,    // the matrix to be canonicalize
    int64_t *perm   // the permuation vector applied on each vector of A,
                    // considered as identity if input as NULL
)
{
    SPEX_info info;
    int64_t i, j, p, diag;
    for (j = 0; j < A->n; j++)
    {
        diag = (perm == NULL) ? j : perm[j];
        for (p = 0; p < A->v[j]->nz; p++)
        {
            i = A->v[j]->i[p];
            if (i == diag)
            {
                if (p != 0)
                {
                    SPEX_CHECK(SPEX_mpz_swap(A->v[j]->x[0], A->v[j]->x[p]));
                    A->v[j]->i[p] = A->v[j]->i[0];
                    A->v[j]->i[0] = i;
                }
                break;
            }
        }
    }
    return SPEX_OK;
}
