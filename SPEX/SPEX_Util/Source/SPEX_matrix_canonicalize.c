//------------------------------------------------------------------------------
// SPEX_Util/SPEX_matrix_canonicalize.c: canonicalize a SPEX_matrix matrix such
// that the pivot of the vector is found as the first entry of the nnz list.
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2020-2021, Jinhao Chen, Chris Lourenco (US Naval Academy),
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// SPEX_matrix_canonicalize is called to canonicalize a SPEX_matrix matrix in
// SPEX_DYNAMIC_CSC format, such that the pivot of the vector is found as the
// first entry of the nnz list.

#include "spex_util_internal.h"

SPEX_info SPEX_matrix_canonicalize
(
    SPEX_matrix *A,       // the matrix to be canonicalize
    const int64_t *perm,  // the permuation vector applied on each vector of A,
                          // considered as identity if input as NULL
    const SPEX_options *option
)
{
    SPEX_REQUIRE_KIND(A, SPEX_DYNAMIC_CSC);
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
