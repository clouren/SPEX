//------------------------------------------------------------------------------
// SPEX_Update/spex_update_permute_row.c: permute row indices of a matrix.
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Chris Lourenco (US Naval Academy),
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// spex_update_permute_row is called to permute the row indices of all columns
// of a given SPEX_matrix A.

#include "spex_update_internal.h"

SPEX_info spex_update_permute_row
(
    SPEX_matrix *A,     // input matrix
    const int64_t *perm // desire permutation to be applied to A, must be
                        // non-NULL
)
{
    // check inputs
    if (!spex_initialized()) {return SPEX_PANIC;}
    SPEX_REQUIRE_KIND(A, SPEX_DYNAMIC_CSC);
    if (perm == NULL) {return SPEX_INCORRECT_INPUT;}

    // permute each vector of A
    for (int64_t j = 0; j < A->n; j++)
    {
        for (int64_t p = 0; p < A->v[j]->nz; p++)
        {
            int64_t i = A->v[j]->i[p];
            A->v[j]->i[p] = perm[i];
        }
    }

    return SPEX_OK;
}
