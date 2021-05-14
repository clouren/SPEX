//------------------------------------------------------------------------------
// SPEX_Util/SPEX_permute_row.c: permute row indices of a matrix.
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2020-2021, Jinhao Chen, Chris Lourenco (US Naval Academy),
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// SPEX_permute_row is called to permute the row indices of all columns of
// a give SPEX_matrix A.

#include "spex_util_internal.h"

SPEX_info SPEX_permute_row
(
    SPEX_matrix *A,     // input matrix
    const int64_t *perm,// desire permutation to be applied to A
    const SPEX_options *option
)
{
    SPEX_REQUIRE_KIND(A, SPEX_DYNAMIC_CSC);

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
