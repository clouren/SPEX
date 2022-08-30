//------------------------------------------------------------------------------
// SPEX_Utilities/spex_CSC_mpz_to_dynamic.c: convert a SPEX_matrix of CSC x MPZ to
// a dynamic_CSC matrix.
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2020-2021, Jinhao Chen, Chris Lourenco (US Naval Academy),
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// spex_CSC_mpz_to_dynamic create a SPEX_matrix A of dynamic_CSC from the given
// SPEX_matrix B which is in Compressed Sparse Column (CSC) format with entries
// in mpz_t type.

#define SPEX_FREE_ALL                \
    SPEX_matrix_free(&A, option);

#include "spex_util_internal.h"

SPEX_info spex_CSC_mpz_to_dynamic
(
    SPEX_matrix **A_handle,      // converted SPEX_matrix of dynamic_CSC
    // input:
    const SPEX_matrix *B,         // original matrix (unmodified)
    const SPEX_options *option
)
{
    SPEX_REQUIRE (B, SPEX_CSC, SPEX_MPZ) ;
    if (A_handle == NULL)   {return SPEX_INCORRECT_INPUT;}

    SPEX_info info;
    (*A_handle) = NULL;
    SPEX_matrix *A = NULL;
    int64_t *Bp = B->p;

    int64_t i, j, p, Ap = 0;

    // allocate space for A
    SPEX_CHECK(SPEX_matrix_allocate(&A, SPEX_DYNAMIC_CSC, SPEX_MPZ, B->m, B->n,
        0, false, true, option));

    for (j = 0; j < B->n; j++)
    {
        // reallocate space for each column of A
        SPEX_CHECK(SPEX_vector_realloc(A->v[j], Bp[j+1]-Bp[j], option));
        Ap = 0;
        for (p = Bp[j]; p < Bp[j+1]; p++)
        {
            i = B->i[p];
            A->v[j]->i[Ap] = i;
            // A->v[j]->x[Ap] = B->x[p]
            SPEX_CHECK(SPEX_mpz_set(A->v[j]->x[Ap], SPEX_1D(B, p, mpz)));
            Ap++;
        }
        A->v[j]->nz = Ap;
    }
    SPEX_CHECK(SPEX_mpq_set(A->scale, B->scale));

    (*A_handle) = A;
    return SPEX_OK;
}
