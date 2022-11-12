//------------------------------------------------------------------------------
// SPEX_Utilities/spex_dynamic_to_CSC_mpz.c: convert SPEX_DYNAMIC_CSC matrix to
// a SPEX_CSC matrix with mpz_t entries.
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2020-2021, Jinhao Chen, Chris Lourenco,
// Erick Moreno-Centeno, Timothy A. Davis.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// spex_dynamic_to_CSC_mpz create a SPEX_matrix A in Compressed Sparse Column
// (CSC) format with entries in mpz_t type from the give matrix B, which is a
// SPEX_DYNAMIC_CSC matrix with mpz_t entries.

#define SPEX_FREE_ALL                \
    SPEX_matrix_free(&A, option);

#include "spex_util_internal.h"

SPEX_info spex_dynamic_to_CSC_mpz
(
    SPEX_matrix *A_handle,          // converted CSC matrix
    // input:
    const SPEX_matrix B,            // original matrix (not modified)
    const int64_t nnz,              // number of nonzeros in B
    const SPEX_options option
)
{
    SPEX_REQUIRE (B, SPEX_DYNAMIC_CSC, SPEX_MPZ) ;
    if (A_handle == NULL)   {return SPEX_INCORRECT_INPUT;}
    // B has been checked by the caller so no need to check here

    SPEX_info info;
    int sgn;
    int64_t i, j, Ap = 0, Bp;
    (*A_handle) = NULL;
    SPEX_matrix A = NULL;

    // allocate space for A
    SPEX_CHECK(SPEX_matrix_allocate(&A, SPEX_CSC, SPEX_MPZ, B->m, B->n, nnz,
        false, true, option));

    // initialize for A and construct A from B
    A->p[0] = 0;
    for (j = 0 ; j < B->n ; j++)
    {
        SPEX_CHECK(SPEX_mpq_cmp_ui(&sgn, B->v[j]->scale, 1, 1));
        for (Bp = 0 ; Bp < B->v[j]->nz ; Bp++)
        {
            i = B->v[j]->i[Bp];
            A->i[Ap] = i ;
            // A->x[Ap] = B->v[j]->x[Bp]*scale
            if (sgn == 0) // scale == 1
            {
                SPEX_CHECK(SPEX_mpz_set(SPEX_1D(A, Ap, mpz),
                    B->v[j]->x[Bp])) ;
            }
            else
            {
                SPEX_CHECK(SPEX_mpz_divexact(SPEX_1D(A, Ap, mpz),
                    B->v[j]->x[Bp], SPEX_MPQ_DEN(B->v[j]->scale))) ;
                SPEX_CHECK(SPEX_mpz_mul(SPEX_1D(A, Ap, mpz),
                    SPEX_1D(A, Ap, mpz), SPEX_MPQ_NUM(B->v[j]->scale))) ;
            }
            Ap++;
        }
        A->p[j+1] = Ap;
    }
    SPEX_CHECK(SPEX_mpq_set(A->scale, B->scale));

    (*A_handle) = A;
    return SPEX_OK;
}
