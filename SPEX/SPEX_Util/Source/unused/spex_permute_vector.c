//------------------------------------------------------------------------------
// SPEX_Util/spex_permute_dense_matrix_row: permute rows of a dense matrix A,
// as A = P*A
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function permutes a dense matrix A as A = P*A.
copied from SPEX_Left_LU/Sourse/spex_left_lu_permute_b.c
 */

#define SPEX_FREE_ALL        \
    SPEX_matrix_free (&Atmp, NULL) ;

#include "spex_util_internal.h"

SPEX_info spex_permute_dense_matrix_row
(
    SPEX_matrix **A_handle,     // permuted A
    const SPEX_matrix *A,       // unpermuted A (not modified)
    const int64_t *P,           // row permutation
    const SPEX_options* option
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SPEX_info info ;
    SPEX_REQUIRE (A, SPEX_DENSE, SPEX_MPZ) ;

    if (A_handle == NULL || !P) {return SPEX_INCORRECT_INPUT;}
    (*A_handle) = NULL ;

    //--------------------------------------------------------------------------
    // b(pinv) = b2
    //--------------------------------------------------------------------------

    int64_t m = A->m ;
    int64_t n = A->n ;

    // allocate x
    SPEX_matrix *Atmp = NULL ;
    SPEX_CHECK (SPEX_matrix_allocate (&Atmp, SPEX_DENSE, SPEX_MPZ, m, n,
        0, false, true, option)) ;

    // Set Atmp = P*A
    for (int64_t i = 0 ; i < m ; i++)
    {
        for (int64_t j = 0 ; j < n ; j++)
        {
            SPEX_CHECK(SPEX_mpz_set(SPEX_2D(Atmp,  P[i], j, mpz),
                                    SPEX_2D(A,      i,   j, mpz)));
        }
    }

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    (*A_handle) = Atmp ;
    return SPEX_OK;
}

