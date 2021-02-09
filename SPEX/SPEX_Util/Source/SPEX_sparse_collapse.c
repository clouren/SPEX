//------------------------------------------------------------------------------
// SPEX_Util/spex_sparse_collapse: shrink space required by a CSC mpz matrix
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function collapses a SPEX matrix. Essentially it shrinks the
 * size of x and i. so that they only take up the number of elements in the
 * matrix. For example if A->nzmax = 1000 but nnz(A) = 500, A->i and A->x are
 * of size 1000, so this function shrinks them to size 500.  This is only valid
 * in the factorization routines for sparse csc mpz matrices
 */

#include "spex_util_internal.h"

SPEX_info SPEX_sparse_collapse
(
    SPEX_matrix* A // matrix to be shrunk
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SPEX_REQUIRE (A, SPEX_CSC, SPEX_MPZ) ;

    //--------------------------------------------------------------------------

    int64_t anz = SPEX_matrix_nnz (A, NULL) ;

    // Shrink A->i and A->x such that they're of size anz.  These calls to
    // SPEX_realloc cannot fail since the space is shrinking.

    bool ok ;
    A->i = (int64_t *)
        SPEX_realloc (anz, A->nzmax, sizeof (int64_t), A->i, &ok) ;
    ASSERT (ok) ;

    A->x.mpz = (mpz_t *)
        SPEX_realloc (anz, A->nzmax, sizeof (mpz_t), A->x.mpz, &ok) ;
    ASSERT (ok) ;

    A->nzmax = anz ;
    return (SPEX_OK) ;
}

