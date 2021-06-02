//------------------------------------------------------------------------------
// SPEX_QR/Source/SPEX_dense_matrix_dot.c: Dense dot product
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021, Chris Lourenco, US Naval Academy, All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* This code performs a dense matrix dot product for SPEX QR factorization
 */

# include "spex_qr_internal.h"

/* Purpose: Given a matrix A in m*n and B in m*n, compute the dot product of
 * A(:,i) and B(:,j). Assumed to be dense. prod = A(:,i) dot B(:,j)
 */
SPEX_info SPEX_dense_mat_dot
(
    SPEX_matrix* A,
    int64_t i,
    SPEX_matrix* B,
    int64_t j,
    mpz_t prod
)
{
    SPEX_info info;
    ASSERT(A->m == B->m);
    for (int64_t k = 0; k < A->m; k++)
    {
        SPEX_CHECK(SPEX_mpz_addmul(prod, SPEX_2D(A, k, i, mpz),
                        SPEX_2D(B, k, j, mpz)));
    }
    return SPEX_OK;
}
