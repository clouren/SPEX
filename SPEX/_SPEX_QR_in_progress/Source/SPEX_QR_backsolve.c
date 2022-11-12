//------------------------------------------------------------------------------
// SPEX_QR/Source/SPEX_QR_backsolve.c: Solve exactly x = R \ Q^T b
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021, Chris Lourenco, US Naval Academy, All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* This code performs the dense R \ Q^T b. The input assumes Q^T b has already
 * been computed. Thus, it first scales the right hand side vector then does
 * back solve.  Returns x as integer.
 */

# include "spex_qr_internal.h"

SPEX_info SPEX_QR_backsolve
(
    SPEX_matrix R,        // Upper triangular matrix
    SPEX_matrix b,        // Q^T * b
    SPEX_matrix* x_handle // Solution
)
{
    SPEX_info info;
    // Indices
    int64_t i, j, k;
    // Check inputs, the number of rows/columns of R must equal the number of rows of b
    // Also, both matrices must be mpz and dense
    ASSERT( R->m == b->m);
    ASSERT( R->type == SPEX_MPZ);
    ASSERT( b->type == SPEX_MPZ);
    ASSERT( R->kind == SPEX_DENSE);
    ASSERT( b->kind == SPEX_DENSE);

    // Solution vector
    SPEX_matrix x = NULL;

    // Set x = b
    SPEX_matrix_copy(&x, SPEX_DENSE, SPEX_MPZ, b, NULL);
    // Scale x by determinant of A'*A (R(n,n))
    for (i = 0; i < x->m*x->n; i++)
    {
        SPEX_mpz_mul( x->x.mpz[i], x->x.mpz[i], SPEX_2D(R, R->n-1, R->n-1, mpz));
    }

    // Solve R x = x for each RHS vector

    // Iterate across each RHS vector
    for (k = 0; k < b->n; k++)
    {
        for (i = R->n-1; i >= 0; i--)
        {
            // Compute x[i][k]
            for (j = i+1; j < R->n; j++)
            {
                SPEX_mpz_submul( SPEX_2D(x, i, k, mpz), SPEX_2D(x, j, k, mpz), SPEX_2D(R, i, j, mpz));
            }
            SPEX_mpz_divexact( SPEX_2D(x, i, k, mpz), SPEX_2D(x, i, k, mpz), SPEX_2D(R, i, i, mpz));
        }
    }

    // Return solution as integer
    (*x_handle) = x;
    return SPEX_OK;
}
