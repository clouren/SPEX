//------------------------------------------------------------------------------
// SPEX_QR/Source/SPEX_Qtb.c: Compute dense Q'*b
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021-2022, Chris Lourenco, Lorena Mejia Domenzain,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* This code performs the dense Q'*b
 */

# include "spex_qr_internal.h"

SPEX_info SPEX_Qtb
(
    SPEX_matrix Q,        // Q matrix, want Q'
    SPEX_matrix b,        // Original RHS Vector
    SPEX_matrix* b_handle // Null on input. Contains Q'*b on output
)
{
    SPEX_info info;
    // Check inputs, the number of columns of Q' must equal the number of rows of b
    // Also, both matrices must be mpz and dense
    ASSERT( Q->m == b->m);
    ASSERT( Q->type == SPEX_MPZ);
    ASSERT( b->type == SPEX_MPZ);
    ASSERT( Q->kind == SPEX_DENSE);
    ASSERT( b->kind == SPEX_DENSE);

    SPEX_matrix b_new = NULL;

    // b->new has Q->n rows and b->n columns
    SPEX_CHECK(SPEX_matrix_allocate(&b_new, SPEX_DENSE, SPEX_MPZ, Q->n, b->n, Q->n*b->n,
        false, true, NULL));

    // Indices
    int64_t i, j, k;
    // Need to compute b_new[i] = Q'[i,:] dot b[i]
    // This is equivalent to b_new[i] = Q[:,i] dot b[i]

    // Iterate across every RHS vector
    for (k = 0; k < b->n; k++)
    {
        // Compute b[i,k]
        for (i = 0; i < Q->n; i++)
        {
            for (j = 0; j < Q->m; j++)
            {
                SPEX_CHECK( SPEX_mpz_addmul( SPEX_2D(b_new,i,k,mpz),
                                             SPEX_2D(b, j, k, mpz),
                                             SPEX_2D(Q, j, i, mpz)));
            }
        }
    }

    (*b_handle) = b_new;
    return SPEX_OK;
}
