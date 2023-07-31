//------------------------------------------------------------------------------
// SPEX_QR/Source/SPEX_dot.c: Dense dot product
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021-2023, Chris Lourenco, Lorena Mejia Domenzain,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* This code performs a dense dot product for SPEX QR factorization
 */

# include "spex_qr_internal.h"

/* Compute the dot product of two integer vectors x,y and return in z */
SPEX_info SPEX_dot
(
    SPEX_matrix x, // First vector
    SPEX_matrix y, // Second vector
    mpz_t z         // Dot product of x and y. Should be initialized on input
)
{
    SPEX_info info;
    // Check inputs, x and y must both be column vectors of identical size and
    // stored as dense
    ASSERT( x->n == y->n);
    ASSERT( x->m == y->m);
    ASSERT( x->type == SPEX_MPZ);
    ASSERT( y->type == SPEX_MPZ);
    ASSERT( x->kind == SPEX_DENSE);
    ASSERT( y->kind == SPEX_DENSE);

    int64_t k;

    // Set z = 0
    SPEX_MPZ_SET_UI(z, 0);
    for (k = 0; k < x->m; k++)
    {
        // z += x[k]*y[k]
        SPEX_MPZ_ADDMUL(z, x->x.mpz[k], y->x.mpz[k]);
    }
    return SPEX_OK;
}