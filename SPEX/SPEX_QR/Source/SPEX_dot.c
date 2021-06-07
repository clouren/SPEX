//------------------------------------------------------------------------------
// SPEX_QR/Source/SPEX_dot.c: Dense dot product
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021, Chris Lourenco, US Naval Academy, All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* This code performs a dense dot product for SPEX QR factorization
 */

# include "spex_qr_internal.h"

/* Compute the dot product of two integer vectors x,y and return in z */
SPEX_info SPEX_dot
(
    SPEX_matrix* x,
    SPEX_matrix* y,
    mpz_t z
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
    SPEX_CHECK(SPEX_mpz_set_ui(z, 0));
    for (k = 0; k < x->m; k++)
    {
        SPEX_CHECK(SPEX_mpz_addmul(z, x->x.mpz[k], y->x.mpz[k]));
    }
    return SPEX_OK;
}
