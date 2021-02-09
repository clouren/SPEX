//------------------------------------------------------------------------------
// SPEX_Util/spex_expand_mpq_array: convert mpq array to mpz
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function converts a mpq array of size n into an appropriate
 * mpz array of size n. To do this, the lcm of the denominators is found as a
 * scaling factor. This function allows mpq arrays to be used in SPEX LU.
 */

#define SPEX_FREE_ALL               \
    SPEX_MPZ_CLEAR(temp);           \
    SPEX_matrix_free(&x3, NULL);    \
    SPEX_matrix_free(&x4, NULL);

#include "spex_util_internal.h"

SPEX_info spex_expand_mpq_array
(
    mpz_t* x_out,        // mpz array, on output x_out = x*scale
    mpq_t* x,            // mpq array that needs to be converted
    mpq_t scale,         // scaling factor. x_out = scale*x
    int64_t n,           // size of x
    const SPEX_options* option // Command options
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------
    // inputs have checked in the only caller spex_cast_array
    ASSERT(n >= 0);    
    SPEX_info info ;

    //--------------------------------------------------------------------------

    // Define temporary matrices
    mpz_t temp;
    SPEX_matrix *x3 = NULL;
    SPEX_matrix *x4 = NULL;;
    SPEX_MPZ_SET_NULL(temp);
    SPEX_CHECK (SPEX_mpz_init(temp)) ;

    // Allocate memory
    SPEX_CHECK (SPEX_matrix_allocate(&x3, SPEX_DENSE, SPEX_MPZ, n, 1, n,
        false, true, option));
    SPEX_CHECK (SPEX_matrix_allocate(&x4, SPEX_DENSE, SPEX_MPQ, n, 1, n,
        false, true, option));

    ASSERT( x3 != NULL);
    ASSERT( x4 != NULL);
    // x3 = denominators of (mpq_t) x
    for (int64_t i = 0; i < n; i++)
    {
        SPEX_CHECK(SPEX_mpq_get_den(x3->x.mpz[i], x[i]));
    }

    // Find LCM of denominators of x
    SPEX_CHECK(SPEX_mpz_set(temp,x3->x.mpz[0]));
    for (int64_t i = 1; i < n; i++)
    {
        SPEX_CHECK(SPEX_mpz_lcm(temp, x3->x.mpz[i], temp));
    }
    SPEX_CHECK(SPEX_mpq_set_z(scale,temp));

    for (int64_t i = 0; i < n; i++)
    {
        // x4[i] = x[i]*temp
        SPEX_CHECK(SPEX_mpq_mul(x4->x.mpq[i], x[i], scale));

        // x_out[i] = x4[i]
        SPEX_CHECK(SPEX_mpz_set_q(x_out[i], x4->x.mpq[i]));
    }
    SPEX_FREE_ALL;
    return SPEX_OK;
}

