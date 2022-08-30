//------------------------------------------------------------------------------
// SPEX_Utilities/spex_expand_double_array: convert double vector to mpz
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function converts a double array of size n to an appropriate
 * mpz array of size n. To do this, the number is multiplied by an appropriate
 * power, then, the GCD is found. This function allows the use of matrices in
 * double precision to work with SPEX.
 *
 */

#define SPEX_FREE_WORKSPACE         \
    SPEX_MPZ_CLEAR(gcd);            \
    SPEX_MPZ_CLEAR(one);            \
    SPEX_MPQ_CLEAR(temp);           \
    SPEX_matrix_free(&x3, NULL);    \

#define SPEX_FREE_ALL              \
    SPEX_FREE_WORKSPACE            \

#include "spex_util_internal.h"

SPEX_info spex_expand_double_array
(
    mpz_t* x_out,           // integral final array
    double* x,              // double array that needs to be made integral
    mpq_t scale,            // the scaling factor used (x_out = scale * x)
    int64_t n,              // size of x
    const SPEX_options* option // Command options
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    // inputs have been checked in the only caller spex_cast_array

    //--------------------------------------------------------------------------

    int64_t i, k ;
    int r1, r2 = 1;
    bool nz_found = false;
    SPEX_info info ;
    // Double precision accurate to about 2e-16. We multiply by 1e16 to convert
    // Note that this conversion allows a number like 0.9 to be represented exactly
    // fl(0.9) is not exact in double precision; in fact the exact conversion is
    // fl(0.9) = 45000000000000001 / 50000000000000000. Multiplying by 1e16 gives
    // the actual value of 9/10 when scaled. Note that if this type of conversion
    // is not desired by the user it is suggested they first convert from double to
    // MPFR then from MPFR to MPQ as that will be fully exact.
    double expon = pow(10, 16);
    // Convert the input x into a quad precision matrix. This is to handle the (rare)
    // case that the user gives an input double which is close to DOUBLE_MAX. In that case
    // the multiplication could lead to inf.
    SPEX_matrix* x3 = NULL;
    mpz_t gcd, one; SPEX_MPZ_SET_NULL(gcd); SPEX_MPZ_SET_NULL(one);
    mpq_t temp; SPEX_MPQ_SET_NULL(temp);

    mpfr_rnd_t round = SPEX_OPTION_ROUND (option) ;

    SPEX_CHECK(SPEX_mpq_init(temp));
    SPEX_CHECK(SPEX_mpz_init(gcd));
    SPEX_CHECK(SPEX_mpz_init(one));

    SPEX_CHECK (SPEX_matrix_allocate(&x3, SPEX_DENSE, SPEX_MPFR, n, 1, n,
        false, true, option));

    SPEX_CHECK(SPEX_mpq_set_d(scale, expon));           // scale = 10^16
    for (i = 0; i < n; i++)
    {
        // x3[i] = x[i], cast double to MPFR
        SPEX_CHECK(SPEX_mpfr_set_d(x3->x.mpfr[i], x[i], round));

        // x3[i] = x[i] * 10^16, multiply MPFR by 10^16
        SPEX_CHECK(SPEX_mpfr_mul_d(x3->x.mpfr[i], x3->x.mpfr[i], expon, round));

        // x_out[i] = x3[i], cast MPFR to integer truncating remaining decimal component
        SPEX_CHECK(SPEX_mpfr_get_z(x_out[i], x3->x.mpfr[i], round));
    }

    //--------------------------------------------------------------------------
    // Compute the GCD to reduce the size of scale
    //--------------------------------------------------------------------------

    SPEX_CHECK(SPEX_mpz_set_ui(one, 1));
    // Find an initial GCD
    for (i = 0; i < n; i++)
    {
        if (!nz_found)
        {
            SPEX_CHECK(SPEX_mpz_cmp_ui(&r1, x_out[i], 0)); // Check if x[i] == 0
            if (r1 != 0)
            {
                nz_found = true;
                k = i;
                SPEX_CHECK(SPEX_mpz_set(gcd, x_out[i]));
            }
        }
        else
        {
            // Compute the GCD, stop if gcd == 1
            SPEX_CHECK(SPEX_mpz_gcd(gcd, gcd, x_out[i]));
            SPEX_CHECK(SPEX_mpz_cmp(&r2, gcd, one));
            if (r2 == 0)
            {
                break;
            }
        }
    }

    if (!nz_found)     // Array is all zeros
    {
        SPEX_mpq_set_z(scale, one);
        SPEX_FREE_ALL;
        return SPEX_OK;
    }

    //--------------------------------------------------------------------------
    // Scale all entries to make as small as possible
    //--------------------------------------------------------------------------

    if (r2 != 0)             // If gcd == 1 then stop
    {
        for (i = k; i < n; i++)
        {
            SPEX_CHECK(SPEX_mpz_divexact(x_out[i], x_out[i], gcd));
        }
        SPEX_CHECK(SPEX_mpq_set_z(temp, gcd));
        SPEX_CHECK(SPEX_mpq_div(scale, scale, temp));
    }
    SPEX_FREE_WORKSPACE;
    return SPEX_OK;
}

