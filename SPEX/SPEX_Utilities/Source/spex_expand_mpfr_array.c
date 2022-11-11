//------------------------------------------------------------------------------
// SPEX_Utilities/spex_expand_mpfr_array: convert mpfr aray to mpz
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function converts a mpfr array of size n and precision prec to
 * an appropriate mpz array of size n. To do this, the number is multiplied by
 * the appropriate power of 10 then the gcd is found. This function allows mpfr
 * arrays to be used within SPEX.
 */

#define SPEX_FREE_ALL               \
    SPEX_MPZ_CLEAR(gcd);            \
    SPEX_MPZ_CLEAR(one);            \
    SPEX_MPQ_CLEAR(temp);           \
    if (x_mpq)                      \
    {                               \
        for (i = 0; i < n; i++)     \
        {                           \
            if ( x_mpq[i] != NULL)  \
            {                       \
                SPEX_MPQ_CLEAR(x_mpq[i]);\
            }                       \
        }                           \
    }                               \
    SPEX_FREE (x_mpq);

#include "spex_util_internal.h"

SPEX_info spex_expand_mpfr_array
(
    mpz_t* x_out,         // full precision mpz array
    mpfr_t* x,            // mpfr array to be expanded
    mpq_t scale,          // scaling factor used (x_out = scale*x)
    int64_t n,            // size of x
    const SPEX_options option  // command options containing the prec
                          // and rounding for mpfr
)
{

    //--------------------------------------------------------------------------
    // Input has already been checked
    //--------------------------------------------------------------------------
    ASSERT(n >= 0);
    SPEX_info info ;

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    int64_t i, k ;
    int r1, r2 = 1 ;
    bool nz_found = false;
    mpz_t gcd, one;
    mpq_t* x_mpq = NULL;
    SPEX_MPZ_SET_NULL(gcd);
    SPEX_MPZ_SET_NULL(one);
    mpq_t temp; SPEX_MPQ_SET_NULL(temp);

    SPEX_CHECK(SPEX_mpq_init(temp));
    SPEX_CHECK(SPEX_mpz_init(gcd));
    SPEX_CHECK(SPEX_mpz_init(one));

    x_mpq = spex_create_mpq_array (n) ;
    if (x_mpq == NULL)
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }

    SPEX_CHECK(spex_cast_array(x_mpq, SPEX_MPQ,   x  , SPEX_MPFR,n, NULL,  NULL,
        option));
    SPEX_CHECK(spex_cast_array(x_out, SPEX_MPZ, x_mpq, SPEX_MPQ, n, scale, NULL,
        option));

    //--------------------------------------------------------------------------
    // Find the gcd to reduce scale
    //--------------------------------------------------------------------------

    SPEX_CHECK(SPEX_mpz_set_ui(one, 1));
    // Find an initial GCD
    for (i = 0; i < n; i++)
    {
        if (!nz_found)
        {
            SPEX_CHECK(SPEX_mpz_cmp_ui(&r1, x_out[i], 0));
            if (r1 != 0)
            {
                nz_found = true;
                k = i;
                SPEX_CHECK(SPEX_mpz_set(gcd, x_out[i]));
            }
        }
        else
        {
            // Compute the GCD of the numbers, stop if gcd == 1
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

    if (r2 != 0)  // If gcd == 1 stop
    {
        for (i = k; i < n; i++)
        {
            SPEX_CHECK(SPEX_mpz_divexact(x_out[i],x_out[i],gcd));
        }
        SPEX_CHECK(SPEX_mpq_set_z(temp,gcd));
        SPEX_CHECK(SPEX_mpq_div(scale,scale,temp));
    }
    SPEX_FREE_ALL;
    return SPEX_OK;
}

