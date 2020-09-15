//------------------------------------------------------------------------------
// SPEX_Util/spex_expand_mpfr_array: convert mpfr aray to mpz
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SPEX/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function converts a mpfr array of size n and precision prec to
 * an appropriate mpz array of size n. To do this, the number is multiplied by
 * the appropriate power of 10 then the gcd is found. This function allows mpfr
 * arrays to be used within SPEX LU.
 */

#define SPEX_FREE_ALL               \
    SPEX_MPFR_CLEAR(expon);         \
    SPEX_MPZ_CLEAR(temp_expon);     \
    SPEX_MPZ_CLEAR(gcd);            \
    SPEX_MPZ_CLEAR(one);            \
    SPEX_MPQ_CLEAR(temp);           \
    SPEX_matrix_free(&x3, NULL);    \

#include "spex_util_internal.h"

SPEX_info spex_expand_mpfr_array
(
    mpz_t* x_out,         // full precision mpz array
    mpfr_t* x,            // mpfr array to be expanded
    mpq_t scale,          // scaling factor used (x_out = scale*x)
    int64_t n,            // size of x
    const SPEX_options *option  // command options containing the prec
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
    mpfr_t expon; SPEX_MPFR_SET_NULL(expon);
    mpz_t temp_expon, gcd, one;
    SPEX_matrix* x3 = NULL;
    SPEX_MPZ_SET_NULL(temp_expon);
    SPEX_MPZ_SET_NULL(gcd);
    SPEX_MPZ_SET_NULL(one);
    mpq_t temp; SPEX_MPQ_SET_NULL(temp);

    uint64_t prec = SPEX_OPTION_PREC (option) ;
    mpfr_rnd_t round = SPEX_OPTION_ROUND (option) ;

    SPEX_CHECK(SPEX_mpq_init(temp));
    SPEX_CHECK(SPEX_mpfr_init2(expon, prec));
    SPEX_CHECK(SPEX_mpz_init(temp_expon));
    SPEX_CHECK(SPEX_mpz_init(gcd));
    SPEX_CHECK(SPEX_mpz_init(one));

    SPEX_CHECK (SPEX_matrix_allocate(&x3, SPEX_DENSE, SPEX_MPFR, n, 1, n,
        false, true, option));

    // expon = 2^prec (overestimate)
    SPEX_CHECK(SPEX_mpfr_ui_pow_ui(expon, 2, prec, round)) ;

    for (i = 0; i < n; i++)
    {
        // x3[i] = x[i]*2^prec
        SPEX_CHECK(SPEX_mpfr_mul(x3->x.mpfr[i], x[i], expon, round));

        // x_out[i] = x3[i]
        SPEX_CHECK(SPEX_mpfr_get_z(x_out[i], x3->x.mpfr[i], round));
    }
    SPEX_CHECK(SPEX_mpfr_get_z(temp_expon, expon, round));
    SPEX_CHECK(SPEX_mpq_set_z(scale, temp_expon));

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
        SPEX_FREE_ALL;
        SPEX_mpq_set_z(scale, one);
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

