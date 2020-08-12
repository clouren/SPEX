//------------------------------------------------------------------------------
// SPEX_Util/spex_create_mpz_array: create a dense mpz array
//------------------------------------------------------------------------------

// SPEX: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SPEX/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function creates an mpz array of size n.
 * Utilized internally for creating SPEX_MPZ matrices
 */

#include "spex_util_internal.h"

mpz_t* spex_create_mpz_array
(
    int64_t n            // size of the array
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (n <= 0) {return NULL;}

    //--------------------------------------------------------------------------

    // Malloc space
    mpz_t* x = (mpz_t*) SPEX_calloc(n, sizeof(mpz_t));
    if (!x) {return NULL;}
    for (int64_t i = 0; i < n; i++)
    {
        if (SPEX_mpz_init(x[i]) != SPEX_OK)
        {
            // Out of memory
            SPEX_MPZ_SET_NULL(x[i]);
            for (int64_t j = 0; j < n; j++)
            {
                if ( x[j] != NULL)
                {
                    SPEX_MPZ_CLEAR( x[j]);
                }
            }
            SPEX_FREE(x);
            return NULL;
        }
    }
    return x;
}

