//------------------------------------------------------------------------------
// SPEX_Utilities/spex_create_mpz_array: create a dense mpz array
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2019-2022, Chris Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function creates an mpz array of size n.
 * Utilized internally for creating SPEX_MPZ matrices
 */

#include "spex_util_internal.h"

mpz_t *spex_create_mpz_array
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
    mpz_t *x = (mpz_t*) SPEX_calloc(n, sizeof(mpz_t));
    if (!x) {return NULL;}
    for (int64_t i = 0; i < n; i++)
    {
        #if __GNU_MP_RELEASE < 60200
        SPEX_info info =
        #endif
        SPEX_mpz_init (x [i]);
        #if __GNU_MP_RELEASE < 60200
        if (info != SPEX_OK)
        {
            // out of memory.  NOTE: This can be triggered only when using GMP
            // v6.1.2 or earlier versions. For GMP v6.2.0 or later versions,
            // there is no memory allocation, and thus such failure will never
            // occur.  As a result, this code cannot be untested by the tests
            // in SPEX/Tcov, when using GMP v6.2.0 or later.
            SPEX_MPZ_SET_NULL(x[i]);
            for (int64_t j = 0; j < i; j++)
            {
                if ( x[j] != NULL)
                {
                    SPEX_MPZ_CLEAR( x[j]);
                }
            }
            SPEX_FREE(x);
            return NULL;
        }
        #endif
    }
    return x;
}

