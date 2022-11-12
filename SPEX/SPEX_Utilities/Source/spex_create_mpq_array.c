//------------------------------------------------------------------------------
// SPEX_Utilities/spex_create_mpq_array: create a dense mpq array
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2019-2021, Chris Lourenco, Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function creates an mpq array of size n.
 * This function is used internally for creating SPEX_MPQ matrices
 */

#include "spex_util_internal.h"

mpq_t *spex_create_mpq_array
(
    int64_t n              // size of the array
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (n <= 0) {return NULL;}

    //--------------------------------------------------------------------------

    // calloc space
    mpq_t *x = (mpq_t*) SPEX_calloc(n, sizeof(mpq_t));
    if (!x) {return NULL;}
    for (int64_t i = 0; i < n; i++)
    {
        if (SPEX_mpq_init(x[i]) != SPEX_OK)
        {
            // Out of memory
            SPEX_MPQ_SET_NULL(x[i]);
            for (int64_t j = 0; j < i; j++)
            {
                SPEX_MPQ_CLEAR( x[j]);
            }
            SPEX_FREE(x);
            return NULL;
        }
    }
    return x;
}

