//------------------------------------------------------------------------------
// SPEX_Utilities/spex_delete_mpz_array: clear the memory used for a mpz array
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2020-2023, Chris Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function clears the memory used for an mpz vector of size n.
 * Call this function for all mpz vectors when done.
 *
 * Input is a mpz_t** array which is destroyed upon function completion
 *
 */

#include "spex_util_internal.h"

void spex_delete_mpz_array
(
    mpz_t **x,      // mpz array to be deleted
    int64_t n       // Size of x
)
{

    if (x == NULL || (*x) == NULL) {return;}
    for (int64_t i = 0; i < n; i++)
    {
        if ( (*x) [i] != NULL)
        {
            SPEX_MPZ_CLEAR((*x) [i]);
        }
    }
    SPEX_FREE ((*x));
}
