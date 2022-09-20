//------------------------------------------------------------------------------
// SPEX_Utilities/spex_create_mpz: create an mpz_t entry
//------------------------------------------------------------------------------

// FIXME: delete this file
#if 0

// SPEX_Utilities: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function safely creates and initializes an mpz_t entry.
 */

#include "spex_util_internal.h"

SPEX_info spex_create_mpz
(
    mpz_t x                  // mpz_t entry to be initialized
)
{
    #if __GNU_MP_RELEASE < 60200
    SPEX_info info =
    #endif
    SPEX_mpz_init(x);
    #if __GNU_MP_RELEASE < 60200
    if (info != SPEX_OK)
    {
        // Out of memory  NOTE: This can be triggered only when using GMP
        // v6.1.2 or earlier versions. For GMP v6.2.0 or later versions,
        // there is no memory allocation, and thus such failure will never
        // occur.  As a result, this code cannot be untested by the tests
        // in SPEX/Tcov, when using GMP v6.2.0 or later.
        SPEX_MPZ_SET_NULL(x);
        return info;
    }
    #endif
    return SPEX_OK;
}

#endif
