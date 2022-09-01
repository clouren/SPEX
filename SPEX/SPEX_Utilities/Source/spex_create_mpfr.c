//------------------------------------------------------------------------------
// SPEX_Utilities/spex_create_mpfr: create an mpfr_t entry
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function safely creates and initializes an mpfr_t entry.
 */

#include "spex_util_internal.h"

SPEX_info spex_create_mpfr
(
    mpfr_t x,                  // mpfr_t entry to be initialized
    const SPEX_options *option // command options containing the prec for mpfr
)
{
    uint64_t prec = SPEX_OPTION_PREC (option) ;
    SPEX_info info = SPEX_mpfr_init2(x, prec);
    if (info != SPEX_OK)
    {
        // Out of memory
        SPEX_MPFR_SET_NULL(x);
        return info;
    }
    return SPEX_OK;
}

