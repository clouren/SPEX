//------------------------------------------------------------------------------
// SPEX_Util/spex_factorization_basic_check.c: perform basic check for a given
// factorization.
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Chris Lourenco
// (US Naval Academy), Erick Moreno-Centeno, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------



#include "spex_util_internal.h"

SPEX_info spex_factorization_basic_check
(
    SPEX_factorization *F // The factorization to check
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (!spex_initialized()) {return SPEX_PANIC;}

    //--------------------------------------------------------------------------
    // basic checks
    //--------------------------------------------------------------------------

 
    return (SPEX_OK) ;
}

