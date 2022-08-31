//------------------------------------------------------------------------------
// SPEX_Utilities/SPEX_create_default_options: set defaults
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: Create SPEX_options pointer with default parameters
 * upon successful allocation, which are defined in spex_internal.h
 */

#include "spex_util_internal.h"


SPEX_info SPEX_create_default_options (SPEX_options **option_handle )
{

    if (!spex_initialized ( )) return (SPEX_PANIC) ;

    //--------------------------------------------------------------------------
    // allocate the option struct
    //--------------------------------------------------------------------------

    *option_handle = SPEX_malloc(sizeof(SPEX_options)) ;
    if (!(*option_handle))
    {
        // out of memory
        return (SPEX_OUT_OF_MEMORY) ;
    }

    //--------------------------------------------------------------------------
    // set defaults
    //--------------------------------------------------------------------------

    (*option_handle)->pivot       = SPEX_DEFAULT_PIVOT ;
    // TODO Discuss: Erick/Tim what should we do here
    (*option_handle)->order       = SPEX_DEFAULT_ORDER ;
    (*option_handle)->print_level = SPEX_DEFAULT_PRINT_LEVEL ;
    (*option_handle)->prec        = SPEX_DEFAULT_PRECISION ;
    (*option_handle)->tol         = SPEX_DEFAULT_TOL ;
    (*option_handle)->round       = SPEX_DEFAULT_MPFR_ROUND ;
    (*option_handle)->algo        = SPEX_DEFAULT_ALGORITHM ;

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    return SPEX_OK ;
}
