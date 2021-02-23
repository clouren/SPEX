//------------------------------------------------------------------------------
// SPEX_Util/SPEX_create_default_options: set defaults
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: Create and return SPEX_options pointer with default parameters
 * upon successful allocation, which are defined in spex_internal.h
 */

#include "spex_util_internal.h"

// TODO redesign this API

SPEX_options* SPEX_create_default_options ( void )
{

    if (!spex_initialized ( )) return (NULL) ;

    //--------------------------------------------------------------------------
    // allocate the option struct
    //--------------------------------------------------------------------------

    SPEX_options* option = SPEX_malloc(sizeof(SPEX_options)) ;
    if (!option)
    {
        // out of memory
        return (NULL) ;
    }

    //--------------------------------------------------------------------------
    // set defaults
    //--------------------------------------------------------------------------

    option->pivot       = SPEX_DEFAULT_PIVOT ;
    option->order       = SPEX_DEFAULT_ORDER ;
    option->print_level = SPEX_DEFAULT_PRINT_LEVEL ;
    option->prec        = SPEX_DEFAULT_PRECISION ;
    option->tol         = SPEX_DEFAULT_TOL ;
    option->round       = SPEX_DEFAULT_MPFR_ROUND ;
    option->check       = false ;

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    return option ;
}

