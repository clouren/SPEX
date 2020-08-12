//------------------------------------------------------------------------------
// SPEX_Util/SPEX_create_default_options: set defaults
//------------------------------------------------------------------------------

// SPEX: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SPEX/License for the license.

//------------------------------------------------------------------------------

/* Purpose: Create and return SPEX_options pointer with default parameters
 * upon successful allocation, which are defined in spex_internal.h
 */

#include "spex_util_internal.h"

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

