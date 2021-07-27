//------------------------------------------------------------------------------
// SPEX_Backslash/spex_backslash_set_defaults.c: Set defaults for factorizing matrix
//------------------------------------------------------------------------------

// SPEX_Backslash: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* Purpose: Set default orderings for Cholesky vs left LU
 * 
 * Input/Output arguments:
 *
 * option:      Command options
 *
 * lu:          If true, LU factorization is used, if false, Cholesky is used
 */

#include "SPEX_Backslash.h"

SPEX_info spex_backslash_set_defaults
(
    SPEX_options* option,
    bool lu
)
{
    ASSERT( option != NULL);
    if (lu)
    {
        option->order = SPEX_COLAMD;
        option->pivot = SPEX_TOL_SMALLEST;
    }
    else
    {
        option->order = SPEX_AMD;
        option->pivot = SPEX_DIAGONAL;
    }
}
