//------------------------------------------------------------------------------
// SPEX_print_options: prints out the user specified/default options
//------------------------------------------------------------------------------

// SPEX: (c) 2019-2022, Chris Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function prints out the user specified/default options.
 * this is primarily intended for debugging
 */

#include "demos.h"


void SPEX_print_options     // display specified/default options to user
(
    SPEX_options option     // options (cannot be NULL) FIXME, why??
)
{

    char *piv, *order;
    if (option->order == SPEX_COLAMD || option->order == SPEX_DEFAULT_ORDERING)
    {
        order = "the COLAMD";
    }
    else if (option->order == SPEX_AMD)
    {
        order = "the AMD";
    }
    else if (option->order == SPEX_NO_ORDERING)
    {
        order = "No";
    }
    else
    {
        order = "(undefined)";
    }

    if (option->pivot == SPEX_SMALLEST)
    {
        piv = "smallest";
    }
    else if (option->pivot == SPEX_DIAGONAL)
    {
        piv = "diagonal";
    }
    else if (option->pivot == SPEX_FIRST_NONZERO)
    {
        piv = "first nonzero";
    }
    else if (option->pivot == SPEX_TOL_SMALLEST)
    {
        piv = "diagonal small tolerance";
    }
    else if (option->pivot == SPEX_TOL_LARGEST)
    {
        piv = "diagonal large tolerance";
    }
    else
    {
        piv = "largest";
    }

    printf("\n\n****COMMAND PARAMETERS****");
    printf("\nUsing %s ordering and selecting the %s pivot", order, piv);
    if (option->pivot == SPEX_TOL_SMALLEST ||
        option->pivot == SPEX_TOL_LARGEST)
    {
        printf("\nTolerance used: %lf\n",option->tol);
    }
}
