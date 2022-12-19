//------------------------------------------------------------------------------
// SPEX_determine_error: auxiliary file for test coverage (tcov)
//------------------------------------------------------------------------------

// SPEX: (c) 2019-2022, Chris Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: Determine why a SPEX function failed
 */

#include "SPEX.h"


void SPEX_determine_error
(
    SPEX_info ok
)
{
    if (ok == SPEX_OUT_OF_MEMORY)
    {
        printf("\nOut of memory\n");
    }
    else if (ok == SPEX_SINGULAR)
    {
        printf("\nInput matrix is singular OR no diagonal pivot. Please ensure input is SPD\n");
    }
    else if (ok == SPEX_UNSYMMETRIC)
    {
        printf("\nInput matrix is unsymmetric, please try left LU\n");
    }
    else if (ok == SPEX_INCORRECT_INPUT)
    {
        printf("\nIncorrect input for a SPEX_Chol Function\n");
    }
}
