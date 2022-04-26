//------------------------------------------------------------------------------
// SPEX_Util/SPEX_symbolic_analysis_free: Free memory for the
// SPEX_symbolic_analysis data type.
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2021, Jinhao Chen, Chris Lourenco (US Naval Academy),
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function frees the memory of the SPEX_symbolic_analysis struct
 *
 * Input is the SPEX_symbolic_analysis structure, it is destroyed on function
 * termination.
 */

#include "spex_util_internal.h"

SPEX_info SPEX_symbolic_analysis_free
(
    SPEX_symbolic_analysis **S, // Structure to be deleted
    const SPEX_options *option
)
{
    if (!spex_initialized ( )) return (SPEX_PANIC) ;

    if ((S != NULL) && (*S != NULL))
    {

        SPEX_FREE((*S)->P_perm);
        SPEX_FREE((*S)->Pinv_perm);
        SPEX_FREE((*S)->Q_perm);
        SPEX_FREE((*S)->Qinv_perm);

        SPEX_FREE((*S)->parent);
        SPEX_FREE((*S)->cp);
        SPEX_FREE (*S);
    }

    return (SPEX_OK) ;
}

