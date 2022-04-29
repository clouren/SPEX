//------------------------------------------------------------------------------
// SPEX_Util/SPEX_factorization_free: Free memory for the SPEX_factorization
// data type.
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2021, Jinhao Chen, Chris Lourenco (US Naval Academy),
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function frees the memory of the SPEX_factorization struct
 *
 * Input is the SPEX_factorization structure, it is destroyed on function
 * termination.
 */

#include "spex_util_internal.h"

SPEX_info SPEX_factorization_free
(
    SPEX_factorization **F, // Structure to be deleted
    const SPEX_options *option
)
{
    if (!spex_initialized ( )) return (SPEX_PANIC) ;

    if ((F != NULL) && (*F != NULL))
    {
        SPEX_MPQ_CLEAR((*F)->scale_for_A);

        SPEX_matrix_free(&((*F)->L), option);
        SPEX_matrix_free(&((*F)->U), option);
        SPEX_matrix_free(&((*F)->rhos), option);

        SPEX_FREE((*F)->P_perm);
        SPEX_FREE((*F)->Pinv_perm);
        SPEX_FREE((*F)->Q_perm);
        SPEX_FREE((*F)->Qinv_perm);

        SPEX_FREE (*F) ;
    }

    return (SPEX_OK) ;
}

