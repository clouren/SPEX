//------------------------------------------------------------------------------
// SPEX_Utilities/SPEX_vector_free.c: free a vector.
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2020-2021, Jinhao Chen, Chris Lourenco (US Naval Academy),
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// Purpose: This function is called to delete a vector.


#include "spex_util_internal.h"

SPEX_info SPEX_vector_free
(
    SPEX_vector *v_handle,  // vector to be deleted
    const SPEX_options *option
)
{
    if (!spex_initialized ( )) { return (SPEX_PANIC) ; } ;
    if(v_handle == NULL || (*v_handle) == NULL) {return SPEX_OK;}
    spex_delete_mpz_array(&((*v_handle)->x), (*v_handle)->nzmax);
    SPEX_MPQ_CLEAR((*v_handle)->scale);
    SPEX_FREE((*v_handle)->i);
    SPEX_FREE(*v_handle);

    return SPEX_OK;
}

