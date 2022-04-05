//------------------------------------------------------------------------------
// SPEX_Util/SPEX_vector_free.c: free a vector.
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2020-2021, Jinhao Chen, Chris Lourenco (US Naval Academy),
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// Purpose: This function is called to delete a vector.


#include "spex_util_internal.h"

SPEX_info SPEX_vector_free
(
    SPEX_vector **v,  // vector to be deleted
    const SPEX_options *option
)
{
    if (!spex_initialized ( )) { return (SPEX_PANIC) ; } ;
    if(v == NULL || (*v) == NULL) {return SPEX_OK;}
    spex_delete_mpz_array(&((*v)->x), (*v)->nzmax);
    SPEX_MPQ_CLEAR((*v)->scale);
    SPEX_FREE((*v)->i);
    SPEX_FREE(*v);

    return SPEX_OK;
}

