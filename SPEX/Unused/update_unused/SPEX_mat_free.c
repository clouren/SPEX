//------------------------------------------------------------------------------
// SPEX_Utilities/SPEX_mat_free.c: free a matrix.
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2020-2021, Jinhao Chen, Chris Lourenco (US Naval Academy),
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// Purpose: This function is called to free a matrix.


#include "spex_util_internal.h"

void SPEX_mat_free
(
    SPEX_mat **A  // matrix to be deleted
)
{
    if(A == NULL || (*A) == NULL) {return;}
    for (int64_t i = 0; i < (*A)->n; i++)
    {
        SPEX_vector_free(&((*A)->v[i]), NULL);
    }
    SPEX_FREE((*A)->v);
    SPEX_MPQ_CLEAR((*A)->scale);
    SPEX_FREE(*A);
}

