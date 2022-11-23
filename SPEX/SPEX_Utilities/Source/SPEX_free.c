//------------------------------------------------------------------------------
// SPEX_Utilities/SPEX_free: wrapper for free
//------------------------------------------------------------------------------

//TEST 

// SPEX_Utilities: (c) 2019-2022, Chris Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// Free the memory allocated by SPEX_calloc, SPEX_malloc, or SPEX_realloc.

#include "spex_util_internal.h"

void SPEX_free
(
    void *p         // pointer to memory space to free
)
{
    SuiteSparse_free (p) ;
}

