//------------------------------------------------------------------------------
// SPEX_Util/SPEX_free: wrapper for free
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SPEX/License for the license.

//------------------------------------------------------------------------------

// Free the memory allocated by SPEX_calloc, SPEX_malloc, or SPEX_realloc.

#include "spex_util_internal.h"

void SPEX_free
(
    void *p         // pointer to memory space to free
)
{
    if (!spex_initialized ( )) return ;
    SuiteSparse_free (p) ;
}

