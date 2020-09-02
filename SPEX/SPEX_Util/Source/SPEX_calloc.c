//------------------------------------------------------------------------------
// SPEX_Util/SPEX_calloc: wrapper for calloc
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SPEX/License for the license.

//------------------------------------------------------------------------------

// Allocate and initialize memory space for SPEX functions

#include "spex_util_internal.h"

void *SPEX_calloc
(
    size_t nitems,      // number of items to allocate
    size_t size         // size of each item
)
{
    if (!spex_initialized ( )) return (NULL) ;

    return (SuiteSparse_calloc (nitems, size)) ;
}

