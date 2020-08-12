//------------------------------------------------------------------------------
// SPEX_Util/SPEX_malloc: wrapper for malloc
//------------------------------------------------------------------------------

// SPEX: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SPEX_LU/License for the license.

//------------------------------------------------------------------------------

// Allocate memory space for SPEX functions.

#include "spex_util_internal.h"

void *SPEX_malloc
(
    size_t size        // size of memory space to allocate
)
{
    if (!spex_initialized ( )) return (NULL) ;
    return (SuiteSparse_malloc (1, size)) ;
}

