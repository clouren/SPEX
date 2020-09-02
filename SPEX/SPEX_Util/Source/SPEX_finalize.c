//------------------------------------------------------------------------------
// SPEX_Util/SPEX_finalize: finalize SPEX_LU
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SPEX/License for the license.

//------------------------------------------------------------------------------

// SPEX_finalize frees the working environment for SPEX library.

#include "spex_util_internal.h"

SPEX_info SPEX_finalize
(
    void
)
{
    if (!spex_initialized ( )) return (SPEX_PANIC) ;

    SPEX_mpfr_free_cache ( ) ;    // Free mpfr internal cache
    spex_gmp_finalize ( ) ;       // Reset GMP memory variables

    spex_set_initialized (false) ;
    return (SPEX_OK) ;
}

