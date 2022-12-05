//------------------------------------------------------------------------------
// SPEX_Utilities/SPEX_finalize: finalize SPEX
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2019-2022, Chris Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// SPEX_finalize frees the working environment for SPEX library.

#include "spex_util_internal.h"

SPEX_info SPEX_finalize
(
    void
)
{

    if (!spex_initialized ( )) { return (SPEX_PANIC) ; }

    SPEX_mpfr_free_cache ( ) ;    // Free mpfr internal cache
    SPEX_gmp_finalize ( ) ;       // Reset GMP memory variables

    spex_set_initialized (false) ;
    return (SPEX_OK) ;
}

