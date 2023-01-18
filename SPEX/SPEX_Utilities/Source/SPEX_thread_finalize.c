//------------------------------------------------------------------------------
// SPEX_Utilities/SPEX_thread_finalize: finish SPEX for a single user thread
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2019-2022, Chris Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// SPEX_thread_finalize frees the working evironment for SPEX for a
// single user thread.

// FIXME: add to user guide

#include "spex_util_internal.h"

SPEX_info SPEX_thread_finalize ( void )
{
    if (!spex_initialized ( )) return (SPEX_PANIC) ;
    spex_gmp_finalize ( ) ;
    return (SPEX_OK) ;
}

