//------------------------------------------------------------------------------
// SPEX_Util/SPEX_vector_realloc: realloc the space for a SPEX_vector object
// to given new size.
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2020-2021, Jinhao Chen, Chris Lourenco (US Naval Academy),
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function realloc a SPEX_vector to given size. It will
 * initialize/allocate for the mpz entries.
 */

#include "spex_util_internal.h"

SPEX_info SPEX_vector_realloc
(
    SPEX_vector* v,         // the vector to be expanded
    const int64_t new_size, // desired new size for v
    const SPEX_options *option
)
{
    SPEX_info info;
    if (!spex_initialized ( )) { return (SPEX_PANIC) ; } ;

    int64_t p, old_size = v->nzmax ;
    if (old_size == new_size) {return SPEX_OK;}

    //--------------------------------------------------------------------------
    // free mpz entries before shrinking the mpz vector
    //--------------------------------------------------------------------------

    if (old_size > new_size)
    {
        for (p = new_size; p < old_size; p++)
        {
            if (v->x[p] != NULL)
            {
                SPEX_MPZ_CLEAR(v->x[p]);
            }
        }
    }

    //--------------------------------------------------------------------------
    // expand the size of v->x and v->i to new_size
    //--------------------------------------------------------------------------

    bool okx, oki ;
    v->x = (mpz_t *)
        SPEX_realloc (new_size, old_size, sizeof (mpz_t), v->x, &okx) ;
    v->i = (int64_t *)
        SPEX_realloc (new_size, old_size, sizeof (int64_t), v->i, &oki) ;
    if (!oki || !okx)
    {
        return (SPEX_OUT_OF_MEMORY) ;
    }

    v->nzmax = new_size ;

    //--------------------------------------------------------------------------
    // set newly allocated mpz entries to NULL and initialize if required
    //--------------------------------------------------------------------------

    if (old_size < new_size)
    {
        for (p = old_size ; p < new_size ; p++)
        {
            SPEX_MPZ_SET_NULL (v->x[p]) ;
        }

        for (p = old_size ; p < new_size ; p++)
        {
            SPEX_CHECK(SPEX_mpz_init (v->x[p])) ;
        }
    }

    return (SPEX_OK) ;
}
