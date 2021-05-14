//------------------------------------------------------------------------------
// SPEX_Util/SPEX_vector_allocate.c: create and initialize a vector with given
// size nzmax.
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2020-2021, Jinhao Chen, Chris Lourenco (US Naval Academy),
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// Purpose: This function is called to create and initialize a mpz vector with
// given size nzmax. The mpz_t vector is allocated with length nzmax. If
// IsSparse is true, then i is allocated with length of nzmax. Otherwise,
// the nnz pattern vector i is set to NULL.

#define SPEX_FREE_ALL \
    SPEX_vector_free (&v, option) ;

#include "spex_util_internal.h"

SPEX_info SPEX_vector_allocate
(
    SPEX_vector **v_handle,         // vector to be allocated
    const int64_t nzmax,            // number of nnz entries in v
    const bool IsSparse,            // indicate if the vector is sparse
    const SPEX_options *option
)
{
    if (!spex_initialized()) return (SPEX_PANIC);
    if (v_handle == NULL || nzmax < 0)
    {
        return SPEX_INCORRECT_INPUT;
    }
    *v_handle = NULL;

    SPEX_info info;
    SPEX_vector *v = (SPEX_vector*) SPEX_malloc(sizeof(SPEX_vector));
    if (!v)
    {
        return SPEX_OUT_OF_MEMORY;
    }

    v->x = NULL;
    v->i = NULL;
    v->nzmax = nzmax;
    v->nz = 0;
    SPEX_MPQ_SET_NULL(v->scale);

    // initialize and set v->scale = 1
    SPEX_CHECK(SPEX_mpq_init(v->scale));
    SPEX_CHECK(SPEX_mpq_set_ui(v->scale, 1, 1));

    // if nzmax == 0
    if (nzmax == 0)
    {
        *v_handle = v;
        return SPEX_OK;
    }

    // allocate for v->x
    v->x = spex_create_mpz_array(nzmax);
    if (!(v->x))
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }

    // allocate for v->i if v is sparse
    if (IsSparse)
    {
        v->i = (int64_t*) SPEX_malloc(nzmax*sizeof(int64_t));
        if (!(v->i))
        {
            SPEX_FREE_ALL;
            return SPEX_OUT_OF_MEMORY;
        }
    }

    *v_handle = v;

    return SPEX_OK;
}

