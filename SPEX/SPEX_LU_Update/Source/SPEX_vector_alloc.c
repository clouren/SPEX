//------------------------------------------------------------------------------
// SPEX_CHOLMOD/SPEX_vector_alloc.c: create and initialize a vector with given
// size nzmax.
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is called to create and initialize a mpz vector with
// given size nzmax. The mpz_t vector is allocated with length nzmax. If
// IsSparse is true, then i is allocated with length of nzmax. Otherwise,
// the nnz pattern vector i is set to NULL.

#include "spex_lu_update_internal.h"

SPEX_info SPEX_vector_alloc
(
    SPEX_vector **v_handle,         // vector to be allocated
    const int64_t nzmax,            // number of nnz entries in v
    const bool IsSparse             // indicate if the vector is sparse
)
{
    if (!spex_initialized()) return (SPEX_PANIC);
    if (v_handle == NULL || nzmax < 0)
    {
        return SPEX_INCORRECT_INPUT;
    }
    *v_handle = NULL;

    SPEX_vector *v = (SPEX_vector*) SPEX_malloc(sizeof(SPEX_vector));
    if (!v)
    {
        return SPEX_OUT_OF_MEMORY;
    }

    v->x = NULL;
    v->i = NULL;
    v->nzmax = nzmax;
    v->nz = 0;
    if (nzmax == 0)
    {
        *v_handle = v;
        return SPEX_OK;
    }

    v->x = SPEX_create_mpz_array(nzmax);
    if (!(v->x))
    {
        SPEX_FREE(v);
        return SPEX_OUT_OF_MEMORY;
    }
    if (IsSparse)
    {
        v->i = (int64_t*) SPEX_malloc(nzmax*sizeof(int64_t));
        if (!(v->i))
        {
            SPEX_delete_mpz_array(&(v->x), nzmax);
            SPEX_FREE(v);
            return SPEX_OUT_OF_MEMORY;
        }
    }

    *v_handle = v;

    return SPEX_OK;
}

