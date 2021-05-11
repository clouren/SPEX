//------------------------------------------------------------------------------
// SPEX_Update/SPEX_vector_alloc.c: create and initialize a vector with given
// size nzmax.
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is called to create and initialize a mpz vector with
// given size nzmax. The mpz_t vector is allocated with length nzmax. If
// IsSparse is true, then i is allocated with length of nzmax. Otherwise,
// the nnz pattern vector i is set to NULL.

#include "spex_update_internal.h"

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

    // initialize and set v->scale = 1
    SPEX_MPQ_SET_NULL(v->scale);
    SPEX_info info = SPEX_mpq_init(v->scale);
    if (info != SPEX_OK)
    {
        SPEX_MPQ_CLEAR(v->scale);
        SPEX_FREE(v);
        return info;
    }

    info = SPEX_mpq_set_ui(v->scale, 1, 1);
    if (info != SPEX_OK)
    {
        SPEX_MPQ_CLEAR(v->scale);
        SPEX_FREE(v);
        return info;
    }

    // if nzmax == 0
    if (nzmax == 0)
    {
        *v_handle = v;
        return SPEX_OK;
    }

    // allocate for v->x
    v->x = SPEX_create_mpz_array(nzmax);
    if (!(v->x))
    {
        SPEX_MPQ_CLEAR(v->scale);
        SPEX_FREE(v);
        return SPEX_OUT_OF_MEMORY;
    }

    // allocate for v->i if v is sparse
    if (IsSparse)
    {
        v->i = (int64_t*) SPEX_malloc(nzmax*sizeof(int64_t));
        if (!(v->i))
        {
            SPEX_delete_mpz_array(&(v->x), nzmax);
            SPEX_MPQ_CLEAR(v->scale);
            SPEX_FREE(v);
            return SPEX_OUT_OF_MEMORY;
        }
    }

    *v_handle = v;

    return SPEX_OK;
}

