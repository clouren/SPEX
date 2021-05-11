//------------------------------------------------------------------------------
// SPEX_Update/SPEX_mat_alloc.c: create and initialize a matrix with given
// size n and m.
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is called to create and initialize a matrix with
// given number of vectors n, and number of max entries in each vector m.  If
// IsSparse is true, nzmax in each vector is initialized as 0, and thus further
// reallocation using SPEX_vector_realloc is needed before assigning value to
// any entries. Otherwise (IsSparse is false), each vector will be initialized
// as a dense vector, which means A->v->x will be a mpz_t vector of size m and
// A->v->i will be NULL.  scale will be set to be 1 as default.

#define SPEX_FREE_ALL       \
    SPEX_mat_free(&A);

#include "spex_update_internal.h"

SPEX_info SPEX_mat_alloc
(
    SPEX_mat **A_handle,      // matrix to be allocated
    const int64_t n,             // number of vectors
    const int64_t m,             // number of max entries in each vector
    const bool IsSparse          // Indicate if the matrix is sparse
)
{
    SPEX_info info ;
    if (!spex_initialized ( )) return (SPEX_PANIC) ;
    if (A_handle == NULL || n < 0 || m < 0)
    {
        return SPEX_INCORRECT_INPUT;
    }
    *A_handle = NULL;

    if (n == 0) { return SPEX_OK; }

    SPEX_mat *A = (SPEX_mat*) SPEX_malloc(sizeof(SPEX_mat));
    if (!A)
    {
        return SPEX_OUT_OF_MEMORY;
    }
    A->m = m;
    A->n = n;
    SPEX_MPQ_SET_NULL(A->scale);

    // make sure each A->v[] is initialized as NULL
    A->v = (SPEX_vector**) SPEX_calloc(n, sizeof(SPEX_vector*));
    if (!(A->v))
    {
        SPEX_FREE(A);
        return SPEX_OUT_OF_MEMORY;
    }
    
    for (int64_t i = 0; i < n; i++)
    {
        SPEX_CHECK(SPEX_vector_alloc(&(A->v[i]), (IsSparse ? 0 : m), IsSparse));
    }
    SPEX_CHECK(SPEX_mpq_init(A->scale));
    SPEX_CHECK(SPEX_mpq_set_ui(A->scale, 1, 1));

    *A_handle = A;

    return SPEX_OK;
}

