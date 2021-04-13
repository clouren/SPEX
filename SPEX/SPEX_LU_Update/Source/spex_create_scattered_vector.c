//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_create_scattered_vector.c: create and initialize a
// scattered vector with given size n.
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is called to create and initialize a scattered mpz
// vector with given size n. Both mpz_t vector and the nnz pattern vector are
// allocated with length n.


#include "spex_lu_update_internal.h"

spex_scattered_vector* spex_create_scattered_vector
(
    const int64_t n             // number of entries in v
)
{
    if (n <= 0) { return NULL; }

    spex_scattered_vector *sv = (spex_scattered_vector*)
                                SPEX_malloc(sizeof(spex_scattered_vector));
    if (!sv)
    {
        return NULL;
    }

    sv->x = SPEX_create_mpz_array(n);
    if (!(sv->x))
    {
        SPEX_FREE(sv);
        return NULL;
    }
    sv->i = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    if (!(sv->i))
    {
        SPEX_delete_mpz_array(&(sv->x), n);
        SPEX_FREE(sv);
        return NULL;
    }

    sv->n = n;
    sv->nz = 0;

    return sv;
}

