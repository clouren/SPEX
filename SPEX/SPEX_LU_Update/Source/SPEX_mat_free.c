//------------------------------------------------------------------------------
// SPEX_CHOLMOD/SPEX_mat_free.c: free a matrix.
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is called to free a matrix.


#include "spex_lu_update_internal.h"

void SPEX_mat_free
(
    SPEX_mat **A  // matrix to be deleted
)
{
    if(A == NULL || (*A) == NULL) {return;}
    for (int64_t i = 0; i < (*A)->n; i++)
    {
        SPEX_vector_free(&((*A)->v[i]));
    }
    SPEX_FREE((*A)->v);
    SPEX_MPQ_CLEAR((*A)->scale);
    SPEX_FREE(*A);
}

