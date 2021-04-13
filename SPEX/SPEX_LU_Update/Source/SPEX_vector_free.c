//------------------------------------------------------------------------------
// SPEX_CHOLMOD/SPEX_vector_free.c: free a vector.
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is called to delete a vector.


#include "spex_lu_update_internal.h"

void SPEX_vector_free
(
    SPEX_vector **v  // vector to be deleted
)
{
    if(v == NULL || (*v) == NULL) {return;}
    SPEX_delete_mpz_array(&((*v)->x), (*v)->nzmax);
    SPEX_FREE((*v)->i);
    SPEX_FREE(*v);
}

