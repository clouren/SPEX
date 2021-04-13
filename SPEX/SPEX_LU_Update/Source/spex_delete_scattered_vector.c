//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_delete_scattered_vector.c: delete a scattered vector.
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is called to delete a scattered mpz vector.


#include "spex_lu_update_internal.h"

void spex_delete_scattered_vector
(
    spex_scattered_vector **sv  // scattered vector to be deleted
)
{
    if(sv == NULL || (*sv) == NULL) {return;}
    SPEX_delete_mpz_array(&((*sv)->x), (*sv)->n);
    SPEX_FREE((*sv)->i);
    SPEX_FREE(*sv);
}

