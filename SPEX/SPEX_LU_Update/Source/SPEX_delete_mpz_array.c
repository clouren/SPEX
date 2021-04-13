//------------------------------------------------------------------------------
// SPEX_CHOLMOD/SPEX_delete_mpz_array.c: clear the memory used for a mpz array
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function clears the memory used for an mpz vector of size n. 
 * Call this function for all mpz vectors when done.
 * 
 * Input is a mpz_t** array which is destroyed upon function completion
 * 
 */

#include "spex_lu_update_internal.h"

void SPEX_delete_mpz_array
(
    mpz_t **x,      // mpz array to be deleted
    int64_t n       // Size of x
)
{
    if (x == NULL || (*x) == NULL) {return;}
    for (int64_t i = 0; i < n; i++)
    {
        if ( (*x) [i] != NULL)
        {
            SPEX_MPZ_CLEAR((*x) [i]);
        }
    }
    SPEX_FREE ((*x));
}
