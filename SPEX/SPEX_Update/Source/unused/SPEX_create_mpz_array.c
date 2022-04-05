//------------------------------------------------------------------------------
// SPEX_Update/SPEX_create_mpz_array: create a dense mpz array
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function creates an mpz array of size n, and initialize each of * the mpz_t entries.
 */

#include "spex_update_internal.h"

mpz_t* SPEX_create_mpz_array
(
    int64_t n            // size of the array
)
{
    return spex_create_mpz_array(n);
}
