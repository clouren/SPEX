//------------------------------------------------------------------------------
// SPEX_CHOLMOD/SPEX_create_mpq_array: create a dense mpq array
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function creates an mpq array of size n, and initialize each of
 * the mpq_t entries.
 */

#include "spex_lu_update_internal.h"

mpq_t* SPEX_create_mpq_array
(
    int64_t n            // size of the array
)
{
    return spex_create_mpq_array(n);
}
