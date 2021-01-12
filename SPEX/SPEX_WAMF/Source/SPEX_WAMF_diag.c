//------------------------------------------------------------------------------
// SPEX_WAMF/SPEX_WAMF_diag: Drop diagonal entries when used as a function
//------------------------------------------------------------------------------

// SPEX_WAMF: (c) 2021, Chris Lourenco US Naval Academy, Erick Moreno-Centeno, Texas
// A&M University. All Rights Reserved.  See SPEX_WAMF/License for the license.

//------------------------------------------------------------------------------

#include "SPEX_WAMF.h"

/* Purpose: drop diagonal entries in a matrix when used as a function */

int64_t SPEX_WAMF_diag
(
    int64_t i, 
    int64_t j, 
    double aij, 
    void *other
) 
{ 
    return (i != j) ; 
}
