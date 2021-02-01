//------------------------------------------------------------------------------
// SPEX_WAMF/SPEX_WAMF_wclear: Clear workspace
//------------------------------------------------------------------------------

// SPEX_WAMF: (c) 2021, Chris Lourenco US Naval Academy, Erick Moreno-Centeno, Texas
// A&M University. All Rights Reserved.  See SPEX_WAMF/License for the license.

//------------------------------------------------------------------------------

#include "SPEX_WAMF.h"

/* Purpose: drop diagonal entries in a matrix when used as a function */

int64_t SPEX_WAMF_wclear
(
    int64_t mark, 
    int64_t lemax, 
    int64_t *w, 
    int64_t n
)
{
    int64_t k ;
    if (mark < 2 || (mark + lemax < 0))
    {
        for (k = 0 ; k < n ; k++) if (w [k] != 0) w [k] = 1 ;
        mark = 2 ;
    }
    return (mark) ;     /* at this point, w [0..n-1] < mark holds */
}
