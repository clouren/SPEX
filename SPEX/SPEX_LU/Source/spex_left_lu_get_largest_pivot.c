//------------------------------------------------------------------------------
// SPEX_LU/spex_left_lu_get_largest_pivot: find a pivot entry in a column
//------------------------------------------------------------------------------

// SPEX_LU: (c) 2019-2025, Christopher Lourenco, Jinhao Chen,,
// Erick Moreno-Centeno, and Timothy A. Davis. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function selects the pivot element as the largest in the
 * column This is activated if the user sets option->pivot = SPEX_LARGEST.
 *
 * Note: This pivoting scheme is NOT recommended for SPEX Left LU.  It is provided
 * for comparison with other pivoting options.
 *
 * On output, the index of the largest pivot is returned.
 */

#define SPEX_FREE_ALL          \
    SPEX_mpz_clear (big);

#include "spex_lu_internal.h"

SPEX_info spex_left_lu_get_largest_pivot
(
    int64_t *pivot,         // the row index of largest pivot
    int64_t *p_pivot,       // pivot is located in xi [*p_pivot]
    SPEX_matrix x,          // kth column of L and U
    int64_t *pivs,          // vector which indicates whether each row
                            // has been pivotal
    int64_t n,              // dimension of problem
    int64_t top,            // nonzero pattern is located in xi[top..n-1]
    int64_t *xi             // nonzero pattern of x
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SPEX_REQUIRE(x, SPEX_DENSE, SPEX_MPZ);

    SPEX_info info ;
    if (!pivs || !xi || !pivot || !p_pivot) {return SPEX_INCORRECT_INPUT;}

    //--------------------------------------------------------------------------
    // allocate workspace
    //--------------------------------------------------------------------------

    int64_t p, i ;
    int r ;
    (*pivot) = -1 ;
    (*p_pivot) = -1 ;
    mpz_t big ;
    SPEX_mpz_set_null (big);
    SPEX_MPZ_INIT (big);

    //--------------------------------------------------------------------------
    // Iterate accross the nonzeros in x
    //--------------------------------------------------------------------------

    for (p = top; p < n; p++)
    {
        // Location of the ith nonzero
        i = xi[p];
        // i can be pivotal
        SPEX_MPZ_CMPABS(&r, big, x->x.mpz[i]);
        if (pivs[i] < 0 && r < 0)
        {
            // Current largest pivot location
            (*pivot) = i;
            (*p_pivot) = p ;
            // Current largest pivot value
            SPEX_MPZ_SET(big, x->x.mpz[i]);
        }
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    SPEX_FREE_ALL;
    if ((*pivot) == -1)
    {
        return SPEX_SINGULAR;
    }
    else
    {
        return SPEX_OK;
    }
}

