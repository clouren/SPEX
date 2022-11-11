//------------------------------------------------------------------------------
// SPEX_Update/spex_update_forward_sub.c: sparse forward substitution,
// i.e., compute x = (LD)\v
//------------------------------------------------------------------------------
    
// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------
    
// Purpose: This function is to perform sparse forward substitution, which is
// essentially the same as sparse REF triangular solve for LDx=v, but with v
// as a dense vector. This function assumes v is in the same row permutation as
// L. This function take v as input using x_input, and the solution is stored
// in x_output. In case when v has multiple column, simply iterate this function
// and solve each column at each iteration.

#include "spex_update_internal.h"

#define SL(k) (L->v[(k)]->scale)

SPEX_info spex_update_forward_sub // perform sparse forward substitution
(
    SPEX_vector *x,     // Input: the right-hand-side vector
                        // Output: solution x
    const SPEX_matrix L,  // matrix L
    const int64_t *P,   // row permutation
    const SPEX_matrix rhos,// array of scaled pivots
    int64_t *h          // history vector for x
)
{
    SPEX_info info;
    int sgn;
    int64_t i, n = L->n;

    // reset history vector
    for (i = 0; i < n; i++)
    {
        h[i] = -1; 
    }

    for (i = 0; i < n; i++)
    {
        // skip if x(P[i]) == 0
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, x->x[P[i]]));
        if (sgn == 0)       { continue; }

        // perform i-th IPGE update for x
        SPEX_CHECK(spex_update_ipge(x, h, NULL, L->v[i], P, NULL, rhos, i));
    }

    return SPEX_OK; 
}
