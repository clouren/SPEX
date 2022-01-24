//------------------------------------------------------------------------------
// SPEX_Chol/spex_Chol_etree: Compute the elimination tree of a matrix A
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#include "spex_chol_internal.h"


/* Purpose: Compute the elimination tree of A */

SPEX_info spex_chol_etree 
(
    // Output
    int64_t** tree_handle,      // On output: contains the elimination tree of A
                                // On input: undefined.
    // Input
    const SPEX_matrix* A        // Input matrix (must be SPD). Note to compute
)
{
    // Check input
    ASSERT(A->kind == SPEX_CSC);
    ASSERT(A->type == SPEX_MPZ);
    ASSERT(A->m >= 0);
    ASSERT(A->n >= 0);

    // Declare variables
    int64_t i, k, p, m, n, inext, *w, *parent, *ancestor, *prev ;
    m = A->m ; n = A->n ;
    
    // Allocate parent
    parent = (int64_t*) SPEX_malloc( n * sizeof(int64_t));
    // Allocate workspace
    w = (int64_t*) SPEX_malloc( (n+m) * sizeof(int64_t) );
    if (!parent || !w)
        return SPEX_OUT_OF_MEMORY;
    ancestor = w ; prev = w + n ;
    for (k = 0 ; k < n ; k++)
    {
        parent [k] = -1 ;                           // node k has no parent yet 
        ancestor [k] = -1 ;                         // nor does k have an ancestor 
        for (p = A->p [k] ; p < A->p [k+1] ; p++)
        {
            i = A->i [p] ;
            for ( ; i != -1 && i < k ; i = inext)   // traverse from i to k 
            {
                inext = ancestor [i] ;              // int64_t*ext = ancestor of i 
                ancestor [i] = k ;                  // path compression 
                if (inext == -1) parent [i] = k ;   // no anc., parent is k 
            }
        }
    }
    SPEX_FREE(w);
    (*tree_handle) = parent;
    return SPEX_OK;
}
