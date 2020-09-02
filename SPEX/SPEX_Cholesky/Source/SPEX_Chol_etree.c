//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_etree: Compute the elimination tree of a matrix A
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

#include "SPEX_Chol.h"


/* Purpose: Compute the elimination tree of A */

int64_t* SPEX_Chol_etree 
(
    SPEX_matrix* A // Input matrix (must be SPD)
)
{
    // Check input
    ASSERT(A->kind == SPEX_CSC);
    ASSERT(A->type == SPEX_MPZ);
    
    int64_t i, k, p, m, n, inext, *w, *parent, *ancestor, *prev ;
    m = A->m ; n = A->n ;
    parent = (int64_t*) SPEX_malloc(n * sizeof(int64_t));
    w = (int64_t*) SPEX_malloc( (n+m) * sizeof(int64_t));
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
    return parent;
}
