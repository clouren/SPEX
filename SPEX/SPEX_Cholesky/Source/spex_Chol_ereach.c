//------------------------------------------------------------------------------
// SPEX_Chol/spex_Chol_ereach: Compute reach of an elimination tree
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

#include "spex_chol_internal.h"

/* Purpose: This function computes the reach of the kth row of A onto the graph of L using the 
   elimination tree. It finds the nonzero pattern of row k of L and uses the upper triangular 
   part of A(:,k) */
   
SPEX_info spex_Chol_ereach 
(
    // Output
    int64_t* top_handle,
    int64_t* xi,            // Contains the nonzero pattern in xi[top..n-1]
    // Input
    const SPEX_matrix* A,   // Matrix to be analyzed
    int64_t k,              // Node to start at
    int64_t* parent,        // ELimination Tree
    int64_t* w               // Workspace array
)
{
    int64_t i, p, n, len, top ;
    // Check inputs
    SPEX_REQUIRE(A, SPEX_CSC, SPEX_MPZ);
    if (!parent || !xi || !w) return (SPEX_INCORRECT_INPUT) ;
    ASSERT(A->n >= 0);
    top = n = A->n ; 
    // Mark node k as visited
    SPEX_MARK (w, k) ;
    
    for (p = A->p [k] ; p < A->p [k+1] ; p++)
    {
        // A(i,k) is nonzero
        i = A->i [p] ;
        if (i > k) 
        {
            continue ;  // only use upper triangular part of A 
        }
        for (len = 0 ; !SPEX_MARKED (w,i) ; i = parent [i]) // traverse up etree
        {
            xi [len++] = i ;           // L(k,i) is nonzero 
            SPEX_MARK (w, i) ;        // mark i as visited 
        }
        while (len > 0) xi [--top] = xi [--len] ; // push path onto stack 
    }
    for (p = top ; p < n ; p++) SPEX_MARK (w, xi [p]) ;    // unmark all nodes 
    SPEX_MARK (w, k) ;                // unmark node k 
    (*top_handle) = top;
    return SPEX_OK ;                    // xi [top..n-1] contains pattern of L(k,:)
}
