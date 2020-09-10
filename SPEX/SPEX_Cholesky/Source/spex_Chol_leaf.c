//------------------------------------------------------------------------------
// SPEX_Chol/spex_Chol_leaf: Subroutine for column counts of Cholesky
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

#include "spex_chol_internal.h"

/* Purpose: consider A(i,j), node j in ith row subtree and return lca(jprev,j) 
   Used to determine Column counts of Cholesky factor*/
int64_t spex_Chol_leaf 
(
    int64_t i, 
    int64_t j, 
    int64_t* first, 
    int64_t* maxfirst, 
    int64_t* prevleaf,
    int64_t* ancestor, 
    int64_t* jleaf
)
{
    int64_t q, s, sparent, jprev ;
    // Check inputs
    if (!first || !maxfirst || !prevleaf || !ancestor || !jleaf) return (-1) ;
    
    *jleaf = 0 ;
    if (i <= j || first [j] <= maxfirst [i]) return (-1) ;  // j not a leaf 
    maxfirst [i] = first [j] ;      // update max first[j] seen so far 
    jprev = prevleaf [i] ;          // jprev = previous leaf of ith subtree 
    prevleaf [i] = j ;
    *jleaf = (jprev == -1) ? 1: 2 ; // j is first or subsequent leaf 
    if (*jleaf == 1) return (i) ;   // if 1st leaf, q = root of ith subtree 
    for (q = jprev ; q != ancestor [q] ; q = ancestor [q]) ;
    for (s = jprev ; s != q ; s = sparent)
    {
        sparent = ancestor [s] ;    // path compression 
        ancestor [s] = q ;
    }
    return (q) ;                    // q = least common ancester (jprev,j) 
}
