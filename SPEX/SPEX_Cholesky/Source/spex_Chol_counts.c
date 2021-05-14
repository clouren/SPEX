//------------------------------------------------------------------------------
// SPEX_Chol/spex_Chol_Counts: Column counts for Cholesky factorization
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

#include "spex_chol_internal.h"

#define HEAD(k,j) (j )
#define NEXT(J)   (-1)

/* Purpose: Obtain the column counts of an SPD matrix for Cholesky factorization*
 * This is a modified version of Csparse's cs_chol_counts function
 */
SPEX_info spex_Chol_counts 
(
    int64_t** c_handle,
    const SPEX_matrix *A, 
    int64_t *parent, 
    int64_t *post
)
{
    SPEX_info info;
    int64_t i, j, k, n, m, J, s, p, q, jleaf, *maxfirst, *prevleaf,
        *ancestor, *head = NULL, *next = NULL, *colcount, *w, *first, *delta ;
    if (!A || !parent || !post) return (SPEX_INCORRECT_INPUT) ;    /* check inputs */
    m = A->m ; n = A->n ;
    // Can not have negative m or n
    ASSERT(n >= 0);
    ASSERT(m >= 0);
    // Size of workspace
    s = 4*n ;
    // Allocate result in delta
    delta = colcount = (int64_t*) SPEX_malloc (n* sizeof (int64_t)) ;
    // Create a workspace of size s
    w = (int64_t*) SPEX_malloc (s* sizeof (int64_t)) ;
    ancestor = w ; maxfirst = w+n ; prevleaf = w+2*n ; first = w+3*n ;
    // Clear workspace
    for (k = 0 ; k < s ; k++)
    {
        w [k] = -1 ;
    }
    // Find first j
    for (k = 0 ; k < n ; k++)
    {
        j = post [k] ;
        delta [j] = (first [j] == -1) ? 1 : 0 ;  /* delta[j]=1 if j is a leaf */
        for ( ; j != -1 && first [j] == -1 ; j = parent [j]) 
        {
            first [j] = k ;
        }
    }
    // Initialize ancestor of each node
    for (i = 0 ; i < n ; i++)
    {
        ancestor [i] = i ;
    }
    for (k = 0 ; k < n ; k++)
    {
        j = post [k] ;          /* j is the kth node in postordered etree */
        if (parent [j] != -1) 
        {
            delta [parent [j]]-- ;    /* j is not a root */
        }
        for (J = HEAD (k,j) ; J != -1 ; J = NEXT (J))   /* J=j for LL'=A case */
        {
            for (p = A->p [J] ; p < A->p [J+1] ; p++)
            {
                i = A->i [p] ;
                SPEX_CHECK(spex_Chol_leaf (&q, i, j, first, maxfirst, prevleaf, ancestor, &jleaf));
                if (jleaf >= 1)
                {
                    delta [j]++ ;   /* A(i,j) is in skeleton */
                }
                if (jleaf == 2)
                {
                    delta [q]-- ;   /* account for overlap in q */
                }
            }
        }
        if (parent [j] != -1) 
        {
            ancestor [j] = parent [j] ;
        }
    }
    for (j = 0 ; j < n ; j++)           /* sum up delta's of each child */
    {
        if (parent [j] != -1)
        {
            colcount [parent [j]] += colcount [j] ;
        }
    }
    SPEX_FREE(w);
    (*c_handle) = colcount;
    return SPEX_OK;    
} 
