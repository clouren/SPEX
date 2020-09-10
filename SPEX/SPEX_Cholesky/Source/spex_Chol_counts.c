//------------------------------------------------------------------------------
// SPEX_Chol/spex_Chol_Counts: Column counts for Cholesky factorization
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

#include "spex_chol_internal.h"

//#define HEAD(k,j) (ata ? head [k] : j)
//#define NEXT(J)   (ata ? next [J] : -1)
#define HEAD(k,j) (j )
#define NEXT(J)   (-1)


/* Purpose: Obtaint the column counts of an SPD matrix for Cholesky factorization */
int64_t* spex_Chol_counts 
(
    SPEX_matrix *A, 
    int64_t *parent, 
    int64_t *post
)
{
    int64_t i, j, k, n, m, J, s, p, q, jleaf, *maxfirst, *prevleaf,
        *ancestor, *head = NULL, *next = NULL, *colcount, *w, *first, *delta ;
    if (!A || !parent || !post) return (NULL) ;    /* check inputs */
    m = A->m ; n = A->n ;
    s = 4*n ;
    delta = colcount = (int64_t*) SPEX_malloc (n* sizeof (int64_t)) ;    /* allocate result */
    w = (int64_t*) SPEX_malloc (s* sizeof (int64_t)) ;                   /* get workspace */
    ancestor = w ; maxfirst = w+n ; prevleaf = w+2*n ; first = w+3*n ;
    for (k = 0 ; k < s ; k++) w [k] = -1 ;      /* clear workspace w [0..s-1] */
    for (k = 0 ; k < n ; k++)                   /* find first [j] */
    {
        j = post [k] ;
        delta [j] = (first [j] == -1) ? 1 : 0 ;  /* delta[j]=1 if j is a leaf */
        for ( ; j != -1 && first [j] == -1 ; j = parent [j]) first [j] = k ;
    }
    for (i = 0 ; i < n ; i++) ancestor [i] = i ; /* each node int64_t* its own set */
    for (k = 0 ; k < n ; k++)
    {
        j = post [k] ;          /* j is the kth node in postordered etree */
        if (parent [j] != -1) delta [parent [j]]-- ;    /* j is not a root */
        for (J = HEAD (k,j) ; J != -1 ; J = NEXT (J))   /* J=j for LL'=A case */
        {
            for (p = A->p [J] ; p < A->p [J+1] ; p++)
            {
                i = A->i [p] ;
                q = spex_Chol_leaf (i, j, first, maxfirst, prevleaf, ancestor, &jleaf);
                if (jleaf >= 1) delta [j]++ ;   /* A(i,j) is in skeleton */
                if (jleaf == 2) delta [q]-- ;   /* account for overlap in q */
            }
        }
        if (parent [j] != -1) ancestor [j] = parent [j] ;
    }
    for (j = 0 ; j < n ; j++)           /* sum up delta's of each child */
    {
        if (parent [j] != -1) colcount [parent [j]] += colcount [j] ;
    }
    SPEX_FREE(w);
    return colcount;    
} 
