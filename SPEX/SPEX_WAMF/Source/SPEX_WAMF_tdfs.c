//------------------------------------------------------------------------------
// SPEX_WAMF/SPEX_WAMF_tdfs: Depth first search of a tree
//------------------------------------------------------------------------------

// SPEX_WAMF: (c) 2021, Chris Lourenco US Naval Academy, Erick Moreno-Centeno, Texas
// A&M University. All Rights Reserved.  See SPEX_WAMF/License for the license.

//------------------------------------------------------------------------------


#include "SPEX_WAMF.h"

/* Purpose: depth-first search and postorder of a tree rooted at node j */
int64_t SPEX_WAMF_tdfs 
(
    int64_t j, 
    int64_t k, 
    int64_t *head, 
    int64_t *next, 
    int64_t *post, 
    int64_t *stack
)
{
    int64_t i, p, top = 0 ;
    if (!head || !next || !post || !stack) return (-1) ;    /* check inputs */
    stack [0] = j ;                 /* place j on the stack */
    while (top >= 0)                /* while (stack is not empty) */
    {
        p = stack [top] ;           /* p = top of stack */
        i = head [p] ;              /* i = youngest child of p */
        if (i == -1)
        {
            top-- ;                 /* p has no unordered children left */
            post [k++] = p ;        /* node p is the kth postordered node */
        }
        else
        {
            head [p] = next [i] ;   /* remove i from children of p */
            stack [++top] = i ;     /* start dfs on child node i */
        }
    }
    return (k) ;
}
