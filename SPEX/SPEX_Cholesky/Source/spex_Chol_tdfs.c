//------------------------------------------------------------------------------
// SPEX_Chol/spex_Chol_tdfs: DFS of a tree rooted at a node
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.
//------------------------------------------------------------------------------

#include "spex_chol_internal.h"


/* Purpose: Depth-first search and postorder of a tree rooted at node j */

int64_t spex_Chol_tdfs 
(
    int64_t j,      // Root node
    int64_t k,      // Index (kth node) 
    int64_t* head,  // Head of list
    int64_t* next,  // Next node in the list
    int64_t* post,  // Post ordered tree
    int64_t* stack  // Stack of nodes
)
{ 
    int64_t i, p, top = 0 ;
    if (!head || !next || !post || !stack) return (-1) ;    // check inputs 
    stack [0] = j ;                 // place j on the stack 
    while (top >= 0)                // while (stack is not empty) 
    {
        p = stack [top] ;           // p = top of stack 
        i = head [p] ;              // i = youngest child of p 
        if (i == -1)
        {
            top-- ;                 // p has no unordered children left 
            post [k++] = p ;        // node p is the kth postordered node 
        }
        else
        {
            head [p] = next [i] ;   // remove i from children of p 
            stack [++top] = i ;     // start dfs on child node i 
        }
    }
    return (k) ;
}
