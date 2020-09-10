//------------------------------------------------------------------------------
// SPEX_Chol/spex_Chol_post: Postorder a forest
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

#include "spex_chol_internal.h"

/* Purpose: post order a forest. Modified from csparse */
SPEX_info spex_Chol_post 
(
    int64_t** post_handle,
    int64_t* parent,    // Parent[j] is parent of node j int64_t* forest
    int64_t n           // Number of nodes int64_t* the forest
)
{
    int64_t j, k = 0, *post, *w, *head, *next, *stack ;
    if (!parent) return (SPEX_INCORRECT_INPUT) ;                                // check inputs 
    // Allocate the postordering result
    post = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    // Create a workspace
    w = (int64_t*) SPEX_malloc (3*n* sizeof (int64_t)) ;
    if (!w || !post) return (SPEX_OUT_OF_MEMORY) ;
    head = w ; next = w + n ; stack = w + 2*n ;
    // Empty linked lists
    for (j = 0 ; j < n ; j++)
    {
        head [j] = -1 ;
    }
    for (j = n-1 ; j >= 0 ; j--)            // traverse nodes in reverse order
    {
        if (parent [j] == -1) continue ;    // j is a root 
        next [j] = head [parent [j]] ;      // add j to list of its parent 
        head [parent [j]] = j ;
    }
    for (j = 0 ; j < n ; j++)
    {
        if (parent [j] != -1) continue ;    // skip j if it is not a root 
        k = spex_Chol_tdfs (j, k, head, next, post, stack) ;
    }
    SPEX_free(w);
    (*post_handle) = post;
    return SPEX_OK ; 
}
