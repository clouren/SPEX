//------------------------------------------------------------------------------
// SPEX_Chol/spex_Chol_post: Postorder a forest
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#include "spex_chol_internal.h"
#ifndef SPEX_FREE_WORKSPACE
#define SPEX_FREE_WORKSPACE
#endif

/* Purpose: post order a forest. */
SPEX_info spex_Chol_post 
(
    // Output
    int64_t** post_handle, // On input is NULL. On output is post-order of the forest
    // Input
    const int64_t* parent, // Parent[j] is parent of node j in forest
    const int64_t n        // Number of nodes in the forest
)
{
    SPEX_info info ;

    // Input checks
    ASSERT( n >= 0);
    if (!parent) return (SPEX_INCORRECT_INPUT) ;
    
    // Declare variables
    int64_t j, k = 0, *post, *w, *head, *next, *stack ;
    
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
        SPEX_CHECK(spex_Chol_tdfs (&k, j, head, next, post, stack)) ;
    }
    SPEX_free(w);
    (*post_handle) = post;
    return SPEX_OK ; 
}
