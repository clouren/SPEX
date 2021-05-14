//------------------------------------------------------------------------------
// SPEX_Chol/spex_chol_internal: include file for internal use in SPEX_Cholesky
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

// This file is not intended to be #include'd in user applications.  Use
// SPEX_Chol.h instead.

#ifndef SPEX_CHOL_INTERNAL_H
#define SPEX_CHOL_INTERNAL_H


#include "spex_util_internal.h"
#include "SPEX_Chol.h"

// ============================================================================
//                           Internal Functions
// ============================================================================

//TODO "separate" inputs and outputs in every function. specify "null" if so in comments
//TODO group similar functions together, check order

/* Purpose: Compute the elimination tree of A */
SPEX_info spex_Chol_etree 
(
    int64_t** tree,
    const SPEX_matrix* A // Input matrix (must be SPD)
);

/* Purpose: This function computes the reach of the kth row of A onto the graph of L using the 
   elimination tree. This is more efficient than the SPEX_reach function 
   It finds the nonzero pattern of row k of L and uses the upper triangular 
   part of A(:,k) */
SPEX_info spex_Chol_ereach 
(
    int64_t* top_handle,
    const SPEX_matrix* A,    // Matrix to be analyzed
    int64_t k,          // Node to start at
    int64_t* parent,    // ELimination Tree
    int64_t* xi,         // Contains the nonzero pattern in s[top..n-1] //TODO s to xi (propagate)
    int64_t* w          // Workspace array
);

/* Purpose: Depth-first search and postorder of a tree rooted at node j */
int64_t spex_Chol_tdfs 
(
    int64_t j,      // Root node
    int64_t k,      // TODO what is this?
    int64_t* head,  // Head of list
    int64_t* next,  // Next node in the list
    int64_t* post,  // Post ordered tree
    int64_t* stack  // Stack of nodes
);

/* Purpose: post order a forest */
SPEX_info spex_Chol_post 
(
    int64_t** post_handle, // on input is null on output is post-order of the forest
    int64_t* parent,    // Parent[j] is parent of node j in forest
    int64_t n           // Number of nodes in the forest
);


/* Purpose: consider A(i,j), node j in ith row subtree and return lca(jprev,j) 
   Used to determine Column counts of cholesky factor*/
SPEX_info spex_Chol_leaf //TODO comment everything 
(
    int64_t* lca_handle,
    int64_t i, 
    int64_t j, 
    int64_t* first, 
    int64_t* maxfirst, 
    int64_t* prevleaf,
    int64_t* ancestor, 
    int64_t* jleaf
);

/*Purpose: Obtain the column counts of an SPD matrix for Cholesky factorization
 * This is a modified version of Csparse's cs_chol_counts function
 */
SPEX_info spex_Chol_counts //TODO comment everything 
(
    int64_t** c_handle,
    const SPEX_matrix *A, 
    int64_t* parent, 
    int64_t* post
);

/* Purpose: This function performs the symmetric sparse REF triangular solve. for uplooking
 * Cholesky factorization. i.e., 
 * (LD) x = A(1:k-1,k). 
 * At the given iteration k it computes the k-th column of L' (k-th row of L)
 */
SPEX_info spex_Up_Chol_triangular_solve // performs the sparse REF triangular solve
(
    int64_t* top_output,               // Output the beginning of nonzero pattern //TODO needs better comment
    SPEX_matrix* L,                    // partial L matrix
    const SPEX_matrix* A,              // input matrix
    int64_t k,                         // iteration of algorithm
    int64_t* xi,                       // nonzero pattern vector
    int64_t* parent,                   // Elimination tree
    int64_t* c,                        // Column pointers
    SPEX_matrix* rhos,                 // sequence of pivots
    int64_t* h,                        // history vector
    SPEX_matrix* x                     // solution of system ==> kth row of L
);


/* Purpose: This solves the system L'x = b for Cholesky factorization */
SPEX_info spex_Chol_ltsolve 
(
    const SPEX_matrix* L,    // The lower triangular matrix
    SPEX_matrix* x           // Solution vector
);

/* Purpose: TODO
 */
SPEX_info spex_Chol_Pre_Left_Factor         // performs the Up looking Cholesky factorization
(
    const SPEX_matrix* A,
    SPEX_matrix** L_handle,              // partial L matrix //TODO this is output, reorganize
    int64_t* xi,                  // nonzero pattern vector
    int64_t* parent,              // Elimination tree
    SPEX_Chol_analysis* S,           // stores guess on nnz and column permutation //TODO more descriptive
    int64_t* c                   // Column pointers
);
//TODO think about combining Pre_Left_Factor and Left_Chol_triangular_solve
/* Purpose: This function performs the symmetric sparse REF triangular solve. i.e., 
 * (LD) x = A(:,k). 
 */
SPEX_info spex_Left_Chol_triangular_solve // performs the sparse REF triangular solve
(
    int64_t* top_output,        // Output the beginning of nonzero pattern //TODO better comment
    SPEX_matrix* L,              // partial L matrix
    const SPEX_matrix* A,              // input matrix
    int64_t k,                    // iteration of algorithm
    int64_t* xi,                  // nonzero pattern vector
    SPEX_matrix* rhos,              // sequence of pivots
    int64_t* h,                   // history vector
    SPEX_matrix* x,                  // solution of system ==> kth column of L and U
    int64_t* parent,            //TODO comment
    int64_t* c          //TODO comment
);

/*Purpose: TODO write
*/
SPEX_info spex_Chol_forward_sub
(
    const SPEX_matrix* L,        // lower triangular matrix
    SPEX_matrix* x,              // right hand side matrix of size n*numRHS
    const SPEX_matrix* rhos      // sequence of pivots used in factorization
);

#endif

