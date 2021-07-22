//------------------------------------------------------------------------------
// SPEX_Chol/spex_chol_internal: include file for internal use in SPEX_Cholesky
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// This file is not intended to be #include'd in user applications.  Use
// SPEX_Chol.h instead.

#ifndef SPEX_CHOL_INTERNAL_H
#define SPEX_CHOL_INTERNAL_H

// Definition of SPEX macros, SPEX data structures, etc
#include "spex_util_internal.h"
// SPEX Chol user callable routines
#include "SPEX_Chol.h"

// ============================================================================
//                           Internal Functions
// ============================================================================

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------Routines to compute and anayze the elimination tree------------------
// ----These routines are taken and lightly modified from Tim Davis' Csparse----
// -------------------------www.suitesparse.com---------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Routines to compute and postorder the etree
//------------------------------------------------------------------------------

/* Purpose: Compute the elimination tree of A */
SPEX_info spex_Chol_etree 
(
    // Output
    int64_t** tree,      // On output: contains the elimination tree of A
                         // On input: undefined.
    // Input
    const SPEX_matrix* A // Input matrix (must be SPD). Note to compute
);

/* Purpose: post order a forest */
//TODO: const where appropiate
SPEX_info spex_Chol_post 
(
    // Output
    int64_t** post_handle, // On output: post-order of the forest
                           // On input: undefied
    // Input
    const int64_t* parent,       // Parent[j] is parent of node j in forest
    const int64_t n              // Number of nodes in the forest
);

/* Purpose: Depth-first search and postorder of a tree rooted at node j */
//TODO: const where appropiate
//TODO: This needs to return SPEX_INFO; specifically it now returns -1 if inputs are wrong, but that whouls be fixed to the appropraite SPEX_INFO
SPEX_info spex_Chol_tdfs 
(
    //TODO: Change to int64_t* k (so that it is both input and output). Remember, when calling this function pass in &k (and not k) DONE
    int64_t* k,     // Index (kth node) 
    const int64_t j,// Root node
    int64_t* head,  // Head of list
    int64_t* next,  // Next node in the list
    int64_t* post,  // Post ordered tree
    int64_t* stack  // Stack of nodes
);
 
//------------------------------------------------------------------------------
// Routines to compute the column counts (number of nonzeros per column) of L
//------------------------------------------------------------------------------

/* Purpose: consider A(i,j), node j in ith row subtree and return lca(jprev,j) 
   Used to determine Column counts of cholesky factor*/
//TODO: const where appropiate
SPEX_info spex_Chol_leaf
(
    int64_t* lca_handle,    // Least common ancestor (jprev,j) 
    const int64_t i,        // Index (subtree i)
    const int64_t j,        // Index (node j)
    const int64_t* first,   // first[j] is the first descendant of node j
    int64_t* maxfirst,      // maxfirst[j] is the maximum first descendant of node j
    int64_t* prevleaf,      // prevleaf[i] is the previous leaf of ith subtree 
    int64_t* ancestor,      // ancestor[i] is the ancestor of ith subtree
    int64_t* jleaf          // indicates whether j is the first leaf (value of 1) or not (value of 2) //output
);

/* Purpose: Obtain the column counts of an SPD matrix for Cholesky factorization
 * This is a modified version of Csparse's cs_chol_counts function
 */
//TODO: const where appropiate
SPEX_info spex_Chol_counts
(
    // Output
    int64_t** c_handle,     // On ouptut: column counts
                            // On input: undefined
    // Input
    const SPEX_matrix *A,   // Input matrix
    int64_t* parent,        // Elimination tree
    int64_t* post           // Post-order of the tree
);


//------------------------------------------------------------------------------
// Routine to compute the reach (nonzeros of L) using the etree
//------------------------------------------------------------------------------

/* Purpose: This function computes the reach of the kth row of A on the 
 * elimination tree of A.
 * On input, k is the iteration of the algorithm, parent contains the 
 * elimination tree and w is workspace.
 * On output, xi[top_handle..n-1] contains the nonzero pattern of the 
 * kth row of L (or the kth column of L')
 */
//TODO: const where appropiate
SPEX_info spex_Chol_ereach 
(
    // Output
    int64_t* top_handle,    // On output: starting point of nonzero pattern
                            // On input: undefined
    int64_t* xi,            // On output: contains the nonzero pattern in xi[top..n-1] 
                            // On input: undefined
    // Input
    const SPEX_matrix* A,   // Matrix to be analyzed
    int64_t k,              // Node to start at
    const int64_t* parent,  // Elimination tree of A
    int64_t* w              // Workspace array
);


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------Internal REF Chol Factorization Routines-------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------



/* Purpose: This function performs a symbolic left-looking factorization.
 * On input, A is the matrix to be factored, parent contains the elimination tree
 * and S contains the row/column permutations and number of nonzeros in L
 * On output, L_handle is allocated to contain the nonzero pattern of L and 
 * memory for the values.
 */
//TODO: const where appropiate (ALL input are const, just need to change and propagate) DONE
SPEX_info spex_Chol_Pre_Left_Factor
(
    // Output
    SPEX_matrix** L_handle,       // On output: partial L matrix 
                                  // On input: undefined
    // Input
    int64_t* xi,                  // Workspace nonzero pattern vector
    const SPEX_matrix* A,         // Input Matrix
    const int64_t* parent,        // Elimination tree
    const SPEX_Chol_analysis* S,  // Symbolic analysis struct containing the
                                  // number of nonzeros in L, and the
                                  // row/coluimn permutation and its inverse  
    int64_t* c                    // Column pointers
);
// TODO think about combining Pre_Left_Factor and Left_Chol_triangular_solve
// Chris TODO Can probably eliminate the prealllocation in the left factorization

/* Purpose: This function performs the symmetric sparse REF triangular solve. i.e., 
 * (LD) x = A(:,k). 
 */
//TODO: const where appropiate
SPEX_info spex_Left_Chol_triangular_solve
(
    //Output
    int64_t* top_output,        // On output: the beginning of nonzero pattern
                                // The nonzero pattern is contained in xi[top_output...n-1]
                                // On input: undefined
    SPEX_matrix* x,             // On output: solution of system ==> kth column of L and U
    int64_t* xi,                // On output: nonzero pattern vector
    // Input
    const SPEX_matrix* L,       // Partial L matrix
    const SPEX_matrix* A,       // Input matrix
    const int64_t k,            // Iteration of algorithm
    const SPEX_matrix* rhos,    // Sequence of pivots
    int64_t* h,                 // History vector
    const int64_t* parent,      // Elimination tree
    int64_t* c                  // Column counts of A
);

/* Purpose: This function performs the symmetric sparse REF triangular solve. for uplooking
 * Cholesky factorization. i.e., 
 * (LD) x = A(1:k-1,k). 
 * At the given iteration k it computes the k-th column of L' (k-th row of L)
 */
SPEX_info spex_Up_Chol_triangular_solve 
(
    //Output
    int64_t* top_output,     // On output: the beginning of nonzero pattern
                             // The nonzero pattern is contained in xi[top_output...n-1]
                             // On input: undefined
    int64_t* xi,             // On output: Nonzero pattern vector
    SPEX_matrix* x,          // On output: solution of system ==> kth row of L
    // Input
    const SPEX_matrix* L,    // Partial L matrix
    const SPEX_matrix* A,    // Input matrix
    const int64_t k,         // Iteration of algorithm
    const int64_t* parent,   // Elimination tree
    int64_t* c,              // Column pointers
    const SPEX_matrix* rhos, // sequence of pivots
    int64_t* h               // History vector
);


/* Purpose: This function performs sparse REF forward substitution for Cholesky
 * factorization. 
 * On input, x contains the righ hand side vectors, L is the Cholesky factor of A
 * and rhos is the sequence of pivots used during factorization. 
 * On output, x contains the solution to LD x = x
 * Note that this function assumes that x is stored as a dense matrix
 */
//TODO: const where appropiate
SPEX_info spex_Chol_forward_sub
(
    // Input/Output
    SPEX_matrix* x,              // Right hand side matrix. 
                                 // On input: contains b
                                 // On output: contains the solution of LD x = x
    // Input
    const SPEX_matrix* L,        // REF Cholesky factor of A (lower triangular)
    const SPEX_matrix* rhos      // Sequence of pivots used in factorization
);

/* Purpose: This solves the system L'x = b for Cholesky factorization 
 * On input, x contains the scaled solution of L D x = b and L is the
 * REF Cholesky factor of A.
 * On output, x is the solution to the linear system Ax = (det A)b.
 */
//TODO: const where appropiate
SPEX_info spex_Chol_backward_sub
(
    // Input/Output
    SPEX_matrix* x,         // Solution vector
    // Input
    const SPEX_matrix* L    // REF Cholesky factor of A
);

#endif
