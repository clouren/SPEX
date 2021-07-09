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

//TODO "separate" inputs and outputs in every functions DONE
//TODO group similar functions together, check order DONE


/* Purpose: This function computes the reach of the kth row of A onto the graph of L using the 
   elimination tree. This is more efficient than the SPEX_reach function 
   It finds the nonzero pattern of row k of L and uses the upper triangular 
   part of A(:,k) */
SPEX_info spex_Chol_ereach 
(
    // Output
    int64_t* top_handle,
    int64_t* xi,            // Contains the nonzero pattern in xi[top..n-1] //TODO s to xi (propagate) DONE
    // Input
    const SPEX_matrix* A,   // Matrix to be analyzed
    int64_t k,              // Node to start at
    int64_t* parent,        // ELimination Tree
    int64_t* w               // Workspace array
);

/* Purpose: Depth-first search and postorder of a tree rooted at node j */
int64_t spex_Chol_tdfs 
(
    int64_t j,      // Root node
    int64_t k,      // Index (kth node) TODO what is this? DONE?? //TOASK is it output? since it is returned?? but the in k and out k are different and they are ints, so idk if the output thing makes much sense 
    int64_t* head,  // Head of list
    int64_t* next,  // Next node in the list
    int64_t* post,  // Post ordered tree
    int64_t* stack  // Stack of nodes
);


/* Purpose: consider A(i,j), node j in ith row subtree and return lca(jprev,j) 
   Used to determine Column counts of cholesky factor*/
SPEX_info spex_Chol_leaf //TODO comment everything //TOASK maxfirst is both in and out? //TODO still needs to be reordered
(
    //Output
    int64_t* lca_handle,    // Least common ancestor (jprev,j) 
    //Input
    int64_t i,              // Index (subtree i)
    int64_t j,              // Index (node j)
    int64_t* first,         // first[j] is the first descendant of node j
    int64_t* maxfirst,      // maxfirst[j] is the maximum first descendant of node j
    int64_t* prevleaf,      // prevleaf[i] is the previous leaf of ith subtree 
    int64_t* ancestor,      // ancestor[i] is the ancestor of ith subtree
    int64_t* jleaf          // indicates whether j is the first leaf (value of 1) or not (value of 2) //output
);

/* Purpose: Compute the elimination tree of A */
SPEX_info spex_Chol_etree 
(
    // Output
    int64_t** tree,      // Elimination tree of A
    // Input
    const SPEX_matrix* A // Input matrix (must be SPD)
);

/* Purpose: post order a forest */
SPEX_info spex_Chol_post 
(
    // Output
    int64_t** post_handle, // On input is NULL. On output is post-order of the forest
    // Input
    int64_t* parent,    // Parent[j] is parent of node j in forest
    int64_t n           // Number of nodes in the forest
);

/*Purpose: Obtain the column counts of an SPD matrix for Cholesky factorization
 * This is a modified version of Csparse's cs_chol_counts function
 */
SPEX_info spex_Chol_counts //TODO comment everything DONE
(
    // Output
    int64_t** c_handle,     // Column counts
    // Input
    const SPEX_matrix *A,   // Input matrix
    int64_t* parent,        // Elimination tree
    int64_t* post           // Post-order of the tree
);

/* Purpose: TODO DONE?
 */
/* Purpose: This function performs a symbolic left-looking factorization
 * It allocates the memory for the L matrix and allocates the individual
 * entries in the matrix.
 */
SPEX_info spex_Chol_Pre_Left_Factor         // performs the Up looking Cholesky factorization
(
    // Output
    SPEX_matrix** L_handle,       // Partial L matrix //TODO this is output, reorganize DONE
    int64_t* xi,                  // Nonzero pattern vector
    // Input
    const SPEX_matrix* A,         // Input Matrix
    int64_t* parent,              // Elimination tree
    SPEX_Chol_analysis* S,        //Symbolic analysis struct that contains column 
                                  //and inverse row permutations, and number of nonzeros in L //TODO more descriptive DONE
    int64_t* c                   // Column pointers
);
//TODO think about combining Pre_Left_Factor and Left_Chol_triangular_solve
/* Purpose: This function performs the symmetric sparse REF triangular solve. i.e., 
 * (LD) x = A(:,k). 
 */
SPEX_info spex_Left_Chol_triangular_solve // performs the sparse REF triangular solve
(
    //Output
    int64_t* top_output,        // On input NULL. On output contains the beginning of nonzero pattern
                                // The nonzero pattern is contained in xi[top_output...n-1] //TODO better comment DONE
    SPEX_matrix* x,             // Solution of system ==> kth column of L and U
    int64_t* xi,                // Nonzero pattern vector
    // Input
    SPEX_matrix* L,             // Partial L matrix
    const SPEX_matrix* A,       // Input matrix
    int64_t k,                  // Iteration of algorithm
    SPEX_matrix* rhos,          // Sequence of pivots
    int64_t* h,                 // History vector
    int64_t* parent,            // Elimination tree
    int64_t* c                  // Column counts of A
);

/* Purpose: This function performs the symmetric sparse REF triangular solve. for uplooking
 * Cholesky factorization. i.e., 
 * (LD) x = A(1:k-1,k). 
 * At the given iteration k it computes the k-th column of L' (k-th row of L)
 */
SPEX_info spex_Up_Chol_triangular_solve // performs the sparse REF triangular solve
(
    //Output
    int64_t* top_output,               // On input NULL. On output contains the beginning of nonzero pattern
                                       // The nonzero pattern is contained in xi[top_output...n-1] //TODO needs better comment DONE
    int64_t* xi,                       // Nonzero pattern vector
    SPEX_matrix* x,                    // Solution of system ==> kth row of L
    // Input
    SPEX_matrix* L,                    // Partial L matrix
    const SPEX_matrix* A,              // Input matrix
    int64_t k,                         // Iteration of algorithm
    int64_t* parent,                   // Elimination tree
    int64_t* c,                        // Column pointers
    SPEX_matrix* rhos,                 // sequence of pivots
    int64_t* hand                      // History vector
);

/*Purpose: TODO write DONE
*/
/* Purpose: This function performs sparse REF forward substitution This is
 * essentially the same as the sparse REF triangular solve applied to each
 * column of the right hand side vectors. Like the normal one, this
 * function expects that the vector x is dense. As a result,the nonzero
 * pattern is not computed and each nonzero in x is iterated across.
 * The system to solve is LDx = x
 *
 * On output, the matrix x structure is modified
 */
SPEX_info spex_Chol_forward_sub
(
    // Output
    SPEX_matrix* x,              // Right hand side matrix of size n*numRHS
    // Input
    const SPEX_matrix* L,        // Lower triangular matrix
    const SPEX_matrix* rhos      // Sequence of pivots used in factorization
);

/* Purpose: This solves the system L'x = b for Cholesky factorization */
SPEX_info spex_Chol_ltsolve 
(
    // Output
    SPEX_matrix* x,           // Solution vector
    // Input
    const SPEX_matrix* L    // The lower triangular matrix
);

#endif

