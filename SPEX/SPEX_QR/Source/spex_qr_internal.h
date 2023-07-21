//------------------------------------------------------------------------------
// SPEX_QR/Source/spex_qr_internal: include file for internal use in
// SPEX_QR
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021-2023, Chris Lourenco, Lorena Mejia Domenzain,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// This file is not intended to be #include'd in user applications.  Use
// SPEX.h instead.

#ifndef SPEX_QR_INTERNAL_H
#define SPEX_QR_INTERNAL_H

#include "spex_util_internal.h"
#include "spex_cholesky_internal.h"

// ============================================================================
//                           Internal Functions
// ============================================================================

SPEX_info spex_qr_preorder
(
    // Output
    SPEX_symbolic_analysis *S_handle,   // Symbolic analysis data structure
                                        // On input: undefined
                                        // On output: contains the
                                        // row/column permutation and its
                                        // inverse.
    // Input
    const SPEX_matrix A,            // Input matrix
    const SPEX_options option       // Control parameters (use default if NULL)
);

SPEX_info spex_qr_permute_A
(
    //Output
    SPEX_matrix* PAQ_handle,   // On input: undefined
                               // On output: contains the permuted matrix
    //Input
    const SPEX_matrix A,       // Input matrix
    const bool numeric,        // True if user wants to permute pattern and
                               // numbers, false if only pattern
    const SPEX_symbolic_analysis S,  // Symbolic analysis struct that contains
                                // row/column permutations
    const SPEX_options option        
);

SPEX_info spex_qr_etree
(
    // Output
    int64_t **tree_handle,      // On output: contains the elimination tree of A
                                // On input: undefined.
    // Input
    const SPEX_matrix A         // Input matrix (must be SPD).
);

/*static void spex_qr_init_ata 
(
    SPEX_matrix AT, 
    const int64_t *post, 
    int64_t *w, 
    int64_t**head, 
    int64_t **next
);*/

SPEX_info spex_qr_counts
(
    // Output
    int64_t **c_handle,     // On ouptut: column counts
                            // On input: undefined
    // Input
    const SPEX_matrix A,    // Input matrix
    const int64_t *parent,  // Elimination tree
    const int64_t *post     // Post-order of the tree
);

SPEX_info spex_qr_ipgs
(
    SPEX_matrix R,    
    SPEX_matrix Q,    
    SPEX_matrix rhos,         // sequence of pivots
    const int64_t j,          // Row of R to compute (col j+1 of Q will also be computed)
    const SPEX_matrix A,      // Matrix to be factored
    /*int64_t *h,         // History vector
    int64_t *xi,        // Nonzero pattern vector of row j of R
    int64_t *c,
    const int64_t *parent,          // Elimination tree*/
    int64_t *Prev,
    int64_t *vec,
    const int64_t *leftmost,
    SPEX_options option
);

SPEX_info spex_qr_ipgsM
(
    SPEX_matrix R,    
    SPEX_matrix Q,    
    SPEX_matrix rhos,         // sequence of pivots
    int64_t *Qj,
    int64_t *col,
    const int64_t j,          // Row of R to compute (col j+1 of Q will also be computed)
    const SPEX_matrix A,      // Matrix to be factored
    int64_t *h,
    SPEX_options option
);

SPEX_info spex_qr_pre_factor
(
    // Output
    SPEX_matrix *R_handle,        // On output: partial L matrix
                                  // On input: undefined
    // Input
    //int64_t *xi,                  // Workspace nonzero pattern vector
    const SPEX_matrix A,          // Input Matrix
    const SPEX_symbolic_analysis S  // Symbolic analysis struct containing the
                                  // number of nonzeros in L, the elimination
                                  // tree, the row/coluimn permutation and its
                                  // inverse
);

SPEX_info spex_qr_pre_Q
(
    SPEX_matrix *Q_handle,
    SPEX_matrix A,
    SPEX_options option
);

SPEX_info spex_qr_pre_factorQR
(
    // Output
    SPEX_matrix *R_handle,        // On output: partial R matrix
                                  // On input: undefined
    SPEX_matrix *Q_handle,        // On output: partial R matrix
                                  // On input: undefined
    // Input
    //int64_t *xi,                  // Workspace nonzero pattern vector
    const SPEX_matrix A,          // Input Matrix
    const SPEX_symbolic_analysis S  // Symbolic analysis struct containing the
                                  // number of nonzeros in L, the elimination
                                  // tree, the row/coluimn permutation and its
                                  // inverse
);


SPEX_info spex_qr_transpose
(
    SPEX_matrix *C_handle,      // C = A'
    SPEX_matrix A,              // Matrix to be transposed
    const SPEX_options option
);

#endif
