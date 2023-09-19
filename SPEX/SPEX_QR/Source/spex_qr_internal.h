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

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------Internal REF QR Analysis Routines--------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

/* Purpose: Matrix preordering for integer-preserving QR factorization. */
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


/* Purpose: Permute the matrix A and return AQ = AQ */
SPEX_info spex_qr_permute_A
(
    //Output
    SPEX_matrix* AQ_handle,   // On input: undefined
                               // On output: contains the permuted matrix
    //Input
    const SPEX_matrix A,       // Input matrix
    const bool numeric,        // True if user wants to permute pattern and
                               // numbers, false if only pattern
    const int64_t *Q_perm,     // column permutation
    const SPEX_options option  // Command options (Default if NULL)    
);

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------Routines to compute and anayze the elimination tree------------------
// ----These routines are taken and lightly modified from Tim Davis' Csparse----
// -------------------------www.suitesparse.com---------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

/* Purpose: Compute the column elimination tree of A */

SPEX_info spex_qr_etree
(
    // Output
    int64_t **tree_handle,      // On output: contains the elimination tree of A
                                // On input: undefined.
    // Input
    const SPEX_matrix A         // Input matrix (must be SPD).
);

/* Purpose: Obtain the column counts for QR factorization */
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

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------Internal REF QR Factorization Routines-------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

/* Purpose: Obtain the nonzero structure of Q and R for QR factorization */
SPEX_info spex_qr_nonzero_structure
(
    // Output
    SPEX_matrix *R_handle,        // On output: partial R matrix
                                  // On input: undefined
    SPEX_matrix *Q_handle,        // On output: partial R matrix
                                  // On input: undefined
    // Input
    const SPEX_matrix A,          // Input Matrix
    const SPEX_symbolic_analysis S, // Symbolic analysis struct containing the
                                  // number of nonzeros in L, the elimination
                                  // tree, the row/coluimn permutation and its
                                  // inverse
    const SPEX_options option     // Command options
);

/* Purpose: Perfmorm one interation of IPGS-QR.
 * Computes one row of R and updates n-j columns of Q (finalizing the j+1th column)*/
SPEX_info spex_qr_ipgs
(
    //Input/Output
    SPEX_matrix R,       // Right triangular matrix
    SPEX_matrix Q,       // Pair-wise orthogonal matrix
    SPEX_matrix rhos,    // sequence of pivots
    int64_t *Qj,         // pointers to elements of the jth column of Q
    int64_t *h,          // History vector
    //Output
    bool *isZeros,       // True if j+1th column of Q is linearly dependent
    //Input
    const int64_t j,     // Row of R to compute (col j+1 of Q will be finalized)
    const SPEX_matrix A, // Matrix to be factored
    const int64_t *Q_perm,     // Column permutation
    const SPEX_options option  // Command options
);


SPEX_info spex_qr_back_sub  // performs sparse REF backward substitution
(
    SPEX_matrix bx,         // right hand side matrix
    const SPEX_matrix R,   // input upper triangular matrix
    const int64_t rank     // rank of right triangular matrix
);



//future TODO merge qr transpose and normal utilities transpose
SPEX_info spex_qr_transpose
(
    SPEX_matrix *C_handle,      // C = A'
    SPEX_matrix A,              // Matrix to be transposed
    const SPEX_options option   // Command options
);

#endif
