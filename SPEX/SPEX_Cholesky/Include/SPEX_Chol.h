//------------------------------------------------------------------------------
// SPEX_Cholesky/Include/SPEX_Cholesky.h: user #include file for SPEX_Cholesky
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//TODO propagate this to every single file :D DOUBLECHECK

//------------------------------------------------------------------------------

#ifndef SPEX_CHOLESKY_H
#define SPEX_CHOLESKY_H


// This software package exactly solves a sparse symmetric positive definite 
// (SPD) system of linear equations using one of two Integer-Preserving Cholesky  
// factorizations. This code accompanies the paper (to be submitted to ACM TOMs)

//    "Algorithm 1xxx: Exactly Solving Sparse Symmetric Positive Definite Linear 
//     Systems via SPEX Cholesky factorization," C. Lourenco, L. Mejia Domenzain,
//     E. Moreno-Centeno, T. Davis, to be submitted ACM TOMS. 

//     The theory associated with this paper is found at:

//    "Exactly Solving Sparse Rational Linear Systems via Roundoff-Error-Free 
//     Cholesky Factorizations", C. Lourenco, E. Moreno-Centeno, 
//     under submission, SIMAX. 

//    To use this code you must first download and install the GMP, 
//    MPFR, SPEX_Util, AMD, and COLAMD libraries. GMP and MPFR can be found at:
//              https://gmplib.org/
//              http://www.mpfr.org/
//
//   SPEX_Util, AMD, and COLAMD are distributed along with SPEX_Cholesky. 
//   The easiest way ensure these dependencies are met is to only access this 
//   package through the SPEX repository.
//
//   All of these codes are components of the SPEX software library. This code 
//   may be found at:
//              https://github.com/clouren/spex
//              www.suitesparse.com
//  
//


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Authors----------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Christopher Lourenco, Lorena Mejia Domenzain, Erick Moreno-Centeno, 
//    Timothy Davis

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Contact Information----------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Please contact Chris Lourenco (chrisjlourenco@gmail.com, lourenco@usna.edu)
//    or Tim Davis (timdavis@aldenmath.com, DrTimothyAldenDavis@gmail.com,
//                  davis@tamu.edu)


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Copyright--------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    SPEX_Cholesky is free software; you can redistribute it and/or modify
//     it under the terms of either:
//
//        * the GNU Lesser General Public License as published by the
//          Free Software Foundation; either version 3 of the License,
//          or (at your option) any later version.
//
//     or
//
//        * the GNU General Public License as published by the Free Software
//          Foundation; either version 2 of the License, or (at your option) any
//          later version.
//
//    or both in parallel, as here.
//
//    See license.txt for license info.
//
// This software is copyright by Christopher Lourenco, Lorena Mejia Domenzain,
// Erick Moreno-Centeno, and Timothy A. Davis. All Rights Reserved.
//

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------------DISCLAIMER-----------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// SPEX Cholesky is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//--------------------------Summary---------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    This software package solves the SPD linear system Ax = b exactly. The key
//    property of this package is that it can exactly solve any SPD input system. 
//    The input matrix and right hand side vectors are stored as either integers,
//    double precision numbers, multiple precision floating points (through the 
//    mpfr library) or as rational numbers (as a collection of numerators and
//    denominators using the GMP mpq_t data structure). Appropriate routines
//    within the code transform the input into an integral matrix in compressed
//    column form.

//    This package computes the factorization PAP' = LDL'. Note that we store
//    the "functional" form of the factorization by only storing the matrix L. 
//    The user is given some freedom to select the permutation matrix P. The
//    recommended default settings select P using the AMD ordering.
//    Alternative strategies allowed to select P include the COLAMD 
//    ordering or no column permutation (P=I).  

//    The factor L is computed via integer preserving operations via
//    integer-preserving Gaussian elimination. The key part of this algorithm
//    is a REF Sparse triangular solve function which exploits sparsity and 
//    symmetry to reduce the number of operations that must be performed.

//    Once L is computed, a simplified version of the triangular solve
//    is performed which assumes the vector b is dense. The final solution
//    vector x is gauranteed to be exact. This vector can be output in one of
//    three ways: 1) full precision rational arithmetic (as a sequence of
//    numerators and denominators) using the GMP mpq_t data type, 2) double
//    precision while not exact will produce a solution accurate to machine
//    roundoff unless the size of the associated solution exceeds double
//    precision (i.e., the solution is 10^500 or something), 3) variable
//    precision floating point using the GMP mpfr_t data type. The associated
//    precision is user defined.

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------Include files required by SPEX Cholesly------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Standard C libraries
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

// GMP and MPFR
#include <gmp.h>
#include <mpfr.h>

// SPEX Utility functions
#include "SPEX_Util.h"

// SuiteSparse headers
#include "SuiteSparse_config.h"
#include "colamd.h"
#include "amd.h"

// Current version of the code
#define SPEX_CHOL_VERSION "0.0.1"
#define SPEX_CHOL_VERSION_MAJOR 0
#define SPEX_CHOL_VERSION_MINOR 0
#define SPEX_CHOL_VERSION_SUB   1



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------Symbolic Analysis Data structure-------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

typedef struct SPEX_Chol_analysis
{
    int64_t* pinv;      // Inverse row/column permutation
    int64_t* p ;        // Row/column permutation, representing //TODO propagate q->p DONE
                        // the permutation matrix P. The matrix P*A*P' is factorized.
                        // If the kth column of L, and P*A*P' is column j of the
                        // unpermuted matrix A, then j = S->p [k].
    int64_t* parent;    // Elimination tree of A for Cholesky
    int64_t* cp;        // Column pointers of L 
    int64_t lnz;        // Number of nonzeros in Cholesky L (might be estimate).
                        // At initialization, if the default column ordering (AMD) is 
                        // used lnz will be exact otherwise lnz will be an estimate.
                        // After the elimination tree is computed lnz will be exact.
} SPEX_Chol_analysis;

/* Purpose: Free the SPEX_Chol_analysis data structure */
void SPEX_Chol_analysis_free
(
    SPEX_Chol_analysis** S
);

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-----------------------Primary SPEX Cholesky routines-------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------    
    
/* Purpose: Compute the exact solution of Ax = b. 
 * On input, A must be SPD and x is NULL
 * On output, x contains the solution of the linear system
 */
SPEX_info SPEX_Chol_backslash
(
    // Output
    SPEX_matrix** X_handle,       // On input: undefined. 
                                  // On output: final solution vector
    // Input
    SPEX_type type,               // Type of output desired
                                  // Must be SPEX_MPQ, SPEX_MPFR, or SPEX_FP64
    const SPEX_matrix* A,         // Input matrix
    const SPEX_matrix* b,         // Right hand side vector(s)
    const SPEX_options* option    // Command options
);
    
/* Purpose: Matrix preordering for integer-preserving Cholesky factorization.
 * On input, S is undefined
 * On output, S contains the row/column permutation of A
 */
SPEX_info SPEX_Chol_preorder
(
    // Output
    SPEX_Chol_analysis** S_handle,  // Symbolic analysis data structure 
                                    // On input: undefined
                                    // On output: contains the 
                                    // row/column permutation and its
                                    // inverse.
    // Input
    const SPEX_matrix* A,           // Input matrix
    const SPEX_options* option      // Control parameters (use default if NULL)
);

/* Purpose: Permute the matrix A and return PAP = PAP' 
 * On input PAP is undefined and A contains the input matrix
 * On output PAP contains the permuted matrix (PAP')
 */
SPEX_info SPEX_Chol_permute_A
(
    //Output
    SPEX_matrix** PAP_handle,  // On input: undefined // TODO Propagate A2 -> PAP DONE TOASK PAP= PAP' <- looks weird and kind of confusing
                               // On output: contains the permuted matrix
    //Input
    SPEX_matrix* A,            // Input matrix
    SPEX_Chol_analysis* S      // Symbolic analysis struct that contains 
                               // column and inverse row permutations
);

/* Purpose: Compute the REF Cholesky factorization A = LDL'
 * only appropriate if A is SPD. 
 * On input A contains the user's matrix, left is a bool
 * where true indicates a left-looking factorization and false
 * indicates the up-looking factorization
 * On output, L contains the REF Cholesky factor of A, rhos contains
 * the REF Cholesky pivot elements and S contains the elimination tree
 * lower triangular matrix and rhos contains the pivots' values
 * used in the factorization 
 *
 * TODO: Write that left/up is based on SPEX_options once done
 */
SPEX_info SPEX_Chol_Factor     
(
    // Output
    SPEX_matrix** L_handle,    // Lower triangular matrix. NULL on input
    SPEX_matrix** rhos_handle, // Sequence of pivots. NULL on input.
    SPEX_Chol_analysis* S,     // Symbolic analysis struct containing the
                               // elimination tree of A, column pointers of L, 
                               // exact number of nonzeros of L and permutation used
    // Input
    const SPEX_matrix* A,      // Matrix to be factored
    
    // TODO: add an option of option->CHOL_TYPE = CHOL_LEFT or CHOL_UP
    // Also could add option->LU_TYPE = LU_LEFT (in future LU_RIGHT)
    // Default should be CHOL_UP
    // TODO Define CHOL_LEFT, CHOL_UP, LU_LEFT, and expand options in SPEX_UTIL
    // Also, upgrade the comments here
    bool left,                 // Set to true if performing a left-looking factorization; 
                               // otherwise perform an up-looking factorization.
    const SPEX_options* option // Command options
                               // Notably, option->CHOL_TYPE indicates whether
                               // CHOL_UP (default) or CHOL_LEFT is used.
);

/* Purpose: After computing the REF Cholesky factorization A = LDL',
 * this function solves the associated linear system LDL' x = b
 * TODO On input blah on output BLAH DONE
 * On input x is undefined, A contains the user's matrix, b contains 
 * the user's right hand side, rhos contains the pivots' values used 
 * in the factorization, L contains the REF Cholesky factor of A, 
 * and S contains the elimination tree
 * On output x contains the rational solution of the system LDL' x = b
 */
//TODO: Revistit arguments after Jinhaos update. (especially A, A_orig etc.)
SPEX_info SPEX_Chol_Solve
(
    // Output
    SPEX_matrix** x_handle,           // On input: undefined.
                                      // On output: Rational solution (SPEX_MPQ) to the system. 
    // Input
    const SPEX_matrix* A,             // Input matrix (permuted)
    // TODO: A_orig is only used for check, the check can be seperated from this function.
    // So remove A_orig and put it seperately
    const SPEX_matrix* A_orig,        // Input matrix (unpermuted)
    const SPEX_matrix* b,             // Right hand side vector
    const SPEX_matrix* rhos,          // Pivots' values 
    const SPEX_matrix* L,             // Lower triangular matrix
    const SPEX_Chol_analysis* S,      // Symbolic analysis struct. contains elimination tree of A,  
                                      // column pointers of L, exact number of nonzeros of L and permutation used
    const SPEX_options* option        // command options
);




#endif
