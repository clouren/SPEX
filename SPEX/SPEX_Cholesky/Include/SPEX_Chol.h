//------------------------------------------------------------------------------
// SPEX_Cholesky/Include/SPEX_Cholesky.h: user #include file for SPEX_Cholesky
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

#ifndef SPEX_CHOLESKY_H
#define SPEX_CHOLESKY_H


// This software package exactly solves a sparse symmetric positive definite (SPD)
// system of linear equations using one of two Integer-Preserving Cholesky factorizations. 
// This code accompanies the paper (to be submitted to ACM TOMs)

//    "Algorithm 1XXX: Exactly Solving Sparse Symmetric Positive Definite Linear Systems
//     via SPEX Cholesky factorization," C. Lourenco, E. Moreno-Centeno, T. Davis,
//     to be submitted ACM TOMS.

//   The theory associated with this paper is found at:

//    "A Sparse Integer-Preserving Cholesky Factorization Computed in Time Proportional
//     to Arithmetic Work", C. Lourenco, E. Moreno-Centeno, under submission, SIMAX.

//    If you use this code, you must first download and install the GMP, 
//    MPFR, SPEX_Left_LU, AMD, and COLAMD libraries. 
//   
//   GMP and MPFR can be found at:
//              https://gmplib.org/
//              http://www.mpfr.org/
//
//   SPEX_Left_LU, AMD, and COLAMD are distributed along with SPEX_Cholesky. The easiest
//   way ensure these dependencies are met is to only access this package through 
//   the SPEX repository.
//
//   All of these codes are components of the SPEX software library. This code may
//   be found at:
//              https://github.com/clouren/spex
//              www.suitesparse.com
//  
//


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Authors----------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Christopher Lourenco, Erick Moreno-Centeno, Timothy Davis, Jinhao Chen

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
// This software is copyright by Christopher Lourenco, Erick
// Moreno-Centeno, Timothy A. Davis and Jinhao Chen. All Rights Reserved.
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
//    property of this package is that it can exactly solve ALL SPD input systems. 
//    The input matrix and right hand side vectors are stored as either integers, double
//    precision numbers, multiple precision floating points (through the mpfr
//    library) or as rational numbers (as a collection of numerators and
//    denominators using the GMP mpq_t data structure). Appropriate routines
//    within the code transform the input into an integral matrix in compressed
//    column form.

//    This package computes the factorization PAP' = LDL'. Note that we store the
//    "functional" form of the factorization by only storing the matrix L. The user
//    is given some freedom to select the permutation matrix P. The
//    recommended default settings select P using the AMD ordering.
//    Alternative strategies allowed to select P include the COLAMD 
//    ordering or no column permutation (P=I).  

//    The factor L is computed via integer preserving operations via
//    integer-preserving Gaussian elimination. The key part of this algorithm
//    is a REF Sparse triangular solve function which exploits sparsity and symmetry to
//    reduce the number of operations that must be performed.

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


// Here are the things to do
// TODO: Replace all prints with SPEX_print
// TODO: Documentation
// TODO: Tcov
// TODO: Change outputs to SPEX_info
// TODO: Make internal functions and SPEX_Chol_internal
// TODO: Put const where appropriate

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------Include files required by SPEX Cholesly------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include "SPEX_Util.h"

// SuiteSparse headers
#include "SuiteSparse_config.h"
#include "colamd.h"
#include "amd.h"

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Default Parameters-----------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Current version of the code
#define SPEX_CHOL_VERSION "1.0.0"
#define SPEX_CHOL_VERSION_MAJOR 1
#define SPEX_CHOL_VERSION_MINOR 0
#define SPEX_CHOL_VERSION_SUB   0



//------------------------------------------------------------------------------
// SPEX_Chol_analysis is the data structure for symbolic analysis in the SPEX_Cholesky 
// factorizations. It includes row permutation, elimination tree, and column
// pointers.
//------------------------------------------------------------------------------

typedef struct SPEX_Chol_analysis
{
    int64_t* pinv;      // Row permutation
    int64_t* parent;    // Elimination tree for Cholesky
    int64_t* cp;        // Column pointers for Cholesky
    int64_t lnz;        // Number of nonzeros in Cholesky L
} SPEX_Chol_analysis;
    
    
/* Purpose: Free the Sym_chol data structure */
void SPEX_Chol_analysis_free
(
    SPEX_Chol_analysis* S
);
   
/* Purpose: Permute the matrix A and return A2 = PAP */
SPEX_info SPEX_Chol_permute_A
(
    SPEX_matrix **A2_handle,// Output permuted matrix
    SPEX_matrix* A,        // Initial input matrix
    int64_t* pinv,             // Row permutation
    SPEX_LU_analysis* S    // Column permutation
);

SPEX_info SPEX_Chol_Solve               //solves the linear system LD^(-1)L' x = b
(
    // Output
    SPEX_matrix** x_handle,     // rational solution to the system
    // Input
    SPEX_matrix *A,             // Input matrix (permuted)
    SPEX_matrix* A_orig,        // Input matrix (unpermuted)
    SPEX_matrix* b,             // right hand side vector
    SPEX_matrix* rhos,          // sequence of pivots
    SPEX_matrix* L,             // lower triangular matrix
    int64_t* pinv,                  // row permutation
    SPEX_LU_analysis* S,
    SPEX_options* option        // command options
);

SPEX_info SPEX_Chol_Factor           // performs the Up lookint64_t*g Cholesky factorization
(
    SPEX_matrix* A,             // matrix to be factored
    SPEX_matrix** L_handle,     // lower triangular matrix
    SPEX_Chol_analysis * S,               // stores guess on nnz and column permutation
    SPEX_matrix ** rhos_handle, // sequence of pivots
    bool left,                  // Set true if performing a left-looking factorization
    SPEX_options* option        // command options
);


SPEX_info SPEX_Chol_analyze
(
    SPEX_LU_analysis** S_handle, // symbolic analysis (row/column perm. and nnz L,U)
    const SPEX_matrix *A,        // Input matrix
    const SPEX_options *option   // Control parameters, if NULL, use default
);

SPEX_info SPEX_Chol_backslash
(
    // Output
    SPEX_matrix **X_handle,       // Final solution vector
    // Input
    SPEX_type type,               // Type of output desired
                                  // Must be SPEX_MPQ, SPEX_MPFR, or SPEX_FP64
    const SPEX_matrix *A,         // Input matrix
    const SPEX_matrix *b,         // Right hand side vector(s)
    const SPEX_options* option    // Command options
);




#endif
