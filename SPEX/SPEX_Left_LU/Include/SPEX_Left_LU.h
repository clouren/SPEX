//------------------------------------------------------------------------------
// SPEX_Left_LU/Include/SPEX_Left_LU.h: user #include file for SPEX_Left_LU.
//------------------------------------------------------------------------------

// SPEX_Left_LU: (c) 2019-2022, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#ifndef SPEX_LEFT_LU_H
#define SPEX_LEFT_LU_H

// This software package exactly solves a sparse system of linear equations
// using the SPEX Left LU factorization. This code accompanies the paper (submitted
// to ACM Transactions on Mathematical Software):

//    "Algorithm 1xxx: SPEX Left LU: Exactly Solving Sparse Linear Systems via
//    A Sparse Left-Looking Integer-Preserving LU Factorization",
//    C. Lourenco, J. Chen, E. Moreno-Centeno, T. Davis, under submission,
//    ACM Trans. Mathematical Software.

//    The theory associated with this software can be found in the paper
//    (published in SIAM journal on matrix analysis and applications):

//    "Exact Solution of Sparse Linear Systems via Left-Looking
//     Roundoff-Error-Free LU Factorization in Time Proportional to
//     Arithmetic Work", C. Lourenco, A. R. Escobedo, E. Moreno-Centeno,
//     T. Davis, SIAM J. Matrix Analysis and Applications.  pp 609-638,
//     vol 40, no 2, 2019.

//    If you use this code, you must first download and install the GMP and
//    MPFR libraries. GMP and MPFR can be found at:
//              https://gmplib.org/
//              http://www.mpfr.org/

//    If you use SPEX Left LU for a publication, we request that you please cite
//    the above two papers.

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Authors----------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Christopher Lourenco, Jinhao Chen, Erick Moreno-Centeno, and Timothy Davis
//

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Contact Information----------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Please contact Chris Lourenco (chrisjlourenco@gmail.com)
//    or Tim Davis (timdavis@aldenmath.com, DrTimothyAldenDavis@gmail.com,
//                  davis@tamu.edu)

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Copyright--------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    SPEX Left LU is free software; you can redistribute it and/or modify
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
// This software is copyright by Christopher Lourenco, Jinhao Chen, Erick
// Moreno-Centeno and Timothy A. Davis. All Rights Reserved.
//

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------------DISCLAIMER-----------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// SPEX Left LU is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//--------------------------Summary---------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    This software package solves the linear system Ax = b exactly. The input
//    matrix and right hand side vectors are stored as either integers, double
//    precision numbers, multiple precision floating points (through the mpfr
//    library) or as rational numbers (as a collection of numerators and
//    denominators using the GMP mpq_t data structure). Appropriate routines
//    within the code transform the input into an integral matrix in compressed
//    column form.

//    This package computes the factorization PAQ = LDU. Note that we store the
//    "functional" form of the factorization by only storing L and U. The user
//    is given some freedom to select the permutation matrices P and Q. The
//    recommended default settings select Q using the COLAMD column ordering
//    and select P via a partial pivoting scheme in which the diagonal entry
//    in column k is selected if it is the same magnitude as the smallest
//    entry, otherwise the smallest entry is selected as the kth pivot.
//    Alternative strategies allowed to select Q include the AMD column
//    ordering or no column permutation (Q=I).  For pivots, there are a variety
//    of potential schemes including traditional partial pivoting, diagonal
//    pivoting, tolerance pivoting etc. This package does not allow pivoting
//    based on sparsity criterion.

//    The factors L and U are computed via integer preserving operations via
//    integer-preserving Gaussian elimination. The key part of this algorithm
//    is a Roundoff Error Free (REF) sparse triangular solve function which
//    exploits sparsity to reduce the number of operations that must be
//    performed.

//    Once L and U are computed, a simplified version of the triangular solve
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
//---------------------Include files required by SPEX Left LU-------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include "SuiteSparse_config.h"
#include "SPEX_Util.h"

//------------------------------------------------------------------------------
// Version
//------------------------------------------------------------------------------

// Current version of the code
#define SPEX_LEFT_LU_VERSION "1.1.1"
#define SPEX_LEFT_LU_VERSION_MAJOR 1
#define SPEX_LEFT_LU_VERSION_MINOR 1
#define SPEX_LEFT_LU_VERSION_SUB   1



//------------------------------------------------------------------------------
// Primary factorization & solve routines
//------------------------------------------------------------------------------

// SPEX_backslash solves the linear system Ax = b. This is the simplest way to
// use the SPEX Left LU package. This function encompasses both factorization and
// solve and returns the solution vector in the user desired type.  It can be
// thought of as an exact version of MATLAB sparse backslash.
SPEX_info SPEX_Left_LU_backslash
(
    // Output
    SPEX_matrix **X_handle,       // Final solution vector
    // Input
    SPEX_type type,               // Type of output desired:
                                  // Must be SPEX_MPQ, SPEX_MPFR,
                                  // or SPEX_FP64
    const SPEX_matrix *A,         // Input matrix
    const SPEX_matrix *b,         // Right hand side vector(s)
    const SPEX_options* option
) ;

// SPEX_Left_LU_factorize performs the SPEX Left LU factorization. 
// This factorization is done via n iterations of the sparse REF 
// triangular solve function. The overall factorization is 
// PAQ = LDU.  The determinant can be obtained as rhos->x.mpz[n-1].
// 
//  L: undefined on input, created on output
//  U: undefined on input, created on output
//  rhos: undefined on input, created on output
//  pinv: undefined on input, created on output
// 
//  A: input only, not modified
//  S: input only, not modified
//  option: input only, not modified
SPEX_info SPEX_Left_LU_factorize
(
    // output:
    SPEX_matrix **L_handle,     // lower triangular matrix
    SPEX_matrix **U_handle,     // upper triangular matrix
    SPEX_matrix **rhos_handle,  // sequence of pivots
    int64_t **pinv_handle,      // inverse row permutation
    // input:
    const SPEX_matrix *A,       // matrix to be factored
    const SPEX_LU_analysis *S,  // column permutation and estimates
                                // of nnz in L and U 
    const SPEX_options* option
) ;

// SPEX_Left_LU_solve solves the linear system LD^(-1)U x = b.
SPEX_info SPEX_Left_LU_solve         // solves the linear system LD^(-1)U x = b
(
    // Output
    SPEX_matrix **X_handle,     // rational solution to the system
    // input:
    const SPEX_matrix *b,       // right hand side vector
    const SPEX_matrix *A,       // Input matrix
    const SPEX_matrix *L,       // lower triangular matrix
    const SPEX_matrix *U,       // upper triangular matrix
    const SPEX_matrix *rhos,    // sequence of pivots
    const SPEX_LU_analysis *S,  // symbolic analysis struct
    const int64_t *pinv,        // inverse row permutation
    const SPEX_options* option
) ;

#endif

