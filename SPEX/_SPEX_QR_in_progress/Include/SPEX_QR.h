//------------------------------------------------------------------------------
// SPEX_QR/SPEX_QR.h: user #include file for SPEX_QR
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021, Chris Lourenco, US Naval Academy, All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#ifndef SPEX_QR_H
#define SPEX_QR_H


// This software performs an exact integer-preserving QR factorization
// WARNING: This code is experimental and developmental, please do not use it.


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Authors----------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Christopher Lourenco

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

//    SPEX_QR is free software; you can redistribute it and/or modify
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
// This software is copyright by Christopher Lourenco
//

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------------DISCLAIMER-----------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// SPEX QR is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//--------------------------Summary---------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Given a rectangular matrix A, SPEX QR factorizes A into the product of A = QDR
// where Q is orthogonal and R is upper trapezoidal. R is also the REF Cholesky 
// factor of A'*A. 

// SPEX QR can be used to solve rectangular linear systems Ax = b or accurately
// determine the singularity of a matrix.


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------Include files required by SPEX QR------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

// GMP and MPFR
#include <gmp.h>
#include <mpfr.h>

// SPEX
#include "SPEX.h"

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Default Parameters-----------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//TODO Seperate and do the SPEX_CHECK

// Current version of the code
#define SPEX_QR_VERSION "0.0.1"
#define SPEX_CHOL_VERSION_MAJOR 0
#define SPEX_CHOL_VERSION_MINOR 0
#define SPEX_CHOL_VERSION_SUB   1

// Should be defined in SPEX_Utilities
//#define ASSERT assert

/* Compute the dot product of two integer vectors x,y and return in z */
SPEX_info SPEX_dot
(
    SPEX_matrix x,
    SPEX_matrix y,
    mpz_t z
);

/* Purpose: Given a matrix A in m*n and B in m*n, compute the dot product of
 * A(:,i) and B(:,j). Assumed to be dense. prod = A(:,i) dot B(:,j)
 */
SPEX_info SPEX_dense_mat_dot
(
    SPEX_matrix A,
    int64_t i,
    SPEX_matrix B,
    int64_t j,
    mpz_t prod
);
    
/* Perform the IPGE version of SPEX QR (aka Algorithm 1 from workpage)
 */
SPEX_info SPEX_QR_IPGE
(
    SPEX_matrix A,            // Matrix to be factored
    SPEX_matrix *R_handle,    // upper triangular matrix
    SPEX_matrix *Q_handle     // orthogonal triangular matrix
);                                 
    
    
    
/* Perform the IPGE version of SPEX QR using Pursell method
 */
SPEX_info SPEX_QR_PURSELL
(
    SPEX_matrix A,            // Matrix to be factored
    SPEX_matrix *R_handle,    // upper triangular matrix
    SPEX_matrix *Q_handle     // orthogonal triangular matrix
);
  
/* Perform the IPGE version of SPEX QR using Pursell method
 */
SPEX_info SPEX_QR_PURSELL2
(
    SPEX_matrix A,            // Matrix to be factored
    SPEX_matrix *R_handle,    // upper triangular matrix
    SPEX_matrix *Q_handle     // orthogonal triangular matrix
);

SPEX_info SPEX_generate_random_matrix
(
    SPEX_matrix *A_handle, // Matrix to be created
    int64_t m,              // Rows of the matrix
    int64_t n,              // Columns of the matrix
    unsigned int seed,      // Random number seed
    int64_t lower,          // Lower bound for numbers to be generated
    int64_t upper           // Upper bound for numbers to be generated
);

SPEX_info SPEX_Qtb
(
    SPEX_matrix Q,        // Q matrix, want Q'
    SPEX_matrix b,        // Original RHS Vector
    SPEX_matrix* b_handle // Q'*b
);
    

SPEX_info SPEX_QR_backsolve
(
    SPEX_matrix R,        // Upper triangular matrix
    SPEX_matrix b,        // Q^T * b
    SPEX_matrix* x_handle // Solution
);
    

#endif
