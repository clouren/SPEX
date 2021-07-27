//------------------------------------------------------------------------------
// SPEX_Backslash/Include/SPEX_Backslash.h: user #include file for SPEX_Backslash
//------------------------------------------------------------------------------

// SPEX_Backslash: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#ifndef SPEX_BACKSLASH_H
#define SPEX_BACKSLASH_H

// This software package exactly solves sparse systems of linear equations 
// using integer-preserving factorization.
// In essence, this code can be thought of as a wrapper for the other SPEX
// packages.

// Current functionality:
//  -Exactly solve unsymmetric linear systems
//  -Exactly solve symmetric positive definite linear systems

// As more packages are added to SPEX, the above functionality will be expanded

//    To use this code you must first download and install:
//      -GMP 
//      -MPFR
//      -SuiteSparseConfig 
//      -SPEX_Util
//      -SPEX_Left_LU
//      -SPEX_Cholesky
//
//    GMP and MPFR can be found at:
//              https://gmplib.org/
//              http://www.mpfr.org/
//
//   The SPEX dependencies, AMD, COLAMD, and SuiteSparse_Config are distributed along 
//   with SPEX_Backslash. The easiest way to satisfy those dependencies is to just do
//   a top level make in the SPEX package.
//
//   All of these codes are components of the SPEX software library. This code 
//   may be found at:
//              https://github.com/clouren/spex
//              www.suitesparse.com
//  


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Authors----------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Christopher Lourenco, Jinhao Chen, Lorena Mejia Domenzain, 
//    Erick Moreno-Centeno, Timothy Davis

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

//    SPEX_Backslash is free software; you can redistribute it and/or modify
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
// This software is copyright by Christopher Lourenco, Jinhao Chen, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, and Timothy A. Davis. 
// All Rights Reserved.
//

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------------DISCLAIMER-----------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// SPEX Backslash is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//--------------------------Summary---------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// SPEX_backslash is a wrapper for the exact routines contained within the 
// SPEX software package.

// SPEX_BACKSLASH: solve Ax=b via sparse integer-preserving factorization.
// SPEX_backslash: computes the exact solution to the sparse linear system Ax =
// b. A and b may be stored as either int64, double precision, 
// arbitrary precision floating point (mpfr_t), arbitrary sized integer (mpz_t), or 
// arbitrary size rational numbers (mpq_t).
// The result x is computed exactly, represented in arbitrary-precision rational values.
// This solution vector may be returned in either this rational form, or in double
// precision or in arbitrary precision floating point. 
//
// A must be square. If A is SPD, an exact up-looking Cholesky factorization is applied.
// Otherwise, an exact left-looking LU functionality is applied.

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

// SuiteSparse headers
#include "SuiteSparse_config.h"
#include "colamd.h"
#include "amd.h"

// GMP and MPFR
#include <gmp.h>
#include <mpfr.h>

// SPEX Utility functions
#include "SPEX_Util.h"
#include "spex_util_internal.h"


// SPEX Left LU
#include "SPEX_Left_LU.h"

// SPEX Cholesky
#include "SPEX_Chol.h"



// Current version of the code
#define SPEX_BACKSLASH_VERSION "0.0.1"
#define SPEX_BACKSLASH_VERSION_MAJOR 0
#define SPEX_BACKSLASH_VERSION_MINOR 0
#define SPEX_BACKSLASH_VERSION_SUB   1

/* Purpose: Solve Ax = b by analyzing the input matrix
 * and applying the appropiate factorization approach
 */
SPEX_info SPEX_Backslash
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

/* Purpose: Modify the options struct based on the chosen factorization
 */
SPEX_info spex_backslash_set_defaults
(
    SPEX_options* option,
    bool lu
);
    
#endif
