//------------------------------------------------------------------------------
// SPEX_Update/Include/SPEX_Update.h: user #include file for SPEX_Update.
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

#ifndef SPEX_UPDATE_H
#define SPEX_UPDATE_H

// This software package exactly updates the Cholesky factorization of a
// sparse matrix. This code accompanies the paper


//    The theory associated with this software can be found in the paper


//    If you use this code, you must first download and install the GMP and
//    MPFR libraries. GMP and MPFR can be found at:
//              https://gmplib.org/
//              http://www.mpfr.org/

//    If you use SPEX UPDATE for a publication, we request that you
//    please cite the above two papers.

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Authors----------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Jinhao Chen, Timothy Davis, and Erick Moreno-Centeno
//

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Contact Information----------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Please contact Jinhao Chen (cjh10644@hotmail.com)
//    or Tim Davis (timdavis@aldenmath.com, DrTimothyAldenDavis@gmail.com,
//                  davis@tamu.edu)

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Copyright--------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    SPEX UPDATE is free software; you can redistribute it and/or modify
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
// This software is copyright by Jinhao Chen, Timothy A. Davis and Erick
// Moreno-Centeno. All Rights Reserved.
//

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------------DISCLAIMER-----------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// SPEX UPDATE is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//--------------------------Summary---------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    This software package update a REF Cholesky factorization A = LDL^T
//    exactly when A is changed with a row and column. The input matrices are
//    stored as either integers, double precision numbers, multiple precision
//    floating points (through the mpfr library) or as rational numbers (as a
//    collection of numerators and denominators using the GMP mpq_t data
//    structure). Appropriate routines within the code transform the input into
//    an integral matrix in compressed column form.


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------Include files required by SPEX UPDATE----------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#include "SPEX_Util.h"

//------------------------------------------------------------------------------
// Version
//------------------------------------------------------------------------------

// Current version of the code
#define SPEX_UPDATE_VERSION "1.0.0"
#define SPEX_UPDATE_VERSION_MAJOR 1
#define SPEX_UPDATE_VERSION_MINOR 0
#define SPEX_UPDATE_VERSION_SUB   0

//------------------------------------------------------------------------------
// LU Update for column replacement
//------------------------------------------------------------------------------

// This function performs LU update for column replacement. The matrices in the
// input factorization can be any type and/or kind and does not have to be in
// updatable format. The function will always first check if the factorization
// is updatable and perform necessary conversion if needed. L and U in the
// output factorization will be updatable. The factorization will be modified
// during the update process.  Therefore, if this function fails for any
// reason, the returned F should be considered as undefined.
//
// The matrix A is not modified during the update. If the updated A is
// needed, user can use the follow code if A is in SPEX_dynamic_CSC form.
//
//       SPEX_vector *Vtmp = A->v[k];
//       A->v[k] = vk->v[0];
//       vk->v[0] = Vtmp;

SPEX_info SPEX_Update_LU_ColRep
(
    SPEX_factorization* F,  // The SPEX factorization of A, including L, U,
                            // rhos, P, Pinv, Q and Qinv. The factorization
                            // will be modified during the update process.
                            // Therefore, if this function fails for any
                            // reason, the returned F should be considered as
                            // undefined.
    SPEX_matrix *vk,        // Pointer to a n-by-1 dynamic_CSC matrix
                            // which contains the column to be inserted.
                            // The rows of vk are in the same order as A.
    int64_t k,              // The column index that vk will be inserted, 0<=k<n
    const SPEX_options *option// Command parameters
);

//------------------------------------------------------------------------------
// Rank-1 Cholesky update/downdate
//------------------------------------------------------------------------------

// This function performs rank-1 Cholersky update/downdate. The matrices in the
// input factorization can be any type and/or kind and does not have to be in
// updatable format. The function will always first check if the factorization
// is updatable and perform necessary conversion if needed. L in the output
// factorization will become updatable.
// 
// The matrix A is not modified during the update. If the updated A is needed,
// user can compute A = A + sigma*w*wT BEFORE using this function (since w
// will be modified).

SPEX_info SPEX_Update_Chol_Rank1
(
    SPEX_factorization *F,  // The SPEX Cholesky factorization of A, including
                            // L, rhos, P and Pinv. This factorization will be
                            // modified during the update process. Therefore,
                            // if this function fails for any reason, the
                            // returned F should be considered as undefined.
    SPEX_matrix *w,         // a n-by-1 dynamic_CSC matrix that contains the
                            // vector to modify the original matrix A, the
                            // resulting A is A+sigma*w*w^T. In output, w is
                            // updated as the solution to L*D^(-1)*w_out = w
    const int64_t sigma,    // a nonzero scalar that determines whether
                            // this is an update or downdate
    const SPEX_options* option // Command options
);

//------------------------------------------------------------------------------
// Function for solving A*x =b via LU or Cholesky factorization
//------------------------------------------------------------------------------

// NOTE: the requirement for L and UT is exactly same as other functions in
// SPEX_Update, which requires the pivot of L->v[j] or UT->v[j] be the
// first entry correspondingly.

SPEX_info SPEX_Update_solve // solves Ax = b via LU or Cholesky factorization
(
    // Output
    SPEX_matrix **x_handle, // a m*n dense matrix contains the solution to
                            // the system.
    // input:
    SPEX_factorization *F,  // The SPEX LU or Cholesky factorization
    const SPEX_matrix *b,   // a m*n dense matrix contains the right-hand-side
                            // vector
    const SPEX_options* option // Command options
);

//------------------------------------------------------------------------------
// Function for solving A^T*x =b via LU or Cholesky factorization of A
//------------------------------------------------------------------------------
SPEX_info SPEX_Update_tsolve // solves A^T*x = b
(
    // Output
    SPEX_matrix **x_handle, // a m*n dense matrix contains the solution to
                            // the system.
    // input:
    SPEX_factorization *F,  // The SPEX LU or Cholesky factorization of A
    const SPEX_matrix *b,   // a m*n dense matrix contains the right-hand-side
                            // vector
    const SPEX_options* option // Command options
);

//------------------------------------------------------------------------------
// Function for converting factorization between updatable (with L and/or U as
// dynamic_CSC MPZ matrices) and non-updatable (with L and/or U as CSC MPZ
// matrices) formats.
//------------------------------------------------------------------------------

// NOTE: if F->updatable == false upon input, F->L (and F->U if exists) must be
// CSC MPZ matrix, otherwise, SPEX_INCORRECT_INPUT will be returned. Likewise,
// if F->updatable == true upon input, F->L (and F->U if exists) must be
// dynamic_CSC MPZ matrix. In addition, both F->L and F->U (if exists) must not
// be shallow matrices. All SPEX functions output either of these two formats
// and non-shallow. Therefore, these input requirements can be met easily if
// users do not try to modify any individual component of F.  The conversion is
// done in place and F->updatable will be set to its complement upon output. In
// case of any error, the returned factorization should be considered as
// undefined.
/*
SPEX_info SPEX_Update_factorization_convert
(
    SPEX_factorization *F,  // The factorization to be converted
    const SPEX_options* option // Command options
);*/

#endif

