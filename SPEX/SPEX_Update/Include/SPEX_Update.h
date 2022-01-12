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

// NOTE: A = LD^(-1)U should hold upon input and will be maintained on sucessful
// output. Additionally, as described in the input below, the first entry of
// L->v[j] and UT->v[j] should be the j-th pivot (use
// SPEX_Update_matrix_canonicalize to ensure this), which must be hold upon
// input and will be maintained on output. vk will be swapped with A->v[k]. If
// the function fails for any reason, L and UT should be considered as
// undefined.

SPEX_info SPEX_Update_LU_ColRep
(
// TODO should I modify A outside this function
    SPEX_matrix *A,         // n-by-n original matrix in dynamic_CSC form
    SPEX_factorization* F,  // The SPEX factorization of A, including L, U,
                            // rhos, P, Pinv, Q and Qinv. The factorization L
                            // and UT must be in dynamic_CSC form. The
                            // factorization will be modified during the update
                            // process. Therefore, if this function fails for
                            // any reason, the returned F should be considered
                            // as undefined.
                            // 
                            // The rows of L are in the same order as the rows
                            // of A, while the columns of L are permuted such
                            // that L->v[j] (i.e., j-th column of L) contains
                            // the j-th pivot, which would be L->v[j]->x[0],
                            // (i.e., L->v[j]->i[0] == P[j]).
                            //
                            // U is stored as UT, i.e., U is in CSR format.
                            // The columns of U (or the rows of UT) are in the
                            // same order as the columns of A, while the rows
                            // of U (or the columns of UT) are permuted such
                            // that UT->v[j] (i.e., j-th column of UT, or j-th
                            // row of U) contains the j-th pivot, which would
                            // be UT->v[j]->x[0], (i.e., UT->v[j]->i[0] ==
                            // Q[j]). 
    SPEX_matrix *vk,        // Pointer to a n-by-1 dynamic_CSC matrix
                            // which contains the column to be inserted.
                            // The rows of vk are in the same order as A. The
                            // contained vector will be swapped with A->v[k] in
                            // the output upon return, regardless of failure.
    int64_t k,              // The column index that vk will be inserted, 0<=k<n
    const SPEX_options *option// Command parameters
);

//------------------------------------------------------------------------------
// Rank-1 Cholesky update/downdate
//------------------------------------------------------------------------------

// NOTE: the requirement for L is exactly same as other functions in
// SPEX_Update, which requires the pivot of L->v[j] be the
// first entry. In addition, A=L*D^(-1)L^T should be hold.

SPEX_info SPEX_Update_Chol_Rank1
(
    SPEX_factorization *F,  // The SPEX Cholesky factorization of A, including
                            // L, rhos, P and Pinv. This factorization will be
                            // modified during the update process. Therefore,
                            // if this function fails for any reason, the
                            // returned F should be considered as undefined.
                            // 
                            // The rows of L are in the same order as the rows
                            // of A, while the columns of L are permuted such
                            // that L->v[j] (i.e., j-th column of L) contains
                            // the j-th pivot, which would be L->v[j]->x[0],
                            // (i.e., L->v[j]->i[0] == P[j]).
    SPEX_matrix *w,         // a n-by-1 dynamic_CSC matrix that contains the
                            // vector to modify the original matrix A, the
                            // resulting A is A+sigma*w*w^T. In output, w is
                            // updated as the solution to L*D^(-1)*w_out = w
    const int64_t sigma,    // a nonzero scalar that determines whether
                            // this is an update or downdate
    const SPEX_options* option // Command options
);

//------------------------------------------------------------------------------
// Function for solving LD^(-1)Ux =b
//------------------------------------------------------------------------------

// NOTE: the requirement for L and UT is exactly same as other functions in
// SPEX_Update, which requires the pivot of L->v[j] or UT->v[j] be the
// first entry correspondingly.

SPEX_info SPEX_Update_LU_Solve // solves Ax = b via REF LU factorization of A
(
    // Output
    SPEX_matrix **x_handle, // a n*m dense matrix contains the solution to
                            // the system. 
    // input:
    SPEX_matrix *b,         // a n*m dense matrix contains the right-hand-side
                            // vector
    const SPEX_factorization *F,// The SPEX LU factorization in dynamic_CSC
                            // format.
    const SPEX_options* option // Command options
);

// SPEX_Update_LU_Solve_mod is similar to SPEX_UPdate_LU_Solve, which is used
// to solve LD^(-1)Ux=b. The only difference is that SPEX_Update_LU_Solve_mod
// will overwrite the solution to the right-hand-side matrix b. If any failure
// occures, each entry value of b should be considered as undefined.

SPEX_info SPEX_Update_LU_Solve_mod // solves LD^(-1)U x_out = x_in
(
    SPEX_matrix **x_handle, // Input as a n*m dense matrix contains the
                            // right-hand-side vectora and output as the
                            // solution to the system. 
    const SPEX_factorization *F,// The SPEX LU factorization in dynamic_CSC
                            // format.
    const SPEX_options* option // Command options
);

//------------------------------------------------------------------------------
// Function to convert the factorization between the updatable and non-updatable
// format. To obtain the updatable format, this function will obtain the
// transpose of U (for LU factorization), permute rows of L and UT, and make
// sure each column of the matrices have cooresponding pivot as the first
// entry. To otain the non-updatable format, this function will transpose UT
// (for LU factorization) and permute rows and L and U.
//------------------------------------------------------------------------------

SPEX_info SPEX_Update_factorization_convert
(
    SPEX_factorization **F_out,// The output factorization with same
                            // factorization kind as F_in
    const SPEX_factorization *F_in, // The factorization to be converted
    const bool updabtable, // true if wish to obtain updatable F.
    const SPEX_options* option // Command options
);

//------------------------------------------------------------------------------
// canonicalize a SPEX_DYNAMIC_CSC matrix such that each column of the input
// matrix have corresponding pivot as the first entry.
//------------------------------------------------------------------------------

// NOTE: This function is used to canonicalize L or UT before they can be used
// in any functions in the SPEX_Update library. perm can be NULL if there is
// no permutation.

SPEX_info SPEX_Update_matrix_canonicalize
(                       
    SPEX_matrix *A,         // the matrix to be canonicalize
    const int64_t *perm,    // the permuation vector applied on each vector of
                            // A, considered as identity if input as NULL
    const SPEX_options *option
);

//------------------------------------------------------------------------------
// permute the row indices of each column of a SPEX_DYNAMIC_CSC matrix such that
// A_out->v[j]->i[p] = perm[A_in->v[j]->i[p]].
//------------------------------------------------------------------------------

// NOTE: The L and U factorization from SPEX_Left_LU or the L from
// SPEX_Cholesky are all in SPEX_CSC format, and their columns and rows are
// permuted to be the same as the permuted matrix A(P,Q), and thus
// A(P,Q)=LD^(-1)U. However, all Update functions requires A=LD^(-1)U.
// Therefore, after performing SPEX_matrix_copy to get the L and/or UT
// SPEX_DYNAMIC_CSC form. Users need to perform this function to permute row
// indices of L or UT such that A=LD^(-1)U. And when all desired updates are
// performed and users wish to obtain L and/or U in the same format as those
// from SPEX_Left_LU or SPEX_Cholesky, this function should be called to
// perform the inverse of the permutation (before calling SPEX_matrix_copy to
// obtain L and/or U in the SPEX_CSC format).

SPEX_info SPEX_Update_permute_row
(
    SPEX_matrix *A,         // input matrix
    const int64_t *perm,    // desire permutation to be applied to A, must be
                            // non-NULL
    const SPEX_options *option
);

#endif

