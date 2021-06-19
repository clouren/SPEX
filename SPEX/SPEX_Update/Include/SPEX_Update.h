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
// the funtion fails for any reason, L and UT should be considered as
// undefined.

SPEX_info SPEX_Update_LU_ColRep
(
    SPEX_matrix *A,         // n-by-n original matrix in dynamic_CSC form
    SPEX_matrix *L,         // n-by-n row-permuted lower triangular
                            // factorization of A(P,Q) in dynamic_CSC form. The
                            // rows of L are in the same order as the rows of
                            // A, while the columns of L are permuted such that
                            // L->v[j] (i.e., j-th column of L) contains the
                            // j-th pivot, which would be L->v[j]->x[0], (i.e.,
                            // L->v[j]->i[0] == P[j]). This matrix will be
                            // modified during the update process. Therefore,
                            // if this function fails for any reason, the
                            // returned L should be considered as undefined.
    SPEX_matrix *UT,        // The transpose of U in dynamic_CSC form, where U
                            // is the n-by-n column-permuted upper triangular
                            // factorization of A(P,Q). The columns of U (or the
                            // rows of UT) are in the same order as the columns
                            // of A, while the rows of U (or the columns of UT)
                            // are permuted such that UT->v[j] (i.e., j-th
                            // column of UT, or j-th row of U) contains the
                            // j-th pivot, which would be UT->v[j]->x[0],
                            // (i.e., UT->v[j]->i[0] == Q[j]). This matrix will
                            // be modified during the update process.
                            // Therefore, if this function fails for any
                            // reason, the returned UT should be considered as
                            // undefined.
    SPEX_matrix *rhos,      // n-by-1 dense matrix for the array of pivots
    int64_t *P,             // Row permutation, P[i]-th row of A is used as the
                            // i-th pivot row
    int64_t *P_inv,         // Inverse of row permutation
    int64_t *Q,             // Column permutation, Q[j]-th column of A is used
                            // as the j-th pivot column
    int64_t *Q_inv,         // Inverse of column permutation
    SPEX_vector **vk,       // Pointer to the inserted column in the compressed
                            // column form, the rows of vk are in the same order
                            // as A. This vector will be swapped with A->v[k]
                            // in the output upon return, regardless of failure.
    int64_t k,              // The column index that vk will be inserted, 0<=k<n
    const SPEX_options *option// Command parameters
);

//------------------------------------------------------------------------------
// Rank-1 Cholesky update/downdate
//------------------------------------------------------------------------------

// NOTE: the requirement for L is exactly same as other functions in
// SPEX_Update, which requires the pivot of L->v[j] be the
// first entry. In addition, A=L*D^(-1)L^T should hold.

SPEX_info SPEX_Update_Chol_Rank1
(
    SPEX_matrix *L,         // n-by-n Cholesky factorization of A(P,P) in
                            // dynamic_CSC form. The rows of L are in the same
                            // order as the rows of A, while the columns of L
                            // are permuted such that L->v[j] (i.e., j-th
                            // column of L) contains the j-th pivot, which
                            // would be L->v[j]->x[0], (i.e., L->v[j]->i[0] ==
                            // P[j]). This matrix will be modified during the
                            // update process. Therefore, if this function
                            // fails for any reason, the returned L should be
                            // considered as undefined.
    SPEX_matrix *rhos,      // n-by-1 dense matrix for the array of pivots
    const int64_t *P,       // row permutation, P[i]-th row of A is used as the
                            // i-th pivot row
    const int64_t *P_inv,   // inverse of row permutation
    SPEX_vector *w,         // a n-by-1 vector in sparse compressed column form
                            // that modifies the original matrix A, the
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
//
// If users wish to overwrite the solution to the right-hand-side matrix b,
// the first input can be provided as &b, then each entry value of b should be
// considered as undefined if any failure occurs.

SPEX_info SPEX_Update_Solve // solves Ax = b via REF LU factorization of A
(
    // Output
    SPEX_matrix **x_handle, // a n*m dense matrix contains the solution to
                            // the system. If users wish to overwrite the
                            // solution to the right-hand-side matrix b, this
                            // can be provided as &b. Otherwise, new space will
                            // be allocated for x_handle
    // input:
    SPEX_matrix *b,         // a n*m dense matrix contains the right-hand-side                              // vector
    const SPEX_matrix *L,   // a n*n dynamic_CSC matrix that gives the lower
                            // triangular matrix
    const SPEX_matrix *UT,  // a n*n dynamic_CSC matrix that gives the transpose
                            // of the upper triangular matrix
    const mpq_t A_scale,    // scale of the input matrix
    int64_t *h,             // history vector
    const SPEX_matrix *rhos,// a n*1 dense matrix that gives the array of pivots
    const int64_t *P,       // row permutation
    const int64_t *Q_inv,   // inverse of column permutation
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

