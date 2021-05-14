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
#define SPEX_UPDATE_VERSION "1.0.0"
#define SPEX_UPDATE_VERSION_MAJOR 1
#define SPEX_UPDATE_VERSION_MINOR 0
#define SPEX_UPDATE_VERSION_SUB   0

//------------------------------------------------------------------------------
// SPEX_matrix: a sparse CSC, sparse triplet, or dense matrix
//------------------------------------------------------------------------------
#if 0
typedef struct
{
    int64_t nz;   // number of nonzeros. nz is meaningless for a dense vector
    int64_t nzmax;// size of array i and x, nz <= nzmax
    int64_t *i;   // array of size nzmax that contains the column/row indices
                  // of each nnz. For a dense vector, i == NULL
    mpz_t *x;     // array of size nzmax that contains the values of each nnz
    mpq_t scale;  // a scale factor that has not applied to entries in this v
} SPEX_vector;

SPEX_info SPEX_vector_alloc 
( 
    SPEX_vector **v_handle,         // vector to be allocated
    const int64_t nzmax,            // number of nnz entries in v
    const bool IsSparse              // indicate if the vector is sparse
);

SPEX_info SPEX_vector_realloc
(
    SPEX_vector* v, // the vector to be expanded
    const int64_t new_size// desired new size for v
);

void SPEX_vector_free
(
    SPEX_vector **v  // vector to be deleted
);

typedef struct
{
    int64_t m ;   // number of entries in each vector
    int64_t n ;   // number of vectors
    SPEX_vector **v;// array of size n, each entry of this array is a column/row
                  // vector.
    mpq_t scale ; // the scale factor applied to original entries to make
                  // them integers. That is, the original value can be obtained
                  // as v->x[*]/scale
} SPEX_mat;

SPEX_info SPEX_mat_alloc
(
    SPEX_mat **A_handle,      // matrix to be allocated
    const int64_t n,             // number of vectors
    const int64_t m,             // number of max entries in each vector
    const bool IsSparse          // indicate if the vector is sparse
);

void SPEX_mat_free
(
    SPEX_mat **A  // matrix to be deleted
);

SPEX_info SPEX_CSC_to_mat
(
    SPEX_mat **A_handle,          // converted SPEX_mat matrix
    const int64_t *perm,          // vector for permutation matrix, consider as
                                  // identity matrix if input as NULL
    const bool A_Is_ColWise,      // true if A->v[i] is the i-th col of A.
                                  // Otherwise, A->v[i] is the i-th row of A
    const SPEX_matrix *B,         // original matrix
    const SPEX_options *option
);

SPEX_info SPEX_mat_to_CSC
(
    SPEX_matrix **A_handle,       // converted CSC matrix
    const SPEX_mat *B,            // original matrix
    const int64_t *perm_inv,      // inverse of permutation, consider as
                                  // identity matrix if input as NULL
    const bool B_Is_ColWise,      // true if B->v[i] is the i-th col of B.
                                  // Otherwise, B->v[i] is the i-th row of B
    const SPEX_options *option
);

SPEX_info SPEX_mat_canonicalize
(
    SPEX_mat *A,    // the matrix to be canonicalize
    int64_t *perm   // the permuation vector applied on each vector of A,
                    // considered as identity if input as NULL
);

mpz_t* SPEX_create_mpz_array
(
    int64_t n            // size of the array
);

void SPEX_delete_mpz_array
(
    mpz_t **x,      // mpz array to be deleted
    int64_t n       // Size of x
);
#endif

//------------------------------------------------------------------------------
// Primary factorization update routine
//------------------------------------------------------------------------------

SPEX_info SPEX_Update_LU_ColRep
(
    SPEX_matrix *A,            // the original matrix in compressed-column form
    SPEX_matrix *L,            // stored in compressed-column form
    SPEX_matrix *U,            // stored in comptessed-row form
    SPEX_matrix *rhos,      // array of scaled pivots
    int64_t *P,             // row permutation
    int64_t *P_inv,         // inverse of row permutation
    int64_t *Q,             // column permutation
    int64_t *Q_inv,         // inverse of column permutation
    SPEX_vector **vk,       // pointer to the inserted column, which will be
                            // swapped with A->v[k] in the output if succeed
    int64_t k,              // the column index that vk will be inserted
    const SPEX_options *option// command parameters
);

SPEX_info SPEX_Update_Chol_Rank1
(
    SPEX_matrix *L,            // Lower-triangular factorization
    SPEX_matrix *rhos,      // array of scaled pivots
    int64_t *P,             // row permutation
    int64_t *P_inv,         // inverse of row permutation
    SPEX_vector *w,         // resulting matrix is A+sigma*w*w^T
    int64_t sigma
);

//------------------------------------------------------------------------------
// Function for solving LDUx =b
//------------------------------------------------------------------------------

SPEX_info SPEX_Update_Solve // solves Ax = b via REF LU factorization of A
(
    // Output
    SPEX_matrix **x_handle,    // solution to the system
    // input:
    SPEX_matrix *b,            // right hand side vector
    const SPEX_matrix *L,      // lower triangular matrix
    const SPEX_matrix *U,      // upper triangular matrix
    const mpq_t A_scale,    // scale of the input matrix
    int64_t *h,             // history vector
    const SPEX_matrix *rhos,// array of scaled pivots
    const int64_t *P,       // row permutation
    const int64_t *Q_inv,   // inverse of column permutation
    const SPEX_options* option // Command options
);


#endif

