//------------------------------------------------------------------------------
// SPEX_Util/Include/SPEX_Util.h: Include file for utility functions for SPEX
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#ifndef SPEX_UTIL_H
#define SPEX_UTIL_H

// SPEX_Util is a collection of utility functions for the SParse EXact package.
// Included are several routines for memory management, matrix operations, and 
// wrappers to the GMP library.
//
// This is the global include file and should be included in all SPEX_* packages
//
//
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Authors----------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Unless otherwise noted all Utility functions are authored by:
//
//    Christopher Lourenco, Jinhao Chen, Erick Moreno-Centeno, and Timothy Davis
//

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

//    SPEX is free software; you can redistribute it and/or modify
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

// SPEX is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------Include files required by SPEX --------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>
#include "SuiteSparse_config.h"

//------------------------------------------------------------------------------
// Version
//------------------------------------------------------------------------------

// Current version of the code
#define SPEX_UTIL_VERSION "1.1.0"
#define SPEX_UTIL_VERSION_MAJOR 1
#define SPEX_UTIL_VERSION_MINOR 1
#define SPEX_UTIL_VERSION_SUB   0

//------------------------------------------------------------------------------
// Error codes
//------------------------------------------------------------------------------

// Most SPEX functions return a code that indicates if it was successful
// or not. Otherwise the code returns a pointer to the object that was created
// or it returns void (in the case that an object was deleted)

typedef enum
{
    SPEX_OK = 0,                  // all is well
    SPEX_OUT_OF_MEMORY = -1,      // out of memory
    SPEX_SINGULAR = -2,           // the input matrix A is singular
    SPEX_INCORRECT_INPUT = -3,    // one or more input arguments are incorrect
    SPEX_INCORRECT = -4,          // The solution is incorrect
    SPEX_UNSYMMETRIC = -5,        // The input matrix is unsymmetric
    SPEX_NOTSPD = -6,             // The input matrix is not SPD
                                  // UNSYMMETRIC and NOTSPD are used for Cholesky
    SPEX_INCORRECT_ALGORITHM = -7,// The algorithm is not compatible with 
                                  // the factorization
    SPEX_PANIC = -8               // SPEX used without proper initialization
}
SPEX_info ;

//------------------------------------------------------------------------------
// Pivot scheme codes
//------------------------------------------------------------------------------

// A code in SPEX_options to tell SPEX what type of pivoting to use for pivoting
// in unsymmetric LU factorization.

typedef enum
{
    SPEX_SMALLEST = 0,      // Smallest pivot
    SPEX_DIAGONAL = 1,      // Diagonal pivoting
    SPEX_FIRST_NONZERO = 2, // First nonzero per column chosen as pivot
    SPEX_TOL_SMALLEST = 3,  // Diagonal pivoting with tol for smallest pivot.
                            //   (Default)
    SPEX_TOL_LARGEST = 4,   // Diagonal pivoting with tol. for largest pivot
    SPEX_LARGEST = 5        // Largest pivot
}
SPEX_pivot ;

//------------------------------------------------------------------------------
// Column ordering scheme codes
//------------------------------------------------------------------------------

// A code in SPEX_options to tell SPEX which column ordering to used prior to 
// exact factorization

typedef enum
{
    SPEX_NO_ORDERING = 0,   // None: A is factorized as-is
    SPEX_COLAMD = 1,        // COLAMD: Default
    SPEX_AMD = 2            // AMD
}
SPEX_col_order ;

//------------------------------------------------------------------------------
// Factorization type codes //TODO we need to change the word type (maybe "direction"?)
//------------------------------------------------------------------------------

// A code in SPEX_options to tell SPEX which "direction" to use when computing 
// the exact factorization

typedef enum
{
    SPEX_ALGORITHM_DEFAULT = 0,    // Defaults: Left for LU, Up for Chol, Gram for QR looking LU factorization //TODO maybe rename default...
    SPEX_LU_LEFT = 1,    // Left looking LU factorization
    SPEX_CHOL_LEFT = 2,  // Left looking Cholesky factorization
    SPEX_CHOL_UP = 3,    // Up looking Cholesky factorization
    SPEX_QR_GRAM = 4     // Default factorization for QR
}
SPEX_factorization_algorithm ; //TODO rename SPEX_factorization_algorithm and Propagate

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Data Structures--------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// This struct serves as a global struct to define all user-selectable options.

typedef struct SPEX_options
{
    SPEX_pivot pivot ;     // row pivoting scheme used.
    SPEX_col_order order ; // column ordering scheme used
    double tol ;           // tolerance for the row-pivotin methods
                           // SPEX_TOL_SMALLEST and SPEX_TOL_LARGEST
    int print_level ;      // 0: print nothing, 1: just errors,
                           // 2: terse (basic stats from COLAMD/AMD and
                           // SPEX Left LU), 3: all, with matrices and results
    int32_t prec ;         // Precision used to output file if MPFR is chosen
    mpfr_rnd_t round ;     // Type of MPFR rounding used
    bool check ;           // Set true if the solution to the system should be
                           // checked.  Intended for debugging only; SPEX is
                           // guaranteed to return the exact solution.
    SPEX_factorization_algorithm algo ; // parameter which tells the function
                           // whether it is performing a left-looking (2) or
                           // up-looking (3) factorization in case of Cholesky
                           // left looking LU (1) or if it is doing Gram-Schmidt
                           // (4). The Default (0) is up looking for Cholesky, 
                           // left looking for LU and Gram-Schmidt for QR
} SPEX_options ;

// Purpose: Create SPEX_options object with default parameters
// upon successful allocation, which are defined in SPEX_util_nternal.h
// To free it, simply use SPEX_FREE (*option).
SPEX_info SPEX_create_default_options (SPEX_options **option) ;

// check if SPEX_initialize* has been called
bool spex_initialized ( void ) ;        // true if called, false if not
void spex_set_initialized (bool s) ;    // set global initialzed flag to s

//------------------------------------------------------------------------------
// SPEX_matrix: a sparse CSC, sparse triplet, or dense matrix
//------------------------------------------------------------------------------

// SPEX uses a single matrix data type, SPEX_matrix, which can be held in
// one of three kinds of formats:  sparse CSC (compressed sparse column),
// sparse triplet, and dense:

typedef enum
{
    SPEX_CSC = 0,           // matrix is in compressed sparse column format
    SPEX_TRIPLET = 1,       // matrix is in sparse triplet format
    SPEX_DENSE = 2          // matrix is in dense format
}
SPEX_kind ;

// Each of the three formats can have values of 5 different data types: mpz_t,
// mpq_t, mpfr_t, int64_t, and double:

typedef enum
{
    SPEX_MPZ = 0,           // matrix of mpz_t integers
    SPEX_MPQ = 1,           // matrix of mpq_t rational numbers
    SPEX_MPFR = 2,          // matrix of mpfr_t
    SPEX_INT64 = 3,         // matrix of int64_t integers
    SPEX_FP64 = 4           // matrix of doubles
}
SPEX_type ;

// This gives a total of 15 different matrix types.  Not all functions accept
// all 15 matrices types, however.

// Suppose A is an m-by-n matrix with nz <= nzmax entries.
// The p, i, j, and x components are defined as:

// (0) SPEX_CSC:  A sparse matrix in CSC (compressed sparse column) format.
//      A->p is an int64_t array of size n+1, A->i is an int64_t array of size
//      nzmax (with nz <= nzmax), and A->x.type is an array of size nzmax of
//      matrix entries ('type' is one of mpz, mpq, mpfr, int64, or fp64).  The
//      row indices of column j appear in A->i [A->p [j] ... A->p [j+1]-1], and
//      the values appear in the same locations in A->x.type.  The A->j array
//      is NULL.  A->nz is ignored; nz is A->p [A->n].

// (1) SPEX_TRIPLET:  A sparse matrix in triplet format.  A->i and A->j are
//      both int64_t arrays of size nzmax, and A->x.type is an array of values
//      of the same size.  The kth tuple has row index A->i [k], column index
//      A->j [k], and value A->x.type [k], with 0 <= k < A->nz.  The A->p array
//      is NULL.

// (2) SPEX_DENSE:  A dense matrix.  The integer arrays A->p, A->i, and A->j
//      are all NULL.  A->x.type is a pointer to an array of size m*n, stored
//      in column-oriented format.  The value of A(i,j) is A->x.type [p]
//      with p = i + j*A->m.  A->nz is ignored; nz is A->m * A->n.

// The SPEX_matrix may contain 'shallow' components, A->p, A->i, A->j, and
// A->x.  For example, if A->p_shallow is true, then a non-NULL A->p is a
// pointer to a read-only array, and the A->p array is not freed by
// SPEX_matrix_free.  If A->p is NULL (for a triplet or dense matrix), then
// A->p_shallow has no effect.

typedef struct
{
    int64_t m ;         // number of rows
    int64_t n ;         // number of columns
    int64_t nzmax ;     // size of A->i, A->j, and A->x
    int64_t nz ;        // # nonzeros in a triplet matrix .
                        // Ignored for CSC and dense matrices.
    SPEX_kind kind ;    // CSC, triplet, or dense
    SPEX_type type ;    // mpz, mpq, mpfr, int64, or fp64 (double)

    int64_t *p ;        // if CSC: column pointers, an array size is n+1.
                        // if triplet or dense: A->p is NULL.
    bool p_shallow ;    // if true, A->p is shallow.

    int64_t *i ;        // if CSC or triplet: row indices, of size nzmax.
                        // if dense: A->i is NULL.
    bool i_shallow ;    // if true, A->i is shallow.

    int64_t *j ;        // if triplet: column indices, of size nzmax.
                        // if CSC or dense: A->j is NULL.
    bool j_shallow ;    // if true, A->j is shallow.

    union               // A->x.type has size nzmax.
    {
        mpz_t *mpz ;            // A->x.mpz
        mpq_t *mpq ;            // A->x.mpq
        mpfr_t *mpfr ;          // A->x.mpfr
        int64_t *int64 ;        // A->x.int64
        double *fp64 ;          // A->x.fp64
    } x ;
    bool x_shallow ;    // if true, A->x.type is shallow.

    mpq_t scale ;       // scale factor for mpz matrices (never shallow)
                        // For all matrices who's type is not mpz,
                        // mpz_scale = 1. 

} SPEX_matrix ;

//------------------------------------------------------------------------------
// SPEX_matrix_allocate: allocate an m-by-n SPEX_matrix
//------------------------------------------------------------------------------

// if shallow is false: All components (p,i,j,x) are allocated and set to zero,
//                      and then shallow flags are all false.

// if shallow is true:  All components (p,i,j,x) are NULL, and their shallow
//                      flags are all true.  The user can then set A->p,
//                      A->i, A->j, and/or A->x accordingly, from their own
//                      arrays.

SPEX_info SPEX_matrix_allocate
(
    SPEX_matrix **A_handle, // matrix to allocate
    SPEX_kind kind,         // CSC, triplet, or dense
    SPEX_type type,         // mpz, mpq, mpfr, int64, or double
    int64_t m,              // # of rows
    int64_t n,              // # of columns
    int64_t nzmax,          // max # of entries
    bool shallow,           // if true, matrix is shallow.  A->p, A->i, A->j,
                            // A->x are all returned as NULL and must be set
                            // by the caller.  All A->*_shallow are returned
                            // as true.
    bool init,              // If true, and the data types are mpz, mpq, or
                            // mpfr, the entries are initialized (using the
                            // appropriate SPEX_mp*_init function). If false,
                            // the mpz, mpq, and mpfr arrays are allocated but
                            // not initialized. Meaningless for data types 
                            // FP64 or INT64
    const SPEX_options *option
) ;

//------------------------------------------------------------------------------
// SPEX_matrix_free: free a SPEX_matrix
//------------------------------------------------------------------------------

SPEX_info SPEX_matrix_free
(
    SPEX_matrix **A_handle, // matrix to free
    const SPEX_options *option
) ;

//------------------------------------------------------------------------------
// SPEX_matrix_nnz: # of entries in a matrix
//------------------------------------------------------------------------------

SPEX_info SPEX_matrix_nnz     // find the # of entries in A
(
    int64_t *nnz,              // # of entries in A, -1 if A is NULL
    const SPEX_matrix *A,      // matrix to query
    const SPEX_options *option
) ;

//------------------------------------------------------------------------------
// SPEX_matrix_copy: makes a copy of a matrix
//------------------------------------------------------------------------------

// SPEX_matrix_copy: make a copy of a SPEX_matrix, into another kind and type.

SPEX_info SPEX_matrix_copy
(
    SPEX_matrix **C_handle, // matrix to create (never shallow)
    // inputs, not modified:
    SPEX_kind C_kind,       // C->kind: CSC, triplet, or dense
    SPEX_type C_type,       // C->type: mpz_t, mpq_t, mpfr_t, int64_t, or double
    SPEX_matrix *A,         // matrix to make a copy of (may be shallow)
    const SPEX_options *option
) ;

//------------------------------------------------------------------------------
// SPEX_matrix macros
//------------------------------------------------------------------------------

// These macros simplify the access to entries in a SPEX_matrix.
// The type parameter is one of: mpq, mpz, mpfr, int64, or fp64.

// To access the kth entry in a SPEX_matrix using 1D linear addressing,
// in any matrix kind (CSC, triplet, or dense), in any type:
#define SPEX_1D(A,k,type) ((A)->x.type [k])

// To access the (i,j)th entry in a 2D SPEX_matrix, in any type:
#define SPEX_2D(A,i,j,type) SPEX_1D (A, (i)+(j)*((A)->m), type)

//------------------------------------------------------------------------------
// SPEX_LU_analysis: symbolic pre-analysis
//------------------------------------------------------------------------------

// data structure for symbolic analysis
typedef enum
{
    SPEX_LU_SYMBOLIC_ANALYSIS = 0,                 // LU analysis
    SPEX_CHOLESKY_SYMBOLIC_ANALYSIS = 1,           // Cholesky analysis
    SPEX_QR_SYMBOLIC_ANALYSIS = 2                  // QR analysis
/* TODO or?
    SPEX_LU_ANALYSIS = 0,                 // LU analysis
    SPEX_CHOLESKY_ANALYSIS = 1,           // Cholesky analysis
    SPEX_QR_ANALYSIS = 2                  // QR analysis
*/
}SPEX_symbolic_analysis_kind ;

// This struct stores the results of symbolic analysis

typedef struct
{
    SPEX_symbolic_analysis_kind kind;             // LU, Cholesky or QR analysis

    //--------------------------------------------------------------------------
    // The permutations of the matrix that are found during the symbolic
    // analysis process.  One or more of these permutations could be NULL for
    // some SPEX_symbolic_analysis_kind. Specifically,
    // For kind == SPEX_LU_ANALYSIS, only Q_perm is not NULL.
    // For kind == SPEX_CHOLESKY_ANALYSIS, both Q_perm and Qinv_perm are
    // NULL.
    // TODO for QR????
    //--------------------------------------------------------------------------
    int64_t *P_perm;                     // row permutation
                                         // should not free by Left_LU_factorize
    int64_t *Pinv_perm;                  // inverse of row permutation

    int64_t *Q_perm;                     // column permutation
    int64_t *Qinv_perm;                  // inverse of column permutation
                                         // should not free by Left_LU_factorize

    //--------------------------------------------------------------------------
    // estimates of nonzeros that will apprear in the factorization
    //--------------------------------------------------------------------------

    int64_t lnz ;                        // Approximate number of nonzeros in L.
                                         // Available only for SPEX_LU_analysis
                                         // or SPEX_CHOLESKY_analysis.
    int64_t unz ;                        // Approximate number of nonzeros in U.
                                         // lnz and unz are used to allocate
                                         // the initial space for L and U; the
                                         // space is reallocated as needed.
                                         // Available only for SPEX_LU_analysis.

    //--------------------------------------------------------------------------
    // These are only used in the Cholesky analysis process
    //--------------------------------------------------------------------------
    int64_t* parent;                     // Elimination tree of target matrix
                                         // for Cholesky factorization.
    int64_t* cp;                         // column pointers of L for Cholesky
                                         // factorization.
} SPEX_symbolic_analysis ;

// The symbolic analysis object is created by SPEX_symbolic_analyze.

// SPEX_symbolic_analysis_free frees the SPEX_symbolic_analysis object.
SPEX_info SPEX_symbolic_analysis_free        
(
    SPEX_symbolic_analysis **S, // Structure to be deleted
    const SPEX_options *option
) ;

// SPEX_symbolic_analyze performs the symbolic ordering and analysis for
// specified kind. The pre-ordering method is specified in the option.
// For the pre-ordering for LU factorization, there are three options:
// no ordering, COLAMD, and AMD.
// TODO? For the pre-ordering for Cholesky factorization, 
SPEX_info SPEX_symbolic_analyze
(
    SPEX_symbolic_analysis **S_handle, // symbolic analysis
    const SPEX_matrix *A, // Input matrix
    const SPEX_symbolic_analysis_kind kind, // LU, Cholesky or QR analysis
    const SPEX_options *option  // Control parameters
) ;


//------------------------------------------------------------------------------
// SPEX_factorization: data structure for factorization
//------------------------------------------------------------------------------
// data structure for factorization
typedef enum
{
    SPEX_LU_FACTORIZATION = 0,            // LU factorization
    SPEX_CHOLESKY_FACTORIZATION = 1,      // Cholesky factorization
    SPEX_QR_FACTORIZATION = 2             // QR factorization
}SPEX_factorization_kind ;

typedef struct
{
    SPEX_factorization_kind kind;         // LU, Cholesky, QR factorization

    mpq_t scale_for_A;                    // the scale of the target matrix

    //--------------------------------------------------------------------------
    // These are used for LU or Cholesky factorization, but ignored for QR
    // factorization. Check L->kind to see if the factorization is updatable.
    //--------------------------------------------------------------------------

    SPEX_matrix *L;                       // The lower-triangular matrix from LU
                                          // or Cholesky factorization.
    SPEX_matrix *U;                       // The upper-triangular matrix from LU
                                          // factorization. NULL for Cholesky
                                          // factorization.
    SPEX_matrix *rhos;                    // A n-by-1 dense matrix for the
                                          // pivot values
    // This is only used in the Cholesky factorization
    int64_t* parent;                     // Elimination tree of target matrix
                                         // for Cholesky factorization.

    //--------------------------------------------------------------------------
    // These are used for QR factorization, but ignored for LU or Cholesky
    // factorization. Check L->kind to see if the factorization is updatable
    //--------------------------------------------------------------------------

    SPEX_matrix *Q;                       // The orthogonal matrix from QR
    SPEX_matrix *R;                       // The upper triangular matrix from QR

    //--------------------------------------------------------------------------
    // The permutations of the matrix that are used during the factorization.
    // These are currently used only for LU or Cholesky factorization.
    // One or more of these permutations could be NULL for some
    // SPEX_factorization_kind. Specifically,
    // For kind == SPEX_LU_FACTORIZATION, Qinv_perm will be NULL, but it will
    // be generated when the factorization is converted to the updatable form.
    // For kind == SPEX_CHOLESKY_FACTORIZATION, both Q_perm and Qinv_perm are
    // NULL.
    // TODO for QR????
    //--------------------------------------------------------------------------

    int64_t *P_perm;                     // row permutation
                                         // should not free by Left_LU_factorize
    int64_t *Pinv_perm;                  // inverse of row permutation

    int64_t *Q_perm;                     // column permutation
    int64_t *Qinv_perm;                  // inverse of column permutation
                                         // should not free by Left_LU_factorize

} SPEX_factorization;

// SPEX_factorization_free frees the SPEX_factorization object.
SPEX_info SPEX_factorization_free        
(
    SPEX_factorization **F, // Factorization to be deleted
    const SPEX_options *option
) ;


// SPEX_factorize performs SPEX factorization based on the available symbolic
// analysis S->kind, which should be either LU, Choleksy or QR.
SPEX_info SPEX_factorize
(
    SPEX_factorization **F_handle, // The resulted factorization as specified
    const SPEX_matrix *A, // Input matrix
    const SPEX_symbolic_analysis *S, // symbolic analysis
    const SPEX_options *option  // Control parameters
) ;

// SPEX_solve solves the linear system Ax=b with available factorization of A.
SPEX_info SPEX_solve
(
    SPEX_matrix **x_handle, // the solution to the system
    const SPEX_matrix *b,   // the right-hand-side vector
    const SPEX_factorization *F, // factorization of A
    const SPEX_options *option  // Control parameters
) ;


/* TODO
// SPEX_solve solves the linear system Ax=b with available factorization of A.
SPEX_info SPEX_solve
(
    SPEX_matrix **x, // input as the right-hand-side vector, and output with 
                     // the solution to the system
    const SPEX_factorization *F, // factorization of A
    const SPEX_options *option  // Control parameters
) ;
*/



//------------------------------------------------------------------------------
// Memory management
//------------------------------------------------------------------------------

// SPEX relies on the SuiteSparse memory management functions,
// SuiteSparse_malloc, SuiteSparse_calloc, SuiteSparse_realloc, and
// SuiteSparse_free.

// Allocate and initialize memory space for SPEX
void *SPEX_calloc
(
    size_t nitems,      // number of items to allocate
    size_t size         // size of each item
) ;

// Allocate memory space for SPEX
void *SPEX_malloc
(
    size_t size        // size of memory space to allocate
) ;

// Free the memory allocated by SPEX_calloc, SPEX_malloc, or SPEX_realloc.
void SPEX_free
(
    void *p         // pointer to memory space to free
) ;

// Free a pointer and set it to NULL.
#define SPEX_FREE(p)                        \
{                                           \
    SPEX_free (p) ;                         \
    (p) = NULL ;                            \
}

// SPEX_realloc is a wrapper for realloc.  If p is non-NULL on input, it points
// to a previously allocated object of size old_size * size_of_item.  The
// object is reallocated to be of size new_size * size_of_item.  If p is NULL
// on input, then a new object of that size is allocated.  On success, a
// pointer to the new object is returned.  If the reallocation fails, p is not
// modified, and a flag is returned to indicate that the reallocation failed.
// If the size decreases or remains the same, then the method always succeeds
// (ok is returned as true).

// Typical usage:  the following code fragment allocates an array of 10 int's,
// and then increases the size of the array to 20 int's.  If the SPEX_malloc
// succeeds but the SPEX_realloc fails, then the array remains unmodified,
// of size 10.
//
//      int *p ;
//      p = SPEX_malloc (10 * sizeof (int)) ;
//      if (p == NULL) { error here ... }
//      printf ("p points to an array of size 10 * sizeof (int)\n") ;
//      bool ok ;
//      p = SPEX_realloc (20, 10, sizeof (int), p, &ok) ;
//      if (ok) printf ("p has size 20 * sizeof (int)\n") ;
//      else printf ("realloc failed; p still has size 10 * sizeof (int)\n") ;
//      SPEX_free (p) ;

void *SPEX_realloc      // pointer to reallocated block, or original block
                        // if the realloc failed
(
    int64_t nitems_new,     // new number of items in the object
    int64_t nitems_old,     // old number of items in the object
    size_t size_of_item,    // sizeof each item
    void *p,                // old object to reallocate
    bool *ok                // true if success, false on failure
) ;

//------------------------------------------------------------------------------
// SPEX memory environment routines
//------------------------------------------------------------------------------

// SPEX_initialize: initializes the working evironment for SPEX library.
// It must be called prior to calling any other SPEX_* function.
SPEX_info SPEX_initialize (void) ;

// SPEX_initialize_expert is the same as SPEX_initialize, except that it allows
// for a redefinition of custom memory functions that are used for SPEX and
// GMP.  The four inputs to this function are pointers to four functions with
// the same signatures as the ANSI C malloc, calloc, realloc, and free.
SPEX_info SPEX_initialize_expert
(
    void* (*MyMalloc) (size_t),             // user-defined malloc
    void* (*MyCalloc) (size_t, size_t),     // user-defined calloc
    void* (*MyRealloc) (void *, size_t),    // user-defined realloc
    void  (*MyFree) (void *)                // user-defined free
) ;

// SPEX_finalize: This function finalizes the working evironment for SPEX
// library, and frees any internal workspace created by SPEX.  It must be
// called as the last SPEX_* function called.
SPEX_info SPEX_finalize (void) ;


// SPEX_matrix_check: check and print a SPEX_sparse matrix
SPEX_info SPEX_matrix_check     // returns a SPEX status code
(
    const SPEX_matrix *A,       // matrix to check
    const SPEX_options* option  // defines the print level
) ;


/* Purpose: This function takes as input a mpz_t SPEX_matrix and divides
 * it by an mpz_t constant storing the solution in a mpq_t dense SPEX_matrix
 * array. 
 */
SPEX_info SPEX_matrix_div // divides the x matrix by a scalar
(
    SPEX_matrix **x2_handle,    // x2 = x/scalar
    SPEX_matrix* x,             // input vector x
    const mpz_t scalar,         // the scalar
    const SPEX_options *option
) ;

/* Purpose: This function multiplies matrix x a scalar
 */
SPEX_info SPEX_matrix_mul   // multiplies x by a scalar
(
    SPEX_matrix *x,         // matrix to be multiplied
    const mpz_t scalar      // scalar to multiply by
) ;


/* SPEX_check_solution: checks the solution of the linear system.  Performs a
 * quick rational arithmetic check of A*x=b.
 */
SPEX_info SPEX_check_solution
(
    const SPEX_matrix *A,          // input matrix
    const SPEX_matrix *x,          // solution vector
    const SPEX_matrix *b,          // right hand side
    const SPEX_options* option           // Command options
);

/* Purpose: p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1]
 * into c.  This function is lightly modified from CSparse.
 */
SPEX_info SPEX_cumsum
(
    int64_t *p,          // vector to store the sum of c
    int64_t *c,          // vector which is summed
    int64_t n           // size of c
);


//------------------------------------------------------------------------------
//---------------------------SPEX GMP/MPFR Functions----------------------------
//------------------------------------------------------------------------------

// The following functions are the SPEX interface to the GMP/MPFR libary.
// Each corresponding GMP/MPFR function is given a wrapper to ensure that no
// memory leaks or crashes occur. All covered GMP functions can be found in
// SPEX_gmp.c

// The GMP library does not handle out-of-memory failures.  However, it does
// provide a mechanism for passing function pointers that replace GMP's use of
// malloc, realloc, and free.  This mechanism is used to provide a try/catch
// mechanism for memory allocation errors, using setjmp and longjmp.

// When a GMP function is called, this wrapper keeps track of a list of objects
// allocated by that function.  The list is started fresh each time a GMP
// function is called.  If any allocation fails, the NULL pointer is not
// returned to GMP.  Instead, all allocated blocks in the list are freed,
// and spex_gmp_allocate returns directly to wrapper.

SPEX_info SPEX_mpfr_asprintf (char **str, const char *format, ... ) ;

SPEX_info SPEX_gmp_fscanf (FILE *fp, const char *format, ... ) ;

SPEX_info SPEX_mpz_init (mpz_t x) ;

SPEX_info SPEX_mpz_init2(mpz_t x, const size_t size) ;

SPEX_info SPEX_mpz_set (mpz_t x, const mpz_t y) ;

SPEX_info SPEX_mpz_set_ui (mpz_t x, const uint64_t y) ;

SPEX_info SPEX_mpz_set_si (mpz_t x, const int64_t y) ;

SPEX_info SPEX_mpz_swap (mpz_t x, mpz_t y);

SPEX_info SPEX_mpz_get_d (double *x, const mpz_t y) ;

SPEX_info SPEX_mpz_get_si (int64_t *x, const mpz_t y) ;

SPEX_info SPEX_mpz_set_q (mpz_t x, const mpq_t y) ;

SPEX_info SPEX_mpz_mul (mpz_t a, const mpz_t b, const mpz_t c) ;

SPEX_info SPEX_mpz_mul_si (mpz_t a, const mpz_t b, const int64_t c) ;

SPEX_info SPEX_mpz_sub (mpz_t a, const mpz_t b, const mpz_t c) ;

SPEX_info SPEX_mpz_add (mpz_t a, const mpz_t b, const mpz_t c) ;

SPEX_info SPEX_mpz_addmul (mpz_t x, const mpz_t y, const mpz_t z) ;

SPEX_info SPEX_mpz_submul (mpz_t x, const mpz_t y, const mpz_t z) ;

SPEX_info SPEX_mpz_fdiv_q (mpz_t q, const mpz_t n, const mpz_t d) ;

SPEX_info SPEX_mpz_cdiv_q (mpz_t q, const mpz_t n, const mpz_t d) ;

SPEX_info SPEX_mpz_divexact (mpz_t x, const mpz_t y, const mpz_t z) ;

SPEX_info SPEX_mpz_gcd (mpz_t x, const mpz_t y, const mpz_t z) ;

SPEX_info SPEX_mpz_lcm (mpz_t lcm, const mpz_t x, const mpz_t y) ;

SPEX_info SPEX_mpz_neg (mpz_t x, const mpz_t y) ;

SPEX_info SPEX_mpz_abs (mpz_t x, const mpz_t y) ;

SPEX_info SPEX_mpz_cmp (int *r, const mpz_t x, const mpz_t y) ;

SPEX_info SPEX_mpz_cmpabs (int *r, const mpz_t x, const mpz_t y) ;

SPEX_info SPEX_mpz_cmp_ui (int *r, const mpz_t x, const uint64_t y) ;

SPEX_info SPEX_mpz_cmpabs_ui (int *r, const mpz_t x, const uint64_t y) ;

SPEX_info SPEX_mpz_sgn (int *sgn, const mpz_t x) ;

SPEX_info SPEX_mpz_sizeinbase (size_t *size, const mpz_t x, int64_t base) ;

SPEX_info SPEX_mpq_init (mpq_t x) ;

SPEX_info SPEX_mpq_set (mpq_t x, const mpq_t y) ;

SPEX_info SPEX_mpq_set_z (mpq_t x, const mpz_t y) ;

SPEX_info SPEX_mpq_canonicalize (mpq_t x);

SPEX_info SPEX_mpq_set_d (mpq_t x, const double y) ;

SPEX_info SPEX_mpq_swap (mpq_t x, mpq_t y);

SPEX_info SPEX_mpq_set_ui (mpq_t x, const uint64_t y, const uint64_t z) ;

SPEX_info SPEX_mpq_set_si (mpq_t x, const int64_t y, const uint64_t z) ;

SPEX_info SPEX_mpq_set_num (mpq_t x, const mpz_t y) ;

SPEX_info SPEX_mpq_set_den (mpq_t x, const mpz_t y) ;

SPEX_info SPEX_mpq_get_den (mpz_t x, const mpq_t y) ;

SPEX_info SPEX_mpq_get_d (double *x, const mpq_t y) ;

SPEX_info SPEX_mpq_neg (mpq_t x, const mpq_t y) ;

SPEX_info SPEX_mpq_abs (mpq_t x, const mpq_t y) ;

SPEX_info SPEX_mpq_add (mpq_t x, const mpq_t y, const mpq_t z) ;

SPEX_info SPEX_mpq_mul (mpq_t x, const mpq_t y, const mpq_t z) ;

SPEX_info SPEX_mpq_div (mpq_t x, const mpq_t y, const mpq_t z) ;

SPEX_info SPEX_mpq_cmp (int *r, const mpq_t x, const mpq_t y) ;

SPEX_info SPEX_mpq_cmp_ui (int *r, const mpq_t x,
                    const uint64_t num, const uint64_t den) ;

SPEX_info SPEX_mpq_cmp_z (int *r, const mpq_t x, const mpz_t y) ;

SPEX_info SPEX_mpq_sgn (int *sgn, const mpq_t x) ;

SPEX_info SPEX_mpq_equal (int *r, const mpq_t x, const mpq_t y) ;

SPEX_info SPEX_mpfr_init2(mpfr_t x, const uint64_t size) ;

SPEX_info SPEX_mpfr_set (mpfr_t x, const mpfr_t y, const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_set_d (mpfr_t x, const double y, const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_set_si (mpfr_t x, int64_t y, const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_set_q (mpfr_t x, const mpq_t y, const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_set_z (mpfr_t x, const mpz_t y, const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_get_z (mpz_t x, const mpfr_t y, const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_get_q (mpq_t x, const mpfr_t y, const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_get_d (double *x, const mpfr_t y, const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_get_si (int64_t *x, const mpfr_t y, const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_mul (mpfr_t x, const mpfr_t y, const mpfr_t z,
                    const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_mul_d (mpfr_t x, const mpfr_t y, const double z,
                    const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_div_d (mpfr_t x, const mpfr_t y, const double z,
                    const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_ui_pow_ui (mpfr_t x, const uint64_t y, const uint64_t z,
                    const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_sgn (int *sgn, const mpfr_t x) ;

SPEX_info SPEX_mpfr_free_cache (void) ;

SPEX_info SPEX_mpfr_free_str (char *str) ;

SPEX_info SPEX_gmp_printf (const char *format, ... ) ;
#if 0
// These functions are currently unused, but kept here for future reference.
SPEX_info SPEX_gmp_asprintf (char **str, const char *format, ... ) ;
SPEX_info SPEX_mpfr_printf ( const char *format, ... ) ;
SPEX_info SPEX_gmp_fprintf (FILE *fp, const char *format, ... ) ;
SPEX_info SPEX_mpfr_fprintf (FILE *fp, const char *format, ... ) ;
SPEX_info SPEX_mpz_set_d (mpz_t x, const double y) ;

SPEX_info SPEX_mpfr_log2(mpfr_t x, const mpfr_t y, const mpfr_rnd_t rnd) ;
#endif

/* WARNING: These functions have not been test covered!*/


/* Purpose: This function sets C = A' 
 * C_handle is NULL on input. On output, C_handle contains a pointer to A'
 */
SPEX_info SPEX_transpose
(
    SPEX_matrix **C_handle,     // C = A'
    SPEX_matrix *A              // Matrix to be transposed
);

SPEX_info SPEX_determine_symmetry
(
    SPEX_matrix* A,
    bool check_if_numerically_symmetric
            // if true, check A=A' (pattern & values). if false,
            // only check if the pattern of A is symmetric, not
            // the values
);

#endif

